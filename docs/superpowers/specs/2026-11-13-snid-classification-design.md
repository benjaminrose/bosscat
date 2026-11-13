# Design: SNID-SAGE classification of BOSS ToO spectra

Date: 2026-11-13. Status: draft for approval.

## Intent

Classify BOSS Target-of-Opportunity (ToO) transient spectra with SNID-SAGE
and write results to a table. SNID-SAGE is the first of several
classification tools (NGSF is already run for the current data; others will
follow). The design must make adding a tool cheap and keep every tool's
output traceable to its source.

## Current state (verified 2026-11-13)

- Working: FITS loader (`bosscat/io/boss.py`), CSV writer, pipeline, CLI,
  synthetic-FITS tests, CI. All pipeline tests mock the classifier.
- `bosscat/classify.py` calls SNID-SAGE 1.2.5 but does not work end to end:
  - `run_snid_analysis` defaults to `show_plots=True`; its plotting step
    raises `ImportError` here. With `show_plots=False` it runs (~11 s per
    spectrum) and returns a result.
  - `results.py` reads `result.age` / `result.age_error`, which do not exist
    (real: `consensus_age`, `consensus_age_error`).
  - `TEMPLATES_DIR` is hardcoded to a personal path outside the repo.
  - `ivar == 0` pixels are dropped before the call, leaving wavelength gaps.
- Data: 105 spectra in `data/spectra/` (read-only; do not modify). `data/manifest.csv`
  holds the NGSF results for the same 105 objects (`file_name` has no `.fits`).
  Only ~19 have a TNS type. Median `snr_median` is ~3.6 (range 1 to 20). 5 of 105
  have NGSF `z_fit` > 0.3 (max 0.67). BOSS `Z`/`CLASS` in the FITS files come
  from that NGSF run.

## Goals

1. `bosscat` runs SNID-SAGE on any BOSS coadd and writes one row per spectrum.
2. Production input has no redshift. SNID searches `-0.01 <= z <= 0.7`
   (upper limit set to cover the largest redshift in the NGSF manifest, 0.667).
3. Optional inverse-variance-weighted rebinning before SNID, with the factor
   configurable and its effect measured.
4. Each tool's output is namespaced and uses that tool's native field names.
5. Results are reproducible: configuration and versions are recorded.
6. Tests include at least one that runs real SNID-SAGE.
7. An evaluation script compares SNID output to the NGSF manifest and to
   TNS types where present.

## Non-goals

- Adding NGSF or any second classifier now (only keep the door open).
- Host-galaxy subtraction or extinction fitting (SNID does not do either).
- Changing `data/` or the FITS loader's contract.
- R code (ask first).

## Anti-goals (failures even if it "works")

- A SNID value appearing under a column named for another tool, or empty
  placeholder columns for fields SNID cannot produce (`A_v`, `host_template`,
  `frac_sn`, `second_delta_chi2`).
- The BOSS/NGSF redshift influencing SNID in any way.
- A confident-looking type on a spectrum that is out of range (z > `zmax`),
  quasar-like, or too noisy, with no flag.
- Mixing `redshift` and `consensus_redshift` or their errors.
- A passing test suite that never touches the real SNID-SAGE API.

## Constraints

- Python >= 3.10 (CI uses 3.14); `uv` for everything; NumPy-style docstrings.
- No new test dependency without asking. Minimal, scoped changes.
- Never hand-write BibTeX; use `ads-cite` if citations are added.
- No commit or push without an explicit request.

## Decisions (made with the user)

- **Blind classification.** The BOSS redshift is a test reference only; it is
  never a classifier input and is not required in production.
- **Redshift range:** `zmax` is raised to the largest redshift in the example
  data (NGSF `z_fit` max 0.667). The user said "the max in the example
  redshift"; the writer rounds to 0.7 (writer's call, configurable). `zmin`
  stays at SNID's default (-0.01).
- **Calibration set:** the 5 manifest objects with NGSF `z_fit` > 0.3 are the
  "high-z" calibration set. They define what a successful fit looks like
  (see Quality flags). Using the NGSF/BOSS redshift for this calibration, and
  for evaluation, is approved; it is never a classifier input.
- **Quality flags only annotate.** `snid_consensus_type` is always reported
  unchanged. No suppressed or "reliable type" column for now.
- **Rebinning may be dropped.** It ships only if it helps; the decision metric
  is the stability of SNID's results against the NGSF and TNS reference
  classifications across rebin factors.
- **Real-SNID tests are marked, skipped in CI, and run locally.** The SNID
  templates are not available to CI.
- **Rebinning** is wanted, to raise per-pixel S/N, because SNID templates are
  coarser than BOSS resolution.
- **Namespaced output.** `boss_*` for data/pipeline values, `snid_*` for
  SNID-SAGE, `ngsf_*` etc. for future tools. Native field names. Both
  `redshift` and `consensus_redshift` reported, each with its own error.
- **Classifier interface:** `f(wave, flux, ivar, meta) -> dict` of
  already-prefixed fields. The pipeline merges dicts.

## Design

### Components

- `bosscat/preprocess.py` (new): `rebin_loglam(wave, flux, ivar, factor)`.
  BOSS pixels are uniform in log10 wavelength, so an integer factor is exact.
  Combine with weights `ivar`, output `flux = sum(w*f)/sum(w)` and
  `ivar = sum(ivar)`. Bins with total `ivar == 0` are omitted. Adjacent-bin
  covariance is zero for integer binning, so error propagation is simple.
  `factor=1` is the identity.
- `bosscat/classifiers/snid.py` (replaces `classify.py`): wraps
  `preprocess_spectrum` and `run_snid_analysis` with `show_plots=False`,
  `zmax` from config, rebin factor from config. Returns the `snid_*` dict and
  quality flags.
- `bosscat/config.py` (new): templates directory (environment variable or
  config file; no hardcoded path), `zmax`, rebin factor, quality thresholds.
- `bosscat/pipeline.py`: loop over files and enabled classifiers, merge row
  dicts, catch errors per spectrum and per classifier (record `<tool>_error`),
  write CSV plus sidecar JSON.
- `bosscat/io/results.py`: shared-column row builder; `write_results` stays
  generic. The SNID flattening moves next to the SNID wrapper.
- `bosscat/evaluate.py` (new, separate from the pipeline): joins a results CSV
  to `manifest.csv` and reports redshift residuals, type agreement with NGSF
  `sn_family`, agreement with `tns_type` where present, and failure/flag rates.

### Output columns

- Shared: `filename`, `filepath`, `field`, `mjd`, `fiberid`, `ra`, `dec`,
  `too_id`, `programname`, `ebv`, `sn_median`, plus `boss_z`, `boss_z_err`,
  `boss_zwarning`, `boss_class`, `boss_subclass` (kept for evaluation only).
- SNID: `snid_success`, `snid_consensus_type`, `snid_best_subtype`,
  `snid_template_name`, `snid_template_type`, `snid_template_subtype`,
  `snid_redshift`, `snid_redshift_err`, `snid_consensus_redshift`,
  `snid_consensus_redshift_err`, `snid_consensus_age`,
  `snid_consensus_age_err`, `snid_lap`, `snid_r`, `snid_hsigma_lap_ccc`,
  `snid_frac_<type>` (from `type_fractions`), `snid_error`, and flags
  `snid_flag_*`.
- Sidecar JSON: package and SNID-SAGE versions, template set path and
  checksum if cheap, all config values, run time, input file list.

### Quality flags (SNID-specific)

Out-of-range is still a concern: a spectrum above `zmax` is forced into a wrong
match (observed with `zmax=1.0`: a BOSS-QSO at z=0.76 returned Ia, z=0.08).
Flags only annotate; they never change `snid_consensus_type`.

Success criterion for calibration (user-defined, with the 5 manifest objects
with NGSF `z_fit` > 0.3): with the extended `zmax`, SNID's `consensus_redshift`
should recover the NGSF `z_fit` for these objects, and the chosen flag
thresholds should separate the 5 from a wrong or low-quality fit. Exact
tolerance on redshift agreement is set by the builder from the data and
stated in the evaluation report. Five objects is a very small calibration
set: report the uncertainty and do not over-tune thresholds. Also check the
other ~100 objects for new false matches caused by the wider range (compare
`zmax` 0.3 versus 0.7). Flags proposed (thresholds calibrated, not hard-coded
guesses):

- `snid_flag_low_quality`: `lap`, `hsigma_lap_ccc`, or best-vs-second type
  margin below a threshold.
- `snid_flag_low_snr`: spectrum median S/N (after rebin) below a threshold.
- `snid_flag_nonsn_template`: best match is AGN/Galaxy/Star-type.
- `snid_flag_zedge`: consensus redshift within a small distance of `zmax`.
- Quasar-like spectra are only caught through `snid_flag_nonsn_template`
  (AGN/Galaxy/Star best match); no dedicated quasar flag for now.

Flags describe SNID's own output and data quality only. They do not use the
BOSS redshift.

## Testing

- Keep existing synthetic-FITS tests.
- Unit tests for `rebin_loglam`: identity at factor 1; flux conservation for
  constant spectra; `ivar` propagation against hand computation; masked-bin
  handling; non-divisible lengths.
- Tests for the SNID wrapper's result flattening use a stub with the **real**
  `SNIDResult` field names (taken from the installed class), so a rename in
  SNID-SAGE fails a test.
- Real-SNID tests carry a pytest marker (e.g. `snid_real`), are skipped in CI
  (templates unavailable there), and run locally. They also skip with a clear
  message if the templates directory is not configured. One such test runs
  real SNID-SAGE on a real spectrum and asserts the result has every field the
  flattener reads. Do not download anything into `data/`. CI keeps passing
  with only the stub-based tests; "done" for the real-SNID tests means "passes
  on the user's machine".
- Test that no BOSS redshift reaches the classifier call.
- Evaluation script tested on a tiny synthetic results/manifest pair.

## Risks and things to verify

- **Rebinning value is unproven.** SNID resamples to its own log grid, so the
  gain is mainly weighted pixel combination. It also smears narrow host lines.
  Measure classification stability against factor (e.g. 1, 2, 4, 8) on the 105
  spectra, judged by agreement with the NGSF and TNS classifications, before
  choosing a default. "Do not ship rebinning" is an acceptable outcome.
  Caveat: only ~19 spectra have a TNS type, so most of the evidence is NGSF
  agreement, which is not truth. Compare SNID's own grid spacing with
  BOSS's ~69 km/s pixels first.
- **Wavelength system.** BOSS wavelengths are vacuum. Check whether SNID-SAGE
  templates are air or vacuum (unverified); if they differ, conversion (about
  1 A at 5000 A) affects redshift at the ~1e-4 level.
- **Galactic extinction.** `ebv` is available but not applied; decide whether
  to deredden (builder: raise with the user if it changes results noticeably).
- **Host contamination** and quasar-like objects can drive confident wrong
  types; NGSF handles hosts, SNID does not. Do not compare SNID to NGSF as if
  they were measuring the same thing.
- **Speed.** ~11 s per spectrum; 105 spectra is ~20 minutes serial. Consider
  per-spectrum parallelism only if needed.
- **Comparison caveat.** The manifest is NGSF output, not truth. Only ~19
  objects have a TNS type.

## Calls made by the writer (not yet agreed with the user)

- Rebinning is a shared helper in `preprocess.py`; each classifier decides
  whether and how much to use it through its config (SNID: factor from
  config).
- Config via environment variable for the templates path plus a small
  config file for numeric settings.
- Sidecar JSON for provenance rather than extra CSV columns.
- Slow real-SNID test is opt-in via a pytest marker.
- The evaluation script lives in the package (`bosscat/evaluate.py`), not
  in `tests/`.

- Air/vacuum conversion and dereddening are NOT applied by default. The builder
  measures their effect and reports it; enabling either needs user sign-off.
- Masked (`ivar == 0`) pixels: the builder tests, on real spectra, whether
  SNID-SAGE copes with gaps (passing only valid pixels) versus filling them,
  and documents what the wrapper does. Not assumed.
- Shared column `sn_median` is the BOSS `SN_MEDIAN_ALL` from the loader, not the
  manifest's `snr_median` (different quantity, possibly different definition).
  The evaluation script joins to the manifest on file name (manifest
  `file_name` + `.fits`) and reports both columns.
- No known consumers of the old `classify.py` import path or CSV columns
  (repo is new, 8 commits). Replace outright. Existing CSV columns keep their
  names except the SNID age fields, which are corrected.
- All writer calls above are planned as written unless the user rejects them.

## Left to the builder

Exact flag thresholds (calibrate on data), config file format, rebin default,
CLI flag spelling (`--classifier snid`), CI strategy for templates, and
whether to parallelise.
