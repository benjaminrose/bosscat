# bosscat: project notes for Claude

BOSS Classification Analysis Toolkit. Python package `bosscat` (src layout not used; package lives at `bosscat/`). Research code for supernova/spectra classification.

## Goals and context

- Classify BOSS ToO (transient) coadd spectra (`ToO-specLite-FIELD-MJD-CATALOGID.fits`) with multiple tools. SNID-SAGE is first; NGSF and others are expected later.
- `data/manifest.csv` (105 rows) is the existing **NGSF** output for the 105 spectra in `data/spectra/` (columns `name, file_name, snr_median, tns_type, sn_template, sn_family, z_fit, ...`). Use it as a comparison, not ground truth. Only ~19 rows have a TNS type. `file_name` omits `.fits`.
- BOSS `Z`/`CLASS` in SPALL are the BOSS pipeline's own values (QSO/GALAXY/STAR), NOT NGSF's: only 45 of 105 agree with the manifest `z_fit` within 0.01. The NGSF reference is the manifest `z_fit`/`sn_family`; `boss_*` columns are context only. Never pass any redshift into a classifier (`forced_redshift`); production inputs will not have one.
- Typical data are low S/N (manifest `snr_median` median ~3.6, range 1-20). Redshifts: 5 of 105 have NGSF z > 0.3 (max 0.67).
- Calibration results: `docs/superpowers/reports/2026-11-13-snid-calibration.md`. Run outputs go in `runs/` (gitignored, never committed).
- Design doc: `docs/superpowers/specs/2026-11-13-snid-classification-design.md`.

## Design rules (decided with the user)

- Output columns are namespaced by their source: `boss_*` for BOSS pipeline values, `snid_*` for SNID-SAGE, `ngsf_*` for any future NGSF. Never put a SNID value in a column named for another tool, and never leave placeholder columns for fields a tool cannot produce (e.g. SNID has no `A_v`, `host_template`, `frac_sn`).
- Use each tool's native field names and report both of its values when it has two (e.g. SNID `redshift` and `consensus_redshift`), each with its own error.
- Each classifier takes `(wave, flux, ivar, meta)` and returns a dict of already-prefixed fields; the pipeline merges dicts. Record each tool's configuration (version, templates, z range, rebin factor) in a sidecar `<output file name>.provenance.json` (`complete`, `n_rows`). On Ctrl-C or a crash the rows so far are written and `complete` is false; output names ending `.provenance.json` are rejected.
- Blind rule is enforced at the interface: the pipeline strips `z`, `z_err`, `class`, `subclass`, `zwarning` from the `meta` passed to classifiers (the `boss_*` columns still come from the full meta).
- SNID redshift search `zmax` is 0.7 (covers the manifest max, 0.667); the 5 manifest objects with z > 0.3 were the calibration set. Spectra that may lie beyond `zmax` or fit poorly get `snid_flag_*` annotations. Flags only annotate: never suppress or alter `snid_consensus_type`.
- Rebinning ships only if it improves stability against NGSF/TNS. Default is `rebin_factor = 1`; factor 4/8 (higher NGSF/TNS agreement) is pending user decision, see the calibration report.
- Air/vacuum conversion and Galactic dereddening are measured but not enabled; enabling either needs user sign-off.
- Using the NGSF/BOSS redshift for calibration and evaluation is approved; using it as a classifier input is not.

## SNID-SAGE facts (verified by running v1.2.5, 2026-11-13)

- `run_snid_analysis` defaults to `show_plots=True`, and its plotting step raises `ImportError` in this install. Always pass `show_plots=False`; `_snid_api` defaults `MPLBACKEND=Agg`. ~11 s per spectrum serial; a full 105-spectrum run is ~11 min with 3 runs in parallel.
- `redshift` (best template) differs from `consensus_redshift`. `consensus_age`, `consensus_age_error` and `r` are always 0.0 in v1.2.5 (placeholders), so `snid_consensus_age`, `snid_consensus_age_err` and `snid_r` are NOT emitted. `template_age` is populated; emitting it as `snid_template_age` is a candidate that needs the user's sign-off.
- SNID-SAGE already discards matches with lap < 0.3 or HσLAP-CCC < 1.5 (its defaults) and then returns `success=False`, so `snid_flag_low_quality` at those thresholds marks only failed fits. All SNID failures here have S/N < 3.
- SNID resamples to 1024 log bins over 2500-10000 Å (406 km/s per bin, ~6 BOSS pixels) and takes flux only, no ivar. Template air/vacuum convention is undocumented.
- Templates: `/Users/ben_rose/Documents/1 Projects/as5-classifier/SNID-SAGE/templates` (read-only), passed via `--templates-dir` or `BOSSCAT_SNID_TEMPLATES`.

## Environment

- Python >=3.10 (see `pyproject.toml`). Local venv is `.venv/` (gitignored).
- Dependencies are declared in `pyproject.toml`; the lock file is `uv.lock`. Use `uv` for installs and runs, not bare `pip`.
- Entry point: `bosscat` -> `bosscat.__main__:main`.

## Data

- `data/` holds the spectra and `manifest.csv`. It is gitignored and large.
- Never modify, move, delete, or commit anything in `data/` without asking.
- Read data through the loaders in `bosscat/io/` rather than ad hoc paths.

## Code conventions

- Python is the primary language. R is for learning; ask before adding R code to the package.
- Use NumPy-style docstrings.
- Keep changes minimal and scoped to the request; do not refactor unrelated modules.

## Layout

- `bosscat/io/boss.py`: `read_boss_spectrum` (COADD + SPALL -> wave, flux, ivar, meta).
- `bosscat/io/results.py`: `base_row`, `write_results` (CSV), `write_provenance`/`sidecar_path` (sidecar JSON).
- `bosscat/config.py`: `SnidConfig` (templates dir, z range, rebin factor, flag thresholds).
- `bosscat/preprocess.py`: `rebin_loglam`.
- `bosscat/classifiers/snid.py`: `make_snid`, `snid_result_to_fields`, `snid_flags` (wraps `preprocess_spectrum` then `run_snid_analysis`).
- `bosscat/pipeline.py`: `run_pipeline`; `bosscat/__main__.py`: CLI; `bosscat/evaluate.py`: `python -m bosscat.evaluate RESULTS.csv` vs the manifest (NGSF `super_chandra*` counts as Ia, `TDE*` as TDE; other unmapped reference labels are listed, not silently dropped).

## Testing

- pytest + pytest-cov (dev dependency group); run with `uv run pytest`, add `--cov` for coverage (branch coverage, config in `pyproject.toml`). CI (`.github/workflows/tests.yml`) runs `uv run pytest --cov --cov-fail-under=95`; the threshold lives in the workflow, not `pyproject.toml`. Coverage was 96.4% (CI-like and with real SNID alike) when set; raise the threshold when coverage rises, never lower it to make CI pass.
- `tests/conftest.py` builds synthetic BOSS-style FITS files. SNID-SAGE is stubbed in most tests; `tests/test_snid_real.py` runs real SNID-SAGE (marker `snid_real`, skipped without templates and in CI).
- Ask before adding further test dependencies.

## Bibliography

- Never hand-write BibTeX. Every `.bib` entry comes verbatim from NASA ADS via `/ads-cite` or the `ads-cite` CLI, using `--rekey` (citekeys like `Surname_Subject_Year`).
- Use `natbib` or `biblatex` with `.bib` files; do not hardcode citations in `.tex`.

## Git

- SNID work is on branch `feature/snid-classification`. Untracked: `.claude/`, `tns_2025ahsk_output_example.ascii`, `runs/` (ignored).
- Do not commit or push without an explicit request.
