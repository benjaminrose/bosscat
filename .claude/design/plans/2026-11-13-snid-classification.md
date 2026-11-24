# SNID-SAGE Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run SNID-SAGE blind on BOSS ToO coadd spectra and write namespaced, reproducible results, with an evaluation against the NGSF manifest.

**Architecture:** Each classifier is a `Classifier` object (name, `run(wave, flux, ivar, meta) -> dict` of already-prefixed fields, `provenance() -> dict`). The pipeline reads a spectrum, builds the shared `boss_*` row, merges each classifier's dict, and writes a CSV plus a sidecar JSON. Optional inverse-variance rebinning and quality flags live in/around the SNID wrapper; a separate evaluation module compares results to `data/manifest.csv`.

**Tech Stack:** Python >= 3.10 (CI 3.14), `uv`, numpy, astropy, snid-sage 1.2.5, pytest.

**Spec:** `.claude/design/specs/2026-11-13-snid-classification-design.md` (read it first; also read `CLAUDE.md`).

## Global Constraints

- Production classification is blind: no BOSS/NGSF redshift ever reaches a classifier (never pass `forced_redshift`; never read `meta["z"]` inside a classifier). `boss_*` columns are for evaluation only.
- SNID search range: `zmin = -0.01`, `zmax = 0.7` (writer's rounding of the manifest max 0.667; configurable). Always call `run_snid_analysis(..., show_plots=False)`.
- Output columns are namespaced by source (`boss_*`, `snid_*`, future `ngsf_*`) using the tool's native field names. No placeholder columns for fields SNID cannot produce (`A_v`, `host_template`, `frac_sn`, `second_delta_chi2`). Report both SNID `redshift` and `consensus_redshift`, each with its own error. Age comes from `consensus_age`/`consensus_age_error`.
- Quality flags only annotate; `snid_consensus_type` is never changed or blanked.
- No air/vacuum conversion and no dereddening by default.
- Never modify, move, delete, or commit anything in `data/`. Do not commit CSV/JSON run outputs (put them in `runs/`, gitignored).
- Real-SNID tests are marked `snid_real`, skipped in CI, and run locally.
- NumPy-style docstrings; minimal scoped changes; no new dependencies (ask first); no R code.
- Commit steps below only happen when the user has asked for commits (CLAUDE.md). Otherwise leave changes uncommitted and say so.

## Review Focus

Inputs the spec implies but no happy-path test would catch, most likely first:

1. A spectrum with all `ivar == 0` (or too few valid pixels): row gets `snid_success=False` and `snid_error`, batch continues (Task 3).
2. NaN/inf in `flux` where `ivar > 0`: treated as masked, never passed to SNID (Tasks 2, 3).
3. Spectrum length not a multiple of the rebin factor, or too short after rebinning for SNID's 2000 A overlap requirement: error recorded, no crash (Tasks 2, 3).
4. SNID result with `success=False` or `None` fields (e.g. `consensus_age`): CSV blanks, flags set True for missing evidence (Tasks 3, 4).
5. Two classifiers where one raises: the other's columns survive and the CSV column union stays consistent (Task 5).

---

## File Structure

- Create `bosscat/config.py`: `SnidConfig`.
- Create `bosscat/preprocess.py`: `rebin_loglam`.
- Create `bosscat/classifiers/__init__.py`: `Classifier` dataclass.
- Create `bosscat/classifiers/snid.py`: SNID wrapper, flattening, flags, `make_snid`.
- Modify `bosscat/io/results.py`: `base_row`, `write_results`, `write_provenance` (remove `snid_result_to_row`).
- Modify `bosscat/pipeline.py`, `bosscat/__main__.py`; delete `bosscat/classify.py`.
- Create `bosscat/evaluate.py`.
- Modify `pyproject.toml` (pytest marker), `tests/conftest.py` (skip hook), `.gitignore` (`runs/`), `README.md`, `CLAUDE.md`.
- Tests: `tests/test_config.py`, `tests/test_preprocess.py`, `tests/test_snid_classifier.py`, `tests/test_snid_flags.py`, `tests/test_evaluate.py`, `tests/test_snid_real.py`; update `tests/test_cli.py`.

### Task 1: SnidConfig

**Files:** Create `bosscat/config.py`; Test `tests/test_config.py`

**Interfaces:**
- Produces: frozen dataclass `SnidConfig(templates_dir: str | None = None, zmin: float = -0.01, zmax: float = 0.7, rebin_factor: int = 1, min_lap: float = 0.3, min_hsigma_lap_ccc: float = 1.5, min_snr: float = 3.0, zedge_margin: float = 0.02)` with `resolved_templates_dir(self) -> str` and `as_dict(self) -> dict`. `min_lap`/`min_hsigma_lap_ccc` default to SNID-SAGE's own defaults; `min_snr` and `zedge_margin` are provisional and replaced in Task 8. Env var: `BOSSCAT_SNID_TEMPLATES`.

- [ ] **Step 1: Write failing tests** in `tests/test_config.py`: `test_templates_from_env` (monkeypatch env to `tmp_path`, `resolved_templates_dir() == str(tmp_path)`); `test_explicit_overrides_env`; `test_missing_templates_raises` (env unset, `templates_dir=None` -> `RuntimeError` whose message contains `BOSSCAT_SNID_TEMPLATES`); `test_nonexistent_dir_raises`; `test_defaults` (`zmax == 0.7`, `zmin == -0.01`, `rebin_factor == 1`).
- [ ] **Step 2:** Run `uv run pytest tests/test_config.py -v`; expect FAIL (import error).
- [ ] **Step 3:** Implement `SnidConfig` as specified; `as_dict` uses `dataclasses.asdict`.
- [ ] **Step 4:** Run the tests; expect PASS.
- [ ] **Step 5:** Commit `bosscat/config.py tests/test_config.py` (if the user asked for commits).

### Task 2: Inverse-variance log-lambda rebinning

**Files:** Create `bosscat/preprocess.py`; Test `tests/test_preprocess.py`

**Interfaces:**
- Produces: `rebin_loglam(wave, flux, ivar, factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]` returning `(wave, flux, ivar)`.

Behavior: pixels with `ivar <= 0` or non-finite `flux` or `ivar` are masked (weight 0). Integer `factor` groups consecutive pixels; trailing pixels that do not fill a bin are dropped. Per bin: `flux = sum(ivar*flux)/sum(ivar)`, `ivar = sum(ivar)`, `wave = 10**mean(log10(wave))`. Bins with `sum(ivar) == 0` are omitted. `factor == 1` returns only the valid pixels (same masking rule). Raises `ValueError` if `factor < 1` or `log10(wave)` spacing is not uniform within 1e-3 relative tolerance.

- [ ] **Step 1: Write failing tests**: `test_factor_one_keeps_valid_pixels` (masked pixel removed, others equal); `test_hand_computed_bin` (`flux=[1,3]`, `ivar=[1,3]`, factor 2 -> flux `2.5`, ivar `4.0`); `test_constant_flux_preserved`; `test_fully_masked_bin_omitted` (`ivar=[0,0,2,2]`, factor 2 -> length 1); `test_partially_masked_bin_uses_valid_pixel`; `test_remainder_dropped` (length 5, factor 2 -> 2 bins); `test_nan_flux_with_positive_ivar_is_masked`; `test_non_log_uniform_raises`; `test_bad_factor_raises`.
- [ ] **Step 2:** Run `uv run pytest tests/test_preprocess.py -v`; expect FAIL.
- [ ] **Step 3:** Implement with reshape-and-sum on the truncated arrays.
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit (if asked).

### Task 3: SNID wrapper and flattening

**Files:** Create `bosscat/classifiers/__init__.py`, `bosscat/classifiers/snid.py`; Test `tests/test_snid_classifier.py`; Delete `bosscat/classify.py` in Task 5 (leave it until the pipeline stops importing it).

**Interfaces:**
- Consumes: `SnidConfig` (Task 1), `rebin_loglam` (Task 2).
- Produces:
  - `Classifier(name: str, run: Callable[[ndarray, ndarray, ndarray, dict], dict], provenance: Callable[[], dict])` frozen dataclass in `bosscat/classifiers/__init__.py`.
  - `_snid_api() -> tuple[Callable, Callable]` returning `(preprocess_spectrum, run_snid_analysis)` via lazy import (so importing `bosscat` does not need SNID-SAGE).
  - `snid_result_to_fields(result) -> dict` with keys: `snid_success`, `snid_consensus_type`, `snid_best_subtype`, `snid_template_name`, `snid_template_type`, `snid_template_subtype`, `snid_redshift`, `snid_redshift_err` (from `redshift_error`), `snid_consensus_redshift`, `snid_consensus_redshift_err`, `snid_consensus_age`, `snid_consensus_age_err`, `snid_lap`, `snid_r`, `snid_hsigma_lap_ccc`, and `snid_frac_<type>` per `type_fractions` (`/` and spaces replaced by `_`).
  - `make_snid(cfg: SnidConfig) -> Classifier` with `name == "snid"`. `run` rebins (`cfg.rebin_factor`), calls `preprocess_spectrum(input_spectrum=(wave, flux))` then `run_snid_analysis(processed, templates_dir, zmin=cfg.zmin, zmax=cfg.zmax, show_plots=False)`, merges `snid_result_to_fields` with `snid_snr_median` and the flags from Task 4. On no valid pixels returns `{"snid_success": False, "snid_error": "no valid pixels"}` without calling SNID. `provenance()` returns `{"snid_sage_version": ..., "config": cfg.as_dict(), "templates_dir": ...}`.

- [ ] **Step 1: Write failing tests** (monkeypatch `bosscat.classifiers.snid._snid_api` with recorders; build stub results as `SimpleNamespace(**{f.name: None for f in dataclasses.fields(SNIDResult)})` so any attribute not on the real class fails): `test_flattener_uses_consensus_age_and_both_redshifts`; `test_type_fraction_keys_sanitized` (`"Ia-norm/91T"` -> `snid_frac_Ia-norm_91T`; `"II P"` -> `snid_frac_II_P`); `test_calls_snid_blind_with_plots_off` (kwargs contain `show_plots=False`, `zmax=0.7`, `zmin=-0.01`, no `forced_redshift`, even when `meta={"z": 0.5}`); `test_only_valid_pixels_reach_snid` (NaN flux with `ivar>0` and `ivar==0` pixels excluded); `test_rebin_factor_applied` (factor 2 halves array length); `test_all_masked_returns_error_without_calling_snid`; `test_short_spectrum_error_propagates` (the stub SNID raises `ValueError`; `run` lets it propagate, since the pipeline records it in Task 5); `test_failed_result_keeps_columns` (`success=False`, `consensus_age=None` -> `snid_success False`, value `None`, no exception).
- [ ] **Step 2:** Run `uv run pytest tests/test_snid_classifier.py -v`; expect FAIL.
- [ ] **Step 3:** Implement per Interfaces. Compute `snid_snr_median` as the median of `flux*sqrt(ivar)` over valid rebinned pixels.
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit (if asked).

### Task 4: SNID quality flags (annotate only)

**Files:** Modify `bosscat/classifiers/snid.py`; Test `tests/test_snid_flags.py`

**Interfaces:**
- Consumes: `SnidConfig` thresholds; result stub from Task 3.
- Produces: `snid_flags(result, snr_median: float, cfg: SnidConfig) -> dict[str, bool]` with keys `snid_flag_low_quality`, `snid_flag_low_snr`, `snid_flag_nonsn_template`, `snid_flag_zedge`, merged into the `run` output of Task 3.

Rules: `low_quality` if `lap < min_lap` or `hsigma_lap_ccc < min_hsigma_lap_ccc` (a `None` value counts as flagged); `low_snr` if `snr_median < min_snr`; `nonsn_template` if `consensus_type` or `template_type` is in `{"AGN", "Galaxy", "Star", "CV"}` (verify type names against the template index in the installed SNID-SAGE and fix the set if they differ); `zedge` if `consensus_redshift >= cfg.zmax - cfg.zedge_margin` (`None` counts as flagged). Flags never alter any other key.

- [ ] **Step 1: Write failing tests**: one per flag for true/false cases; `test_none_fields_flag_true`; `test_flags_do_not_change_type` (type key identical with and without flags firing).
- [ ] **Step 2:** Run `uv run pytest tests/test_snid_flags.py -v`; expect FAIL.
- [ ] **Step 3:** Implement `snid_flags` and call it from `make_snid`'s `run`.
- [ ] **Step 4:** Run the whole suite `uv run pytest -v`; expect PASS.
- [ ] **Step 5:** Commit (if asked).

### Task 5: Pipeline, results, CLI

**Files:** Modify `bosscat/io/results.py`, `bosscat/pipeline.py`, `bosscat/__main__.py`, `tests/test_cli.py`, `README.md`; Delete `bosscat/classify.py`

**Interfaces:**
- Consumes: `Classifier`, `make_snid`, `SnidConfig`.
- Produces:
  - `base_row(meta: dict) -> dict` (the shared and `boss_*` columns of the old `snid_result_to_row`, no SNID columns).
  - `write_provenance(output_path, classifiers: Sequence[Classifier], input_files: Sequence[str]) -> Path` writing `<output stem>.json` next to the CSV with `bosscat_version`, `python_version`, UTC `created`, `input_files`, and `classifiers: {name: provenance()}`.
  - `run_pipeline(input_files, output_path, classifiers: Sequence[Classifier] = ()) -> None`. Per classifier exceptions are caught and recorded as `<name>_error` and `<name>_success=False`; read errors keep the existing behavior (row with `filepath` and `error`).
  - CLI: `--classifier [NAME ...]` (default `snid`; given with no names runs none), `--zmax` (float), `--rebin-factor` (int), `--templates-dir`, `-o/--output`. Replaces `--snid/--no-snid`.

- [ ] **Step 1: Write/update failing tests** in `tests/test_cli.py` and a new `tests/test_pipeline.py`: replace the `--no-snid` tests with `--classifier` (no names) and assert no `snid_*` columns; `test_classifier_columns_merged` (fake `Classifier` returning `{"fake_a": 1}`); `test_one_classifier_fails_other_survives` (two fakes, one raises; the other's column and the failing one's `_error` both present, CSV columns consistent); `test_read_error_row_kept`; `test_sidecar_json_written` (keys above, classifier provenance included); `test_boss_redshift_not_passed_to_classifier` (fake records its `meta` argument and the pipeline-provided `wave, flux, ivar`; assert the classifier is called with 4 arguments and that production code never calls `forced_redshift`).
- [ ] **Step 2:** Run `uv run pytest -v`; expect FAIL.
- [ ] **Step 3:** Implement; delete `bosscat/classify.py`; update the README usage section (CLI flags, `BOSSCAT_SNID_TEMPLATES`, output columns, sidecar).
- [ ] **Step 4:** Run `uv run pytest -v`; expect all PASS.
- [ ] **Step 5:** Commit (if asked).

### Task 6: Real-SNID tests (marked, skipped in CI)

**Files:** Modify `pyproject.toml`, `tests/conftest.py`; Create `tests/test_snid_real.py`

- [ ] **Step 1:** Register the marker in `pyproject.toml` (`[tool.pytest.ini_options] markers = ["snid_real: runs real SNID-SAGE; needs templates; skipped in CI"]`). In `tests/conftest.py` add `pytest_collection_modifyitems` that skips `snid_real` tests when `os.environ.get("CI")` is set or `BOSSCAT_SNID_TEMPLATES` is unset, with a message naming the reason.
- [ ] **Step 2: Write the test** `test_real_snid_end_to_end` in `tests/test_snid_real.py`: pick the first file in `data/spectra/` (skip if none), read with `read_boss_spectrum`, run `make_snid(SnidConfig()).run(wave, flux, ivar, meta)`, assert `snid_success is True`, every key in the Task 3 Interfaces list is present, and `snid_consensus_redshift` is a float within `[-0.01, 0.7]`. A second test `test_snid_result_fields_exist` asserts every attribute read by `snid_result_to_fields` is in `dataclasses.fields(SNIDResult)`.
- [ ] **Step 3:** Run `uv run pytest -v` in CI-like mode (`CI=1 uv run pytest -v`): real tests show as skipped. Then locally with templates: `BOSSCAT_SNID_TEMPLATES=<path> uv run pytest -m snid_real -v`; expect PASS (~15 s).
- [ ] **Step 4:** Commit (if asked).

### Task 7: Evaluation against NGSF and TNS

**Files:** Create `bosscat/evaluate.py`; Test `tests/test_evaluate.py`

**Interfaces:**
- Produces:
  - `broad_type(label: str | None) -> str | None`: strip a leading `"SN "`; `Ia*` -> `Ia`, `Ib*` -> `Ib`, `Ic*` -> `Ic`, `II*` (including IIn, IIb, IIP) -> `II`, `SLSN*` -> `SLSN`; anything else or blank -> `None` (excluded from agreement).
  - `evaluate(results_csv, manifest_csv, high_z: float = 0.3) -> dict` joining on `results.filename == manifest.file_name + ".fits"`. Returns counts (`n_results`, `n_matched`, `n_snid_success`, `n_flagged` per flag), redshift residuals `snid_consensus_redshift - z_fit` (median and MAD) overall and for `z_fit > high_z`, and broad-type agreement fractions versus NGSF `sn_family` and versus `tns_type` (non-blank only), each with its denominator.
  - CLI `python -m bosscat.evaluate RESULTS.csv [--manifest data/manifest.csv]` printing the dict as readable text.

- [ ] **Step 1: Write failing tests** with tiny synthetic CSVs: `test_broad_type_mapping` (`"SN Ia-91T-like"` -> `Ia`, `"IIn"` -> `II`, `"SLSN-I"` -> `SLSN`, `""` -> `None`, `"Ia 99aa-like"` -> `Ia`); `test_join_on_filename_suffix`; `test_residual_stats_overall_and_high_z`; `test_agreement_denominators_exclude_blank_tns`; `test_failed_snid_rows_counted_not_scored`.
- [ ] **Step 2:** Run `uv run pytest tests/test_evaluate.py -v`; expect FAIL.
- [ ] **Step 3:** Implement with the standard library `csv` and numpy (no pandas).
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit (if asked).

### Task 8: Calibration experiments and defaults

**Files:** Modify `.gitignore` (`runs/`), `bosscat/config.py` (defaults), `CLAUDE.md`; Create `.claude/design/reports/2026-11-13-snid-calibration.md`

Runs use the user's machine and templates (about 11 s per spectrum, ~20 min per full run; do the runs in the background and record the time). Outputs go to `runs/` and are not committed.

- [ ] **Step 1:** Add `runs/` to `.gitignore`. Run the CLI on all 105 spectra for: `zmax=0.3, rebin 1`; `zmax=0.7` with `rebin` in `1, 2, 4, 8`. Run `python -m bosscat.evaluate` on each.
- [ ] **Step 2:** Decide, from the evaluations: (a) whether `zmax=0.7` introduces new false matches among the ~100 low-z objects (compare against `zmax=0.3`); (b) the rebin default. Rule: ship a factor above 1 only if agreement with NGSF/TNS is more stable than at factor 1 (state the numbers); otherwise keep `rebin_factor=1` and record that rebinning was dropped. Also compare SNID's log-grid spacing with BOSS's ~69 km/s pixels and note it.
- [ ] **Step 3:** Calibrate flags on the 5 objects with NGSF `z_fit > 0.3`: success means `snid_consensus_redshift` recovers `z_fit` (state the tolerance used and why) and the flags separate these from wrong or low-quality fits. Set `min_snr`, `zedge_margin`, `min_lap`, `min_hsigma_lap_ccc` defaults accordingly; state the uncertainty from n=5 and do not over-tune.
- [ ] **Step 4 (measure only, do not enable):** For air/vacuum and for dereddening with `boss_ebv`, report the size of the effect on one or two spectra. Enabling either needs user sign-off.
- [ ] **Step 5:** Write the report in `.claude/design/reports/` (results tables, decisions, caveats: NGSF is not truth, only ~19 TNS types). Update `CLAUDE.md` layout/rules to match the final code. Re-run `uv run pytest -v`; expect PASS.
- [ ] **Step 6:** Commit config, report, `.gitignore`, `CLAUDE.md` (if asked).

---

## Self-Review

- **Spec coverage:** blind classification (Tasks 3, 5 tests); `zmax` 0.7 and `show_plots=False` (Task 3); rebinning with drop option (Tasks 2, 8); namespaced output and sidecar (Tasks 3, 5); flags annotate-only (Task 4); real-SNID marked tests (Task 6); evaluation (Task 7); calibration, air/vacuum and dereddening measurement, defaults (Task 8); templates path out of code (Task 1).
- **Step scan:** thresholds `min_snr` and `zedge_margin` start as provisional defaults and are set in Task 8, not left as TBD.
- **Type consistency:** `Classifier.run` signature, `make_snid`, `snid_flags`, `rebin_loglam`, `base_row`, `write_provenance`, `run_pipeline`, and `broad_type` are named identically across tasks.
- **Review Focus:** items 1-5 map to tests in Tasks 3 and 4; item 5 maps to the Task 5 test.
- **Known judgement calls to verify while building:** `SNIDResult` has `redshift_error` but I have not checked its `lap`/`hsigma_lap_ccc` types or whether non-SN type names are exactly `AGN`/`Galaxy`/`Star`/`CV`; the real-SNID test in Task 6 and the check in Task 4 settle this.
