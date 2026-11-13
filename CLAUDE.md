# bosscat: project notes for Claude

BOSS Classification Analysis Toolkit. Python package `bosscat` (src layout not used; package lives at `bosscat/`). Research code for supernova/spectra classification.

## Goals and context

- Classify BOSS ToO (transient) coadd spectra (`ToO-specLite-FIELD-MJD-CATALOGID.fits`) with multiple tools. SNID-SAGE is first; NGSF and others are expected later.
- `data/manifest.csv` (105 rows) is the existing **NGSF** output for the 105 spectra in `data/spectra/` (columns `name, file_name, snr_median, tns_type, sn_template, sn_family, z_fit, ...`). Use it as a comparison, not ground truth. Only ~19 rows have a TNS type. `file_name` omits `.fits`.
- BOSS `Z`/`CLASS` in SPALL are the BOSS pipeline's own values (QSO/GALAXY/STAR), NOT NGSF's: only 45 of 105 agree with the manifest `z_fit` within 0.01. The NGSF reference is the manifest `z_fit`/`sn_family`; `boss_*` columns are context only. Never pass any redshift into a classifier (`forced_redshift`); production inputs will not have one.
- Typical data are low S/N (manifest `snr_median` median ~3.6, range 1-20). Redshifts: 5 of 105 have NGSF z > 0.3 (max 0.67).
- Design doc: `docs/superpowers/specs/2026-11-13-snid-classification-design.md`.

## Design rules (decided with the user)

- Output columns are namespaced by their source: `boss_*` for BOSS pipeline values, `snid_*` for SNID-SAGE, `ngsf_*` for any future NGSF. Never put a SNID value in a column named for another tool, and never leave placeholder columns for fields a tool cannot produce (e.g. SNID has no `A_v`, `host_template`, `frac_sn`).
- Use each tool's native field names and report both of its values when it has two (e.g. SNID `redshift` and `consensus_redshift`), each with its own error.
- Each classifier takes `(wave, flux, ivar, meta)` and returns a dict of already-prefixed fields; the pipeline merges dicts. Record each tool's configuration (version, templates, z range, rebin factor) in a sidecar JSON.
- SNID redshift search `zmax` is ~0.7 (covers the manifest max, 0.667); the 5 manifest objects with z > 0.3 are the calibration set for what counts as a good fit. Spectra that may lie beyond `zmax` or fit poorly get `snid_flag_*` annotations. Flags only annotate: never suppress or alter `snid_consensus_type`.
- Rebinning ships only if it improves stability against NGSF/TNS; dropping it is acceptable.
- Using the NGSF/BOSS redshift for calibration and evaluation is approved; using it as a classifier input is not.

## SNID-SAGE facts (verified by running v1.2.5, 2026-11-13)

- `run_snid_analysis` defaults to `show_plots=True`, and its plotting step raises `ImportError` in this install. Always pass `show_plots=False`. A run takes ~11 s per spectrum.
- `SNIDResult` has no `age`, `age_error`; use `consensus_age`, `consensus_age_error`. Also `redshift` (best template) differs from `consensus_redshift`.
- `TEMPLATES_DIR` in `bosscat/classify.py` is hardcoded to a path outside this repo; replace with config.
- Mocked tests hid these API mismatches. Real-SNID tests get a pytest marker, are skipped in CI (templates unavailable there), and run locally.

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
- `bosscat/io/results.py`: `snid_result_to_row`, `write_results` (CSV).
- `bosscat/classify.py`: `classify_spectrum` wraps SNID-SAGE (`preprocess_spectrum` then `run_snid_analysis`).
- `bosscat/pipeline.py`: `run_pipeline`; `bosscat/__main__.py`: CLI.

## Testing

- pytest (dev dependency group); run with `uv run pytest`. CI: `.github/workflows/tests.yml`.
- `tests/conftest.py` builds synthetic BOSS-style FITS files. SNID-SAGE is mocked in pipeline tests; there is no test that runs real SNID-SAGE yet.
- Ask before adding further test dependencies.

## Bibliography

- Never hand-write BibTeX. Every `.bib` entry comes verbatim from NASA ADS via `/ads-cite` or the `ads-cite` CLI, using `--rekey` (citekeys like `Surname_Subject_Year`).
- Use `natbib` or `biblatex` with `.bib` files; do not hardcode citations in `.tex`.

## Git

- Package code, tests, and CI are committed on `main`. Untracked: `.claude/`, `tns_2025ahsk_output_example.ascii`.
- Do not commit or push without an explicit request.
