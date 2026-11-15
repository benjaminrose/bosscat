# bosscat, BOSS Classification Analysis Toolkit

## Usage

```
uv run bosscat FILE.fits [FILE.fits ...] [-o results.csv]
```

| Flag | Meaning |
| --- | --- |
| `--classifier NAME` | Classifier to run; repeatable. Available: `snid` (default). |
| `--no-classifier` | Run no classifier; write BOSS metadata only. |
| `--zmax Z` | SNID maximum redshift (default from `SnidConfig`). |
| `--rebin-factor N` | SNID rebinning factor (default from `SnidConfig`). |
| `--templates-dir DIR` | SNID templates directory. |
| `-o, --output PATH` | Output CSV (default `results.csv`). |

SNID templates: pass `--templates-dir` or set `BOSSCAT_SNID_TEMPLATES`. If
neither points to an existing directory, the run exits with status 2 before
any spectrum is processed. SNID takes about 11 s per spectrum.

Classification is blind: the BOSS redshift is recorded in `boss_z` but never
given to a classifier.

### Output

One CSV row per input file. Columns are the union across rows, in first-seen
order, so the header depends on the first row: shared columns (`filename`,
`filepath`, `field`, `mjd`, `fiberid`, `ra`, `dec`, `too_id`, `programname`,
`ebv`, `sn_median`) do NOT always come first (an unreadable first input, or an
error row, changes the header). A normal row has the shared columns, then
`boss_*` (BOSS pipeline values: `boss_z`, `boss_z_err`, `boss_zwarning`,
`boss_class`, `boss_subclass`), then one namespace per classifier (`snid_*`:
result fields, `snid_snr_median`, `snid_flag_*`, and `snid_frac_<type>`, which
can appear after the flag columns). A classifier that raises on a spectrum
yields `<name>_success=False` and `<name>_error`. Files that cannot be read
yield a row with `filepath` and `error` only. Absent values are empty.

`snid_consensus_age`, `snid_consensus_age_err` and `snid_r` are intentionally
not emitted: SNID-SAGE 1.2.5 returns constant 0.0 for them.

A sidecar `<output file name>.provenance.json` (`results.csv` gives
`results.csv.provenance.json`) records `bosscat_version`, `python_version`,
UTC `created`, `input_files`, each classifier's provenance (for SNID: the
`snid-sage` version, config, and templates directory), `complete` and
`n_rows`. If a run is interrupted (Ctrl-C) or crashes, the rows collected so
far are still written and `complete` is false. Output names ending in
`.provenance.json` are rejected.

### Evaluation

```
uv run python -m bosscat.evaluate results.csv [--manifest data/manifest.csv] [--json]
```

Compares SNID results with the NGSF manifest (`z_fit`, `sn_family`) and TNS
types: redshift residuals, broad-type agreement, failure and flag counts, and
`unmapped_reference_labels` (reference labels that map to no broad type and are
therefore excluded from agreement).

## Testing and coverage

```
uv run pytest                 # tests only
uv run pytest --cov           # tests plus per-module coverage (branch coverage, missing lines)
```

Real SNID-SAGE tests (marker `snid_real`) run only when `BOSSCAT_SNID_TEMPLATES`
points at the templates and `CI` is unset; to include them locally:

```
MPLBACKEND=Agg BOSSCAT_SNID_TEMPLATES=/path/to/templates uv run pytest --cov
```

CI runs `uv run pytest --cov --cov-fail-under=95`. CI has no templates, so the
real-SNID tests are skipped there and the CI number is the one that gates (96.4%
when the threshold was set; the local run with real SNID gave the same figure).
Coverage settings live in `pyproject.toml` under `[tool.coverage.*]`. The
`.coverage` data file is gitignored. The threshold is deliberately a CI flag, not
a `pyproject.toml` setting, so running `pytest --cov` on a single test file does
not fail locally.
