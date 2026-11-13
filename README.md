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
any spectrum is processed.

Classification is blind: the BOSS redshift is recorded in `boss_z` but never
given to a classifier.

### Output

One CSV row per input file. Shared columns (`filename`, `filepath`, `field`,
`mjd`, `fiberid`, `ra`, `dec`, `too_id`, `programname`, `ebv`, `sn_median`)
come first, then `boss_*` (BOSS pipeline values: `boss_z`, `boss_z_err`,
`boss_zwarning`, `boss_class`, `boss_subclass`), then one namespace per
classifier (`snid_*`, including `snid_frac_<type>`, `snid_flag_*`). A
classifier that raises on a spectrum yields `<name>_success=False` and
`<name>_error`. Files that cannot be read yield a row with `filepath` and
`error` only. Columns are the union across rows; absent values are empty.

A sidecar `<output stem>.json` records `bosscat_version`, `python_version`,
UTC `created`, `input_files`, and each classifier's provenance (for SNID: the
`snid-sage` version, config, and templates directory).
