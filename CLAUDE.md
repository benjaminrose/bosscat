# bosscat: project notes for Claude

BOSS Classification Analysis Toolkit. Python package `bosscat` (src layout not used; package lives at `bosscat/`). Research code for supernova/spectra classification.

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

## Testing

- `tests/` is currently empty. No test runner is configured yet. Ask before adding a test dependency.

## Bibliography

- Never hand-write BibTeX. Every `.bib` entry comes verbatim from NASA ADS via `/ads-cite` or the `ads-cite` CLI, using `--rekey` (citekeys like `Surname_Subject_Year`).
- Use `natbib` or `biblatex` with `.bib` files; do not hardcode citations in `.tex`.

## Git

- Current state: the package code (`bosscat/`, `pyproject.toml`, `uv.lock`) is not yet committed, and `src/` is deleted in the working tree. Confirm with the user before committing or changing history.
- Do not commit or push without an explicit request.
