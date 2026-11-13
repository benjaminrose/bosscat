"""End-to-end pipeline: BOSS FITS → classifier(s) → CSV."""

from dataclasses import replace
from functools import partial
from pathlib import Path

from .io.boss import read_boss_spectrum
from .io.results import (PROVENANCE_SUFFIX, base_row, write_provenance,
                         write_results)


# BOSS pipeline values that must never reach a classifier (blind rule).
BLIND_META_KEYS = ("z", "z_err", "class", "subclass", "zwarning")


def run_pipeline(input_files, output_path, classifiers=()):
    """Classify a list of BOSS FITS files and write results to a CSV.

    Reads each coadded spectrum, runs every classifier, and collects one row
    per file. A JSON provenance sidecar is written next to the CSV. Errors
    at the read stage or inside a classifier are recorded in the row rather
    than aborting the batch. Classifiers receive ``meta`` without the BOSS
    ``z``, ``z_err``, ``class``, ``subclass`` and ``zwarning`` keys; the
    ``boss_*`` columns are still written. If the run is interrupted
    (``KeyboardInterrupt``) or an unexpected exception escapes the loop, the
    rows collected so far are written, the sidecar records
    ``"complete": false``, and the exception is re-raised.

    Parameters
    ----------
    input_files : sequence of str
        Paths to BOSS coadd FITS files.
    output_path : str or Path
        Destination CSV file path.
    classifiers : sequence of Classifier, optional
        Classifiers to run on each spectrum. Default is none.

    Raises
    ------
    RuntimeError
        Propagated from a classifier's ``provenance()`` (for example a
        missing templates directory). It is evaluated once before any file
        is processed, so setup problems fail fast.
    ValueError
        If ``output_path`` ends with ``.provenance.json`` (reserved for
        sidecars).
    """
    # sidecar_path() appends a suffix, so an output can never equal its own
    # sidecar; reject names that look like a sidecar instead.
    if Path(output_path).name.endswith(PROVENANCE_SUFFIX):
        raise ValueError(
            f"output path {output_path} ends with {PROVENANCE_SUFFIX!r}, "
            f"which is reserved for provenance sidecars; choose another name"
        )
    classifiers = list(classifiers)
    # Fail fast, and call each provenance() exactly once: the sidecar reuses it.
    frozen = [replace(c, provenance=partial(dict, c.provenance()))
              for c in classifiers]
    rows = []

    def finish(complete):
        write_results(rows, output_path)
        write_provenance(output_path, frozen, input_files,
                         complete=complete, n_rows=len(rows))

    try:
        for filepath in input_files:
            print(f"Processing {filepath} ...")
            meta = {"filepath": str(filepath)}

            try:
                wave, flux, ivar, meta = read_boss_spectrum(filepath)
            except Exception as exc:
                print(f"  [ERROR] Could not read {filepath}: {exc}")
                rows.append({"filepath": str(filepath), "error": str(exc)})
                continue

            row = base_row(meta)
            blind = {k: v for k, v in meta.items() if k not in BLIND_META_KEYS}
            for clf in classifiers:
                try:
                    row.update(clf.run(wave, flux, ivar, dict(blind)))
                except Exception as exc:
                    print(f"  [ERROR] {clf.name} failed on {filepath}: {exc}")
                    row[f"{clf.name}_error"] = str(exc)
                    row[f"{clf.name}_success"] = False
            rows.append(row)
    except BaseException:
        finish(complete=False)
        raise
    finish(complete=True)
