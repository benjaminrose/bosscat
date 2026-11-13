"""End-to-end pipeline: BOSS FITS → classifier(s) → CSV."""

from dataclasses import replace
from functools import partial

from .io.boss import read_boss_spectrum
from .io.results import base_row, write_provenance, write_results


def run_pipeline(input_files, output_path, classifiers=()):
    """Classify a list of BOSS FITS files and write results to a CSV.

    Reads each coadded spectrum, runs every classifier, and collects one row
    per file. A JSON provenance sidecar is written next to the CSV. Errors
    at the read stage or inside a classifier are recorded in the row rather
    than aborting the batch.

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
    """
    classifiers = list(classifiers)
    # Fail fast, and call each provenance() exactly once: the sidecar reuses it.
    frozen = [replace(c, provenance=partial(dict, c.provenance()))
              for c in classifiers]
    rows = []

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
        for clf in classifiers:
            try:
                row.update(clf.run(wave, flux, ivar, meta))
            except Exception as exc:
                print(f"  [ERROR] {clf.name} failed on {filepath}: {exc}")
                row[f"{clf.name}_error"] = str(exc)
                row[f"{clf.name}_success"] = False
        rows.append(row)

    write_results(rows, output_path)
    write_provenance(output_path, frozen, input_files)
