"""Write classification results to CSV."""

import csv
import datetime
import importlib.metadata
import json
import platform
from pathlib import Path


PROVENANCE_SUFFIX = ".provenance.json"


def base_row(meta):
    """Build the shared and ``boss_*`` columns of a result row.

    Parameters
    ----------
    meta : dict
        Metadata dict returned by
        :func:`~bosscat.io.boss.read_boss_spectrum`.

    Returns
    -------
    row : dict
        Flat dict of metadata columns; classifier columns are added by the
        caller.
    """
    return {
        "filename": meta.get("filename", ""),
        "filepath": meta.get("filepath", ""),
        "field": meta.get("field", ""),
        "mjd": meta.get("mjd", ""),
        "fiberid": meta.get("fiberid", ""),
        "ra": meta.get("ra", ""),
        "dec": meta.get("dec", ""),
        "boss_z": meta.get("z", ""),
        "boss_z_err": meta.get("z_err", ""),
        "boss_zwarning": meta.get("zwarning", ""),
        "boss_class": meta.get("class", ""),
        "boss_subclass": meta.get("subclass", ""),
        "too_id": meta.get("too_id", ""),
        "programname": meta.get("programname", ""),
        "ebv": meta.get("ebv", ""),
        "sn_median": meta.get("sn_median", ""),
    }


def sidecar_path(output_path):
    """Return the provenance sidecar path for a results file.

    Parameters
    ----------
    output_path : str or Path
        Results CSV path.

    Returns
    -------
    Path
        ``<output file name>.provenance.json`` in the same directory, so it
        can never coincide with the results file (``results.csv`` gives
        ``results.csv.provenance.json``).
    """
    output_path = Path(output_path)
    return output_path.with_name(output_path.name + PROVENANCE_SUFFIX)


def write_provenance(output_path, classifiers, input_files, complete=True,
                     n_rows=0):
    """Write a JSON sidecar describing a run, next to the results CSV.

    Parameters
    ----------
    output_path : str or Path
        Results CSV path; see :func:`sidecar_path` for the sidecar name.
    classifiers : sequence of Classifier
        Classifiers used; each ``provenance()`` is recorded under its name.
    input_files : sequence of str
        Input file paths.
    complete : bool, optional
        False when the run stopped early (interrupt or crash).
    n_rows : int, optional
        Number of result rows written.

    Returns
    -------
    Path
        Path of the JSON file written.
    """
    path = sidecar_path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        version = importlib.metadata.version("bosscat")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    record = {
        "bosscat_version": version,
        "python_version": platform.python_version(),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "input_files": [str(f) for f in input_files],
        "classifiers": {c.name: c.provenance() for c in classifiers},
        "complete": bool(complete),
        "n_rows": int(n_rows),
    }
    path.write_text(json.dumps(record, indent=2) + "\n")
    return path


def write_results(rows, output_path):
    """Write a list of row dicts to a CSV file.

    The CSV column set is the union of all keys across every row so that
    the schema is self-consistent within a batch run even when some spectra
    fail and produce only a subset of columns.  Missing values are written
    as empty strings.

    Parameters
    ----------
    rows : list of dict
        Row dicts, typically produced by :func:`base_row` plus classifier fields.
    output_path : str or Path
        Destination CSV path.  Parent directories are created if they do
        not exist.
    """
    if not rows:
        print("No results to write.")
        return

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect all fieldnames in insertion order
    seen = {}
    for row in rows:
        for k in row:
            seen[k] = None
    fieldnames = list(seen.keys())

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,
                                extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} row(s) to {output_path}")
