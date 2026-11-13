"""Write classification results to CSV."""

import csv
from pathlib import Path


def snid_result_to_row(meta, result):
    """Flatten a ``SNIDResult`` and metadata dict into a single CSV row.

    All SNID-SAGE columns are prefixed with ``snid_`` so that future
    classifiers can add their own columns without name collisions.  Type
    fraction columns are named ``snid_frac_<type>`` where ``<type>`` is
    taken directly from ``SNIDResult.type_fractions`` (forward-slashes and
    spaces replaced with underscores).

    Parameters
    ----------
    meta : dict
        Metadata dict returned by
        :func:`~bosscat.io.boss.read_boss_spectrum`.
    result : SNIDResult or None
        Classification result from SNID-SAGE.  Pass ``None`` when the
        classifier did not run or raised an exception; all SNID columns
        will be omitted from the row.

    Returns
    -------
    row : dict
        Flat dict suitable for writing with :class:`csv.DictWriter`.
    """
    row = {
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

    if result is None:
        row["snid_success"] = False
        return row

    row["snid_success"] = result.success
    row["snid_consensus_type"] = result.consensus_type
    row["snid_redshift"] = result.redshift
    row["snid_redshift_err"] = result.redshift_error
    row["snid_age"] = result.age
    row["snid_age_err"] = result.age_error
    row["snid_lap"] = result.lap

    # Type fractions — keyed by whatever types SNID found in this run
    fractions = result.type_fractions or {}
    for t, frac in fractions.items():
        safe_key = t.replace("/", "_").replace(" ", "_")
        row[f"snid_frac_{safe_key}"] = frac

    return row


def write_results(rows, output_path):
    """Write a list of row dicts to a CSV file.

    The CSV column set is the union of all keys across every row so that
    the schema is self-consistent within a batch run even when some spectra
    fail and produce only a subset of columns.  Missing values are written
    as empty strings.

    Parameters
    ----------
    rows : list of dict
        Row dicts, typically produced by :func:`snid_result_to_row`.
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
