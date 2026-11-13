"""Evaluate SNID results against the NGSF manifest and TNS types."""

import argparse
import csv
import json

import numpy as np

REDSHIFT_COLUMNS = ("snid_consensus_redshift", "snid_redshift")


def broad_type(label):
    """Map a detailed supernova label to a broad class.

    Parameters
    ----------
    label : str or None
        Label such as ``"SN Ia-91T-like"``, ``"IIn"`` or ``"SLSN-I"``.

    Returns
    -------
    str or None
        One of ``"Ia"``, ``"Ib"``, ``"Ic"``, ``"II"``, ``"SLSN"``, or None
        when the label is blank or not recognised.
    """
    if label is None:
        return None
    text = str(label).strip()
    if text.startswith("SN "):
        text = text[3:].strip()
    # Prefixes are distinct ("Ia", "Ib", "Ic", "II"), so order is safe
    # except SLSN, which does not start with any of them.
    if text.startswith("SLSN"):
        return "SLSN"
    for prefix, broad in (("Ia", "Ia"), ("Ib", "Ib"), ("Ic", "Ic"),
                          ("II", "II")):
        if text.startswith(prefix):
            return broad
    return None


def _float(value):
    """Parse a CSV string to a finite float, or None if blank/invalid."""
    if value is None:
        return None
    try:
        x = float(str(value).strip())
    except ValueError:
        return None
    return x if np.isfinite(x) else None


def _true(value):
    """Return True when a CSV string is ``"True"`` (case-insensitive)."""
    return str(value).strip().lower() == "true"


def _stats(values):
    """Median and MAD (about the median) of residuals, None when empty."""
    if not values:
        return {"n": 0, "median": None, "mad": None}
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    return {"n": int(arr.size), "median": med,
            "mad": float(np.median(np.abs(arr - med)))}


def _agreement(pairs):
    """Agreement count and fraction from (snid, reference) pairs."""
    agree = sum(1 for s, r in pairs if s == r)
    n = len(pairs)
    return {"n": n, "agree": agree, "fraction": agree / n if n else None}


def _read(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader), list(reader.fieldnames or [])


def evaluate(results_csv, manifest_csv, high_z=0.3):
    """Compare SNID results with NGSF (manifest) and TNS references.

    Results are joined to the manifest on ``filename == file_name + ".fits"``.
    The NGSF references are ``z_fit`` and ``sn_family``; the TNS reference is
    ``tns_type``. Rows whose ``snid_success`` is not true are counted but not
    scored.

    Parameters
    ----------
    results_csv : str or Path
        Results CSV written by :func:`bosscat.io.results.write_results`.
    manifest_csv : str or Path
        Manifest CSV with ``file_name``, ``z_fit``, ``sn_family`` and
        ``tns_type`` columns.
    high_z : float, optional
        Rows with manifest ``z_fit > high_z`` form the high-redshift subset.

    Returns
    -------
    dict
        JSON-serializable summary: ``n_results``, ``n_matched``,
        ``n_unmatched_results``, ``n_unmatched_manifest``, ``n_snid_success``,
        ``n_flagged`` (per ``snid_flag_*`` column), ``high_z``, ``redshift``
        (per SNID redshift column, with ``all`` and ``high_z`` entries of
        ``n``/``median``/``mad`` for ``snid_* - z_fit``) and
        ``type_agreement`` (``vs_ngsf`` and ``vs_tns`` with ``n``, ``agree``,
        ``fraction``).
    """
    results, res_cols = _read(results_csv)
    manifest, _ = _read(manifest_csv)
    by_file = {row["file_name"] + ".fits": row for row in manifest}

    flag_cols = [c for c in res_cols if c.startswith("snid_flag_")]
    n_flagged = {c: 0 for c in flag_cols}
    resid = {c: {"all": [], "high_z": []} for c in REDSHIFT_COLUMNS}
    pairs_ngsf, pairs_tns = [], []
    matched_files = set()
    n_matched = 0
    n_success = 0

    for row in results:
        ref = by_file.get(row.get("filename", ""))
        if ref is not None:
            matched_files.add(row["filename"])
            n_matched += 1
        if not _true(row.get("snid_success")):
            continue
        n_success += 1
        for c in flag_cols:
            if _true(row.get(c)):
                n_flagged[c] += 1
        if ref is None:
            continue
        z_ref = _float(ref.get("z_fit"))
        if z_ref is not None:
            for c in REDSHIFT_COLUMNS:
                z = _float(row.get(c))
                if z is None:
                    continue
                resid[c]["all"].append(z - z_ref)
                if z_ref > high_z:
                    resid[c]["high_z"].append(z - z_ref)
        snid_type = (row.get("snid_consensus_type") or "").strip()
        if snid_type:
            for pairs, label in ((pairs_ngsf, ref.get("sn_family")),
                                 (pairs_tns, ref.get("tns_type"))):
                mapped = broad_type(label)
                if mapped is not None:
                    pairs.append((snid_type, mapped))

    return {
        "n_results": len(results),
        "n_matched": n_matched,
        "n_unmatched_results": len(results) - n_matched,
        "n_unmatched_manifest": len(set(by_file) - matched_files),
        "n_snid_success": n_success,
        "n_flagged": n_flagged,
        "high_z": high_z,
        "redshift": {c: {k: _stats(v) for k, v in resid[c].items()}
                     for c in REDSHIFT_COLUMNS},
        "type_agreement": {"vs_ngsf": _agreement(pairs_ngsf),
                           "vs_tns": _agreement(pairs_tns)},
    }


def _fmt(x, spec=".4f"):
    return "n/a" if x is None else format(x, spec)


def format_report(summary):
    """Render an :func:`evaluate` summary as aligned text.

    Parameters
    ----------
    summary : dict
        Output of :func:`evaluate`.

    Returns
    -------
    str
        Multi-line report.
    """
    lines = []
    for key in ("n_results", "n_matched", "n_unmatched_results",
                "n_unmatched_manifest", "n_snid_success"):
        lines.append(f"{key:<24}{summary[key]:>8}")
    for col, n in summary["n_flagged"].items():
        lines.append(f"{col:<24}{n:>8}")
    lines.append("")
    lines.append(f"Redshift residual (SNID - NGSF z_fit), high_z = "
                 f"z_fit > {summary['high_z']}")
    lines.append(f"{'':<26}{'subset':<8}{'n':>6}{'median':>10}{'MAD':>10}")
    for col, subsets in summary["redshift"].items():
        for name, s in subsets.items():
            lines.append(f"{col:<26}{name:<8}{s['n']:>6}"
                         f"{_fmt(s['median']):>10}{_fmt(s['mad']):>10}")
    lines.append("")
    lines.append("Broad-type agreement (SNID consensus type)")
    lines.append(f"{'':<10}{'n':>6}{'agree':>7}{'fraction':>10}")
    for name, s in summary["type_agreement"].items():
        lines.append(f"{name:<10}{s['n']:>6}{s['agree']:>7}"
                     f"{_fmt(s['fraction'], '.3f'):>10}")
    return "\n".join(lines)


def parse_args(argv=None):
    """Parse command-line arguments.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; defaults to ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Attributes ``results``, ``manifest``, ``high_z`` and ``json``.
    """
    parser = argparse.ArgumentParser(
        description="Compare SNID results with the NGSF manifest and TNS types",
    )
    parser.add_argument("results", help="Results CSV from bosscat")
    parser.add_argument("--manifest", default="data/manifest.csv",
                        help="Manifest CSV (default: %(default)s)")
    parser.add_argument("--high-z", type=float, default=0.3,
                        help="High-redshift cut on z_fit (default: %(default)s)")
    parser.add_argument("--json", action="store_true",
                        help="Print the summary as JSON")
    return parser.parse_args(argv)


def main():
    """Entry point for ``python -m bosscat.evaluate``."""
    args = parse_args()
    summary = evaluate(args.results, args.manifest, high_z=args.high_z)
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(format_report(summary))


if __name__ == "__main__":
    main()
