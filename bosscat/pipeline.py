"""End-to-end pipeline: BOSS FITS → classifier(s) → CSV."""

from .classify import classify_spectrum
from .io.boss import read_boss_spectrum
from .io.results import snid_result_to_row, write_results


def run_pipeline(input_files, output_path, run_snid=True):
    """Classify a list of BOSS FITS files and write results to a CSV.

    Iterates over *input_files*, reading each coadded spectrum, running the
    requested classifier(s), and collecting one row per file.  Errors at the
    read or classify stage are caught, printed, and recorded in the output
    rather than aborting the batch.

    Parameters
    ----------
    input_files : list of str
        Paths to BOSS coadd FITS files.
    output_path : str
        Destination CSV file path.
    run_snid : bool, optional
        Whether to run the SNID-SAGE classifier.  Default is ``True``.
    """
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

        snid_result = None
        if run_snid:
            try:
                good = ivar > 0  # drop pixels masked by the BOSS pipeline
                snid_result, _ = classify_spectrum(wave[good], flux[good])
            except Exception as exc:
                print(f"  [ERROR] SNID-SAGE failed on {filepath}: {exc}")
                meta["snid_error"] = str(exc)

        row = snid_result_to_row(meta, snid_result)
        if "snid_error" in meta:
            row["snid_error"] = meta["snid_error"]

        rows.append(row)

    write_results(rows, output_path)
