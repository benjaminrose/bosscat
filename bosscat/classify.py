"""SNID-SAGE classifier wrapper."""

# TODO: replace with env var or config file before production
# These can be downloaded with the snid-sage provided `snid-sage-download-templates`
TEMPLATES_DIR = "/Users/ben_rose/Documents/1 Projects/as5-classifier/SNID-SAGE/templates"

# NOTE: verify these import paths after `pip install snid-sage`.
# The API docs example uses `from snid_sage.snid.snid import run_snid`;
# preprocess_spectrum and run_snid_analysis are expected in the same module.
from snid_sage.snid.snid import preprocess_spectrum, run_snid_analysis


def classify_spectrum(wave, flux):
    """Run SNID-SAGE on a wavelength/flux pair.

    Calls :func:`snid_sage.snid.snid.preprocess_spectrum` followed by
    :func:`snid_sage.snid.snid.run_snid_analysis`.  The template library
    at :data:`TEMPLATES_DIR` is used for matching.

    Parameters
    ----------
    wave : array-like
        Wavelength array in Angstroms.
    flux : array-like
        Flux array.  Units are not critical; SNID-SAGE normalises
        internally.

    Returns
    -------
    result : SNIDResult
        Classification result containing consensus type, redshift, age,
        type fractions, and match quality metrics.
    trace : dict
        Internal SNID-SAGE diagnostic dictionary, useful for debugging.
    """
    processed, _pre_trace = preprocess_spectrum(input_spectrum=(wave, flux))
    result, trace = run_snid_analysis(processed, TEMPLATES_DIR)
    return result, trace
