"""Read BOSS FITS spectra (coadded spectrum only)."""

import math
from pathlib import Path

import numpy as np
from astropy.io import fits

# SPALL column -> meta key.  A missing column or a blank/NaN value gives None.
_SPALL_OPTIONAL = {
    "Z": "z",
    "Z_ERR": "z_err",
    "CLASS": "class",
    "SUBCLASS": "subclass",
    "ZWARNING": "zwarning",
    "FIELD": "field",
    "FIBERID_LIST": "fiberid",
    "TOO_ID": "too_id",
    "PROGRAMNAME": "programname",
    "EBV": "ebv",
    "SN_MEDIAN_ALL": "sn_median",
}


def _clean(value):
    """Convert a FITS table value to a plain Python value.

    Parameters
    ----------
    value : object
        A scalar read from a FITS table row (NumPy scalar, bytes, str, ...).

    Returns
    -------
    object or None
        The value as a built-in Python type, or ``None`` if it is a blank
        string or a non-finite float.
    """
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode()
    if isinstance(value, str):
        return value.strip() or None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def read_boss_spectrum(filepath):
    """Read the coadded spectrum from a BOSS coadd FITS file.

    Reads the ``COADD`` BinTable extension for the spectrum and the
    single-row ``SPALL`` extension for object-level metadata (e.g.
    ``ToO-specLite-FIELD-MJD-CATALOGID.fits``).

    Parameters
    ----------
    filepath : str or Path
        Path to a BOSS FITS file.

    Returns
    -------
    wave : numpy.ndarray
        Observed-frame wavelength in Angstroms, shape ``(N,)``.
    flux : numpy.ndarray
        Flux in units of 1e-17 erg/s/cm²/Å, shape ``(N,)``.
    ivar : numpy.ndarray
        Inverse variance of ``flux``, shape ``(N,)``.  Pixels with
        ``ivar == 0`` are masked by the BOSS pipeline.
    meta : dict
        Object-level metadata.  Always present: ``filename``, ``filepath``,
        ``mjd``, ``ra``, ``dec`` (from ``FIBER_RA``/``FIBER_DEC``).
        Optional (``None`` if absent or blank in ``SPALL``): ``z``,
        ``z_err`` (also ``None`` if negative, a BOSS "no error" sentinel),
        ``class``, ``subclass``, ``zwarning``, ``field``,
        ``fiberid``, ``too_id``, ``programname``, ``ebv``, ``sn_median``.

    Raises
    ------
    ValueError
        If the file has no ``SPALL`` extension, so no coordinates.
    KeyError
        If ``SPALL`` lacks ``FIBER_RA`` or ``FIBER_DEC``.
    """
    filepath = Path(filepath)
    with fits.open(filepath) as hdul:
        coadd = hdul["COADD"].data
        wave = 10.0 ** coadd["loglam"].astype(float)
        flux = coadd["flux"].astype(float)
        ivar = coadd["ivar"].astype(float)

        if "SPALL" not in hdul:
            raise ValueError(
                f"{filepath.name}: no SPALL extension "
                "(needed for FIBER_RA/FIBER_DEC)"
            )
        spall = hdul["SPALL"].data[0]
        cols = hdul["SPALL"].columns.names

        meta = {
            "filename": filepath.name,
            "filepath": str(filepath),
            "mjd": hdul[0].header.get("MJD", ""),
            "ra": float(spall["FIBER_RA"]),
            "dec": float(spall["FIBER_DEC"]),
        }
        for col, key in _SPALL_OPTIONAL.items():
            meta[key] = _clean(spall[col]) if col in cols else None

        # BOSS uses negative z_err (e.g. -1) as a "no error" sentinel
        if meta["z_err"] is not None and meta["z_err"] < 0:
            meta["z_err"] = None

    return wave, flux, ivar, meta
