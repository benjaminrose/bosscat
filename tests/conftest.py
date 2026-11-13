"""Shared fixtures: build small synthetic BOSS-style coadd FITS files."""

import numpy as np
import pytest
from astropy.io import fits

# (SPALL column, FITS format, default value)
_SPALL_COLUMNS = [
    ("FIBER_RA", "D", 150.25),
    ("FIBER_DEC", "D", -2.5),
    ("Z", "E", 0.5),
    ("Z_ERR", "E", 0.001),
    ("CLASS", "6A", "GALAXY"),
    ("SUBCLASS", "21A", "STARFORMING"),
    ("ZWARNING", "K", 0),
    ("FIELD", "K", 100050),
    ("FIBERID_LIST", "10A", "379"),
    ("TOO_ID", "K", 61120000000435),
    ("PROGRAMNAME", "10A", "TOO_TNS"),
    ("EBV", "E", 0.0125),
    ("SN_MEDIAN_ALL", "E", 6.5),
]


def write_boss_fits(path, *, spall=None, drop=(), with_spall=True,
                    loglam=None, flux=None, ivar=None):
    """Write a minimal BOSS-like coadd FITS file.

    Parameters
    ----------
    path : pathlib.Path
        Output path.
    spall : dict, optional
        Overrides for SPALL column values, keyed by column name.
    drop : iterable of str, optional
        SPALL columns to leave out of the table entirely.
    with_spall : bool, optional
        If False, omit the SPALL extension.
    loglam, flux, ivar : array-like, optional
        COADD columns.  Defaults give 5 pixels with ``ivar[2] == 0``.

    Returns
    -------
    pathlib.Path
        *path*, for convenience.
    """
    loglam = np.array([3.55, 3.60, 3.65, 3.70, 3.75] if loglam is None else loglam)
    flux = np.array([1.0, 2.0, 3.0, 4.0, 5.0] if flux is None else flux)
    ivar = np.array([1.0, 2.0, 0.0, 4.0, 5.0] if ivar is None else ivar)

    primary = fits.PrimaryHDU()
    primary.header["MJD"] = 61120
    # Deliberately different from FIBER_RA/FIBER_DEC: the loader must ignore these.
    primary.header["RA"] = 1.0
    primary.header["DEC"] = 1.0

    coadd = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="FLUX", format="E", array=flux.astype("f4")),
            fits.Column(name="LOGLAM", format="E", array=loglam.astype("f4")),
            fits.Column(name="IVAR", format="E", array=ivar.astype("f4")),
        ],
        name="COADD",
    )
    hdus = [primary, coadd]

    if with_spall:
        overrides = spall or {}
        cols = []
        for name, fmt, default in _SPALL_COLUMNS:
            if name in drop:
                continue
            value = overrides.get(name, default)
            cols.append(fits.Column(name=name, format=fmt, array=np.array([value])))
        hdus.append(fits.BinTableHDU.from_columns(cols, name="SPALL"))

    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


@pytest.fixture
def make_fits(tmp_path):
    """Return a factory ``make_fits(name="a.fits", **kwargs) -> Path``."""

    def _make(name="spec.fits", **kwargs):
        return write_boss_fits(tmp_path / name, **kwargs)

    return _make
