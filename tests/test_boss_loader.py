"""Tests for bosscat.io.boss (FITS -> in-memory loader)."""

import glob
from pathlib import Path

import numpy as np
import pytest

from bosscat.io.boss import _clean, read_boss_spectrum

OPTIONAL_KEYS = [
    "z", "z_err", "class", "subclass", "zwarning", "field", "fiberid",
    "too_id", "programname", "ebv", "sn_median",
]
REAL_FILES = sorted(glob.glob(str(Path(__file__).parent.parent / "data/spectra/*.fits")))


# --- spectrum arrays -------------------------------------------------------

def test_arrays_are_converted(make_fits):
    wave, flux, ivar, _ = read_boss_spectrum(make_fits())
    np.testing.assert_allclose(wave, 10.0 ** np.array([3.55, 3.60, 3.65, 3.70, 3.75]),
                               rtol=1e-5)
    np.testing.assert_allclose(flux, [1, 2, 3, 4, 5])
    np.testing.assert_allclose(ivar, [1, 2, 0, 4, 5])
    for arr in (wave, flux, ivar):
        assert arr.dtype == np.float64
        assert arr.shape == (5,)


def test_zero_ivar_pixels_are_kept_by_loader(make_fits):
    """Masking is the pipeline's job; the loader must return every pixel."""
    wave, flux, ivar, _ = read_boss_spectrum(make_fits())
    assert len(wave) == len(flux) == len(ivar) == 5
    assert (ivar == 0).sum() == 1


# --- required metadata -----------------------------------------------------

def test_required_metadata(make_fits):
    path = make_fits("ToO-specLite-1-2-3.fits")
    _, _, _, meta = read_boss_spectrum(path)
    assert meta["filename"] == "ToO-specLite-1-2-3.fits"
    assert meta["filepath"] == str(path)
    assert meta["mjd"] == 61120


def test_ra_dec_come_from_fiber_columns_not_primary_header(make_fits):
    _, _, _, meta = read_boss_spectrum(make_fits(spall={"FIBER_RA": 274.25,
                                                        "FIBER_DEC": 26.5}))
    assert meta["ra"] == pytest.approx(274.25)
    assert meta["dec"] == pytest.approx(26.5)  # primary header says 1.0


def test_accepts_str_and_path(make_fits):
    path = make_fits()
    a = read_boss_spectrum(path)[3]
    b = read_boss_spectrum(str(path))[3]
    assert a == b


# --- optional metadata -----------------------------------------------------

def test_optional_values_populated_with_plain_python_types(make_fits):
    meta = read_boss_spectrum(make_fits())[3]
    assert meta["z"] == pytest.approx(0.5)
    assert meta["z_err"] == pytest.approx(0.001)
    assert meta["class"] == "GALAXY"
    assert meta["subclass"] == "STARFORMING"
    assert meta["zwarning"] == 0
    assert meta["field"] == 100050
    assert meta["fiberid"] == "379"
    assert meta["too_id"] == 61120000000435
    assert meta["programname"] == "TOO_TNS"
    assert meta["ebv"] == pytest.approx(0.0125)
    assert meta["sn_median"] == pytest.approx(6.5)
    for key in ("ra", "dec", "z", "z_err", "ebv", "sn_median"):
        assert type(meta[key]) is float, key
    for key in ("zwarning", "field", "too_id"):
        assert type(meta[key]) is int, key
    for key in ("class", "subclass", "fiberid", "programname"):
        assert type(meta[key]) is str, key


def test_blank_strings_become_none(make_fits):
    meta = read_boss_spectrum(make_fits(spall={"SUBCLASS": "", "PROGRAMNAME": "  "}))[3]
    assert meta["subclass"] is None
    assert meta["programname"] is None


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_non_finite_floats_become_none(make_fits, bad):
    meta = read_boss_spectrum(make_fits(spall={"Z": bad}))[3]
    assert meta["z"] is None


def test_missing_columns_become_none_and_keys_still_present(make_fits):
    meta = read_boss_spectrum(make_fits(drop=["Z", "CLASS", "TOO_ID"]))[3]
    assert meta["z"] is None
    assert meta["class"] is None
    assert meta["too_id"] is None
    assert all(k in meta for k in OPTIONAL_KEYS)
    assert meta["subclass"] == "STARFORMING"  # others unaffected


def test_all_optional_columns_missing(make_fits):
    drop = ["Z", "Z_ERR", "CLASS", "SUBCLASS", "ZWARNING", "FIELD", "FIBERID_LIST",
            "TOO_ID", "PROGRAMNAME", "EBV", "SN_MEDIAN_ALL"]
    meta = read_boss_spectrum(make_fits(drop=drop))[3]
    assert all(meta[k] is None for k in OPTIONAL_KEYS)
    assert meta["ra"] == pytest.approx(150.25)  # required fields still read


# --- z_err sentinel --------------------------------------------------------

@pytest.mark.parametrize("neg", [-1.0, -0.5, -1e-6])
def test_negative_z_err_becomes_none(make_fits, neg):
    assert read_boss_spectrum(make_fits(spall={"Z_ERR": neg}))[3]["z_err"] is None


def test_zero_and_positive_z_err_are_kept(make_fits):
    assert read_boss_spectrum(make_fits(spall={"Z_ERR": 0.0}))[3]["z_err"] == 0.0
    assert read_boss_spectrum(make_fits(spall={"Z_ERR": 0.02}))[3]["z_err"] == pytest.approx(0.02)


def test_negative_z_does_not_become_none(make_fits):
    """Only z_err uses a negative sentinel; blueshifted stars have z < 0."""
    meta = read_boss_spectrum(make_fits(spall={"Z": -0.004}))[3]
    assert meta["z"] == pytest.approx(-0.004)


# --- failure modes ---------------------------------------------------------

def test_missing_spall_extension_raises_value_error(make_fits):
    with pytest.raises(ValueError, match="SPALL"):
        read_boss_spectrum(make_fits(with_spall=False))


@pytest.mark.parametrize("col", ["FIBER_RA", "FIBER_DEC"])
def test_missing_required_coordinate_column_raises_key_error(make_fits, col):
    with pytest.raises(KeyError):
        read_boss_spectrum(make_fits(drop=[col]))


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_boss_spectrum(tmp_path / "nope.fits")


# --- _clean ----------------------------------------------------------------

@pytest.mark.parametrize(
    "value, expected",
    [
        (np.float32(1.5), 1.5),
        (np.int64(7), 7),
        (b"  abc ", "abc"),
        ("  abc ", "abc"),
        ("", None),
        (b"   ", None),
        (float("nan"), None),
        (np.float64("inf"), None),
        (0, 0),
        (0.0, 0.0),
    ],
)
def test_clean(value, expected):
    result = _clean(value)
    assert result == expected
    if expected is not None:
        assert not isinstance(result, np.generic)


# --- real data (skipped if data/ is absent) --------------------------------

@pytest.mark.skipif(not REAL_FILES, reason="data/spectra not present")
@pytest.mark.parametrize("path", REAL_FILES[:5])
def test_real_spectra_load(path):
    wave, flux, ivar, meta = read_boss_spectrum(path)
    assert len(wave) == len(flux) == len(ivar) > 1000
    assert np.all(np.diff(wave) > 0)
    assert np.all(np.isfinite(wave))
    assert 0 <= meta["ra"] < 360
    assert -90 <= meta["dec"] <= 90
    assert meta["z_err"] is None or meta["z_err"] >= 0
