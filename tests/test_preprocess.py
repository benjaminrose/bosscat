"""Tests for inverse-variance log-lambda rebinning."""

import numpy as np
import pytest

from bosscat.preprocess import rebin_loglam


def _grid(n):
    return 10 ** (3.55 + 1e-4 * np.arange(n))


def test_factor_one_keeps_valid_pixels():
    wave = _grid(4)
    flux = np.array([1.0, 2.0, 3.0, 4.0])
    ivar = np.array([1.0, 0.0, 1.0, 1.0])
    w, f, i = rebin_loglam(wave, flux, ivar, 1)
    np.testing.assert_allclose(w, wave[[0, 2, 3]])
    np.testing.assert_allclose(f, [1.0, 3.0, 4.0])
    np.testing.assert_allclose(i, [1.0, 1.0, 1.0])


def test_hand_computed_bin():
    w, f, i = rebin_loglam(_grid(2), np.array([1.0, 3.0]), np.array([1.0, 3.0]), 2)
    np.testing.assert_allclose(f, [2.5])
    np.testing.assert_allclose(i, [4.0])
    np.testing.assert_allclose(w, [10 ** np.mean(np.log10(_grid(2)))])


def test_constant_flux_preserved():
    ivar = np.linspace(0.5, 2.0, 8)
    _, f, _ = rebin_loglam(_grid(8), np.full(8, 7.0), ivar, 4)
    np.testing.assert_allclose(f, [7.0, 7.0])


def test_fully_masked_bin_omitted():
    w, f, i = rebin_loglam(_grid(4), np.ones(4), np.array([0.0, 0.0, 2.0, 2.0]), 2)
    assert len(w) == len(f) == len(i) == 1
    np.testing.assert_allclose(i, [4.0])


def test_partially_masked_bin_uses_valid_pixel():
    _, f, i = rebin_loglam(_grid(2), np.array([5.0, 9.0]), np.array([0.0, 2.0]), 2)
    np.testing.assert_allclose(f, [9.0])
    np.testing.assert_allclose(i, [2.0])


def test_remainder_dropped():
    w, _, _ = rebin_loglam(_grid(5), np.ones(5), np.ones(5), 2)
    assert len(w) == 2


def test_nan_flux_with_positive_ivar_is_masked():
    flux = np.array([np.nan, 4.0])
    _, f, i = rebin_loglam(_grid(2), flux, np.array([1.0, 1.0]), 2)
    np.testing.assert_allclose(f, [4.0])
    np.testing.assert_allclose(i, [1.0])


def test_non_log_uniform_raises():
    wave = np.linspace(4000.0, 9000.0, 20)
    with pytest.raises(ValueError):
        rebin_loglam(wave, np.ones(20), np.ones(20), 2)


def test_float32_boss_like_grid_accepted():
    n = 200
    loglam = np.float32(3.55) + np.arange(n, dtype=np.float32) * np.float32(1e-4)
    wave = 10 ** loglam.astype(float)
    w, _, _ = rebin_loglam(wave, np.ones(n), np.ones(n), 2)
    assert len(w) == n // 2


def test_bad_factor_raises():
    with pytest.raises(ValueError):
        rebin_loglam(_grid(4), np.ones(4), np.ones(4), 0)
