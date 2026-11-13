"""Spectrum preprocessing helpers."""

from __future__ import annotations

import numpy as np

# BOSS loglam is float32, so spacing jitter is ~5e-3 of the step.
_LOG_UNIFORM_TOL = 2e-2


def rebin_loglam(
    wave: np.ndarray, flux: np.ndarray, ivar: np.ndarray, factor: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rebin a log-lambda spectrum by an integer factor, inverse-variance weighted.

    Pixels with ``ivar <= 0`` or non-finite ``flux`` or ``ivar`` are masked.
    Trailing pixels that do not fill a bin are dropped, and bins with no valid
    pixels are omitted. With ``factor == 1`` only the valid pixels are returned.

    Parameters
    ----------
    wave : numpy.ndarray
        Wavelengths (Angstrom), uniformly spaced in log10.
    flux : numpy.ndarray
        Flux values.
    ivar : numpy.ndarray
        Inverse variances.
    factor : int
        Number of consecutive pixels per output bin (>= 1).

    Returns
    -------
    wave : numpy.ndarray
        ``10**mean(log10(wave))`` per bin.
    flux : numpy.ndarray
        ``sum(ivar*flux)/sum(ivar)`` per bin.
    ivar : numpy.ndarray
        ``sum(ivar)`` per bin.

    Raises
    ------
    ValueError
        If ``factor < 1`` or the log10(wave) spacing varies by more than 2%
        of its median.
    """
    if factor < 1:
        raise ValueError(f"factor must be >= 1, got {factor}")
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    ivar = np.asarray(ivar, dtype=float)

    logw = np.log10(wave)
    step = np.diff(logw)
    if step.size and (step.max() - step.min()) / np.median(step) > _LOG_UNIFORM_TOL:
        raise ValueError("wave is not uniformly spaced in log10")

    good = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)
    weight = np.where(good, ivar, 0.0)
    wflux = np.where(good, flux, 0.0) * weight

    n = (len(wave) // factor) * factor
    shape = (n // factor, factor)
    sum_w = weight[:n].reshape(shape).sum(axis=1)
    sum_wf = wflux[:n].reshape(shape).sum(axis=1)
    mean_logw = logw[:n].reshape(shape).mean(axis=1)

    keep = sum_w > 0
    return (
        10 ** mean_logw[keep],
        sum_wf[keep] / sum_w[keep],
        sum_w[keep],
    )
