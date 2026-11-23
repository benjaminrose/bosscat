"""Classifier wrappers sharing a common interface."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Classifier:
    """A named classifier with a run function and provenance reporter.

    Parameters
    ----------
    name : str
        Short classifier name, used as a column prefix.
    run : callable
        ``run(wave, flux, ivar, meta) -> dict`` returning flat result fields.
    provenance : callable
        ``provenance() -> dict`` describing versions and settings.
    """

    name: str
    run: Callable[[np.ndarray, np.ndarray, np.ndarray, dict], dict]
    provenance: Callable[[], dict]
