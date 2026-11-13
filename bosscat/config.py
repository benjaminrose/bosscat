"""Configuration for SNID-SAGE classification."""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass
from pathlib import Path

TEMPLATES_ENV_VAR = "BOSSCAT_SNID_TEMPLATES"


@dataclass(frozen=True)
class SnidConfig:
    """Settings for the SNID-SAGE classifier.

    Parameters
    ----------
    templates_dir : str or None, optional
        Directory of SNID templates. If None, the ``BOSSCAT_SNID_TEMPLATES``
        environment variable is used.
    zmin, zmax : float, optional
        Redshift search range.
    rebin_factor : int, optional
        Rebinning factor applied to spectra before classification.
    min_lap : float, optional
        Minimum template overlap fraction (SNID-SAGE default).
    min_hsigma_lap_ccc : float, optional
        Minimum hsigma-lap-ccc quality metric (SNID-SAGE default).
    min_snr : float, optional
        Minimum spectral S/N. Provisional.
    zedge_margin : float, optional
        Margin from the redshift range edges used for flagging. Provisional.
    """

    templates_dir: str | None = None
    zmin: float = -0.01
    zmax: float = 0.7
    rebin_factor: int = 1
    min_lap: float = 0.3
    min_hsigma_lap_ccc: float = 1.5
    min_snr: float = 3.0
    zedge_margin: float = 0.02

    def resolved_templates_dir(self) -> str:
        """Return the templates directory, validated.

        Returns
        -------
        str
            Path to an existing templates directory. ``templates_dir`` takes
            precedence over the ``BOSSCAT_SNID_TEMPLATES`` environment variable.

        Raises
        ------
        RuntimeError
            If no directory is configured or it does not exist.
        """
        path = self.templates_dir or os.environ.get(TEMPLATES_ENV_VAR)
        if not path:
            raise RuntimeError(
                f"SNID templates directory not set. Pass templates_dir or "
                f"set the {TEMPLATES_ENV_VAR} environment variable."
            )
        if not Path(path).is_dir():
            raise RuntimeError(f"SNID templates directory does not exist: {path}")
        return str(path)

    def as_dict(self) -> dict:
        """Return the configuration as a plain dictionary.

        Returns
        -------
        dict
            Field names mapped to values.
        """
        return asdict(self)
