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
        Minimum template overlap fraction (SNID-SAGE default ``lapmin``).
    min_hsigma_lap_ccc : float, optional
        Minimum hsigma-lap-ccc quality metric (SNID-SAGE default threshold).
    min_snr : float, optional
        Minimum median S/N per pixel, measured after rebinning.
    zedge_margin : float, optional
        Distance below ``zmax`` within which a consensus redshift is flagged.

    Notes
    -----
    Defaults calibrated in
    ``docs/superpowers/reports/2026-11-13-snid-calibration.md``. The data
    could not discriminate better values: ``min_lap`` and
    ``min_hsigma_lap_ccc`` stay at SNID-SAGE's own defaults (so the
    low-quality flag marks only failed fits); ``min_snr = 3`` brackets all 13
    SNID failures (S/N 0.96 to 2.82, measured after rebinning at factor 1);
    ``zedge_margin`` is untested (no fit came within 0.4 of ``zmax``).
    ``rebin_factor`` defaults to 1; factor 4 or 8 (which raised NGSF/TNS
    agreement) is pending a user decision, see the report.
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
