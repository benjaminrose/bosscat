"""SNID-SAGE classifier wrapper."""

from __future__ import annotations

import importlib.metadata
import os
from collections.abc import Callable

import numpy as np

from bosscat.classifiers import Classifier
from bosscat.config import SnidConfig
from bosscat.preprocess import rebin_loglam


def _snid_api() -> tuple[Callable, Callable]:
    """Import SNID-SAGE lazily.

    Returns
    -------
    preprocess_spectrum : callable
    run_snid_analysis : callable
    """
    # snid_sage imports pyplot at module load; avoid a GUI backend.
    os.environ.setdefault("MPLBACKEND", "Agg")
    from snid_sage.snid.snid import preprocess_spectrum, run_snid_analysis

    return preprocess_spectrum, run_snid_analysis


def _plain(value):
    """Convert numpy scalars to Python scalars; leave other values alone."""
    if isinstance(value, np.generic):
        return value.item()
    return value


def snid_result_to_fields(result) -> dict:
    """Flatten a ``SNIDResult`` into scalar columns.

    Parameters
    ----------
    result : snid_sage.snid.snid.SNIDResult
        SNID-SAGE result object.

    Returns
    -------
    dict
        ``snid_*`` columns, with plain Python scalars (None kept as None) and
        one ``snid_frac_<type>`` per entry of ``type_fractions``. The
        ``consensus_age``, ``consensus_age_error`` and ``r`` attributes are
        deliberately not emitted: SNID-SAGE 1.2.5 returns constant 0.0.
    """
    mapping = {
        "snid_success": "success",
        "snid_consensus_type": "consensus_type",
        "snid_best_subtype": "best_subtype",
        "snid_template_name": "template_name",
        "snid_template_type": "template_type",
        "snid_template_subtype": "template_subtype",
        "snid_redshift": "redshift",
        "snid_redshift_err": "redshift_error",
        "snid_consensus_redshift": "consensus_redshift",
        "snid_consensus_redshift_err": "consensus_redshift_error",
        "snid_lap": "lap",
        "snid_hsigma_lap_ccc": "hsigma_lap_ccc",
    }
    out = {col: _plain(getattr(result, attr)) for col, attr in mapping.items()}
    seen: dict[str, object] = {}
    for key, frac in (result.type_fractions or {}).items():
        name = str(key).replace("/", "_").replace(" ", "_")
        if name in seen:
            raise ValueError(
                f"type_fractions keys {seen[name]!r} and {key!r} both map to "
                f"column 'snid_frac_{name}'"
            )
        seen[name] = key
        out[f"snid_frac_{name}"] = _plain(frac)
    return out


def _missing(value) -> bool:
    """Return True if ``value`` is None or not a finite number."""
    return value is None or not np.isfinite(value)


NON_SN_TYPES = frozenset({"AGN", "Galaxy", "Star", "CV"})


def snid_flags(result, snr_median: float, cfg: SnidConfig) -> dict[str, bool]:
    """Compute quality flags for a SNID-SAGE result (annotation only).

    Parameters
    ----------
    result : snid_sage.snid.snid.SNIDResult
        SNID-SAGE result object.
    snr_median : float
        Median spectral S/N.
    cfg : SnidConfig
        Thresholds (``min_lap``, ``min_hsigma_lap_ccc``, ``min_snr``,
        ``zmax``, ``zedge_margin``).

    Returns
    -------
    dict
        Boolean flags ``snid_flag_low_quality``, ``snid_flag_low_snr``,
        ``snid_flag_nonsn_template`` and ``snid_flag_zedge``. Missing (None)
        lap, hsigma_lap_ccc, consensus_redshift or snr_median (or non-finite
        values) count as flagged.
    """
    lap, hs, z = result.lap, result.hsigma_lap_ccc, result.consensus_redshift
    return {
        "snid_flag_low_quality": bool(
            _missing(lap)
            or _missing(hs)
            or lap < cfg.min_lap
            or hs < cfg.min_hsigma_lap_ccc
        ),
        "snid_flag_low_snr": bool(_missing(snr_median) or snr_median < cfg.min_snr),
        "snid_flag_nonsn_template": bool(
            result.consensus_type in NON_SN_TYPES
            or result.template_type in NON_SN_TYPES
        ),
        "snid_flag_zedge": bool(_missing(z) or z >= cfg.zmax - cfg.zedge_margin),
    }


def make_snid(cfg: SnidConfig) -> Classifier:
    """Build the SNID-SAGE classifier.

    Classification is blind: the catalogue redshift is never used.

    Parameters
    ----------
    cfg : SnidConfig
        SNID-SAGE settings.

    Returns
    -------
    Classifier
        Classifier named ``"snid"``.
    """

    def run(wave, flux, ivar, meta):
        templates_dir = cfg.resolved_templates_dir()
        w, f, iv = rebin_loglam(wave, flux, ivar, cfg.rebin_factor)
        if w.size == 0:
            return {"snid_success": False, "snid_error": "no valid pixels"}
        preprocess_spectrum, run_snid_analysis = _snid_api()
        processed, _ = preprocess_spectrum(input_spectrum=(w, f))
        result, _ = run_snid_analysis(
            processed,
            templates_dir,
            zmin=cfg.zmin,
            zmax=cfg.zmax,
            show_plots=False,
        )
        fields = snid_result_to_fields(result)
        snr_median = float(np.median(f * np.sqrt(iv)))
        fields["snid_snr_median"] = snr_median
        fields.update(snid_flags(result, snr_median, cfg))
        return fields

    def provenance():
        return {
            "snid_sage_version": importlib.metadata.version("snid-sage"),
            "config": cfg.as_dict(),
            "templates_dir": cfg.resolved_templates_dir(),
        }

    return Classifier(name="snid", run=run, provenance=provenance)
