"""SNID-SAGE classifier wrapper."""

from __future__ import annotations

import importlib.metadata
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
        one ``snid_frac_<type>`` per entry of ``type_fractions``.
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
        "snid_consensus_age": "consensus_age",
        "snid_consensus_age_err": "consensus_age_error",
        "snid_lap": "lap",
        "snid_r": "r",
        "snid_hsigma_lap_ccc": "hsigma_lap_ccc",
    }
    out = {col: _plain(getattr(result, attr)) for col, attr in mapping.items()}
    for key, frac in (result.type_fractions or {}).items():
        name = str(key).replace("/", "_").replace(" ", "_")
        out[f"snid_frac_{name}"] = _plain(frac)
    return out


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
        fields["snid_snr_median"] = float(np.median(f * np.sqrt(iv)))
        return fields

    def provenance():
        return {
            "snid_sage_version": importlib.metadata.version("snid-sage"),
            "config": cfg.as_dict(),
            "templates_dir": cfg.resolved_templates_dir(),
        }

    return Classifier(name="snid", run=run, provenance=provenance)
