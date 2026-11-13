"""Tests against the real SNID-SAGE library.

Only ``test_real_snid_end_to_end`` needs templates and is marked ``snid_real``.
"""

import dataclasses
import inspect
import time
from pathlib import Path

import numpy as np
import pytest
from snid_sage.snid.snid import SNIDResult, preprocess_spectrum, run_snid_analysis

from bosscat.classifiers.snid import make_snid
from bosscat.config import SnidConfig
from bosscat.io.boss import read_boss_spectrum

SPECTRA_DIR = Path(__file__).resolve().parent.parent / "data" / "spectra"

# Attributes of SNIDResult read by snid_result_to_fields / snid_flags.
RESULT_ATTRS = [
    "success", "consensus_type", "best_subtype", "template_name",
    "template_type", "template_subtype", "redshift", "redshift_error",
    "consensus_redshift", "consensus_redshift_error", "consensus_age",
    "consensus_age_error", "lap", "r", "hsigma_lap_ccc", "type_fractions",
]

EXPECTED_KEYS = [
    "snid_success", "snid_consensus_type", "snid_best_subtype",
    "snid_template_name", "snid_template_type", "snid_template_subtype",
    "snid_redshift", "snid_redshift_err", "snid_consensus_redshift",
    "snid_consensus_redshift_err", "snid_consensus_age",
    "snid_consensus_age_err", "snid_lap", "snid_r", "snid_hsigma_lap_ccc",
    "snid_snr_median", "snid_flag_low_quality", "snid_flag_low_snr",
    "snid_flag_nonsn_template", "snid_flag_zedge",
]


def test_snid_api_signatures():
    """The SNID-SAGE calls made by the wrapper match the installed library."""
    assert "input_spectrum" in inspect.signature(preprocess_spectrum).parameters

    params = list(inspect.signature(run_snid_analysis).parameters)
    assert params[:2] == ["processed_spectrum", "templates_dir"]
    for name in ("zmin", "zmax", "show_plots"):
        assert name in params


def test_snid_result_fields_exist():
    """Every SNIDResult attribute read by the wrapper is a real field."""
    names = {f.name for f in dataclasses.fields(SNIDResult)}
    missing = [a for a in RESULT_ATTRS if a not in names]
    assert not missing


@pytest.mark.snid_real
def test_real_snid_end_to_end(capsys):
    """Run real SNID-SAGE on the first spectrum in ``data/spectra/``."""
    files = sorted(SPECTRA_DIR.glob("*.fits"))
    if not files:
        pytest.skip(f"no spectra found in {SPECTRA_DIR}")
    wave, flux, ivar, meta = read_boss_spectrum(files[0])

    t0 = time.perf_counter()
    out = make_snid(SnidConfig()).run(wave, flux, ivar, meta)
    elapsed = time.perf_counter() - t0

    with capsys.disabled():
        print(f"\n{files[0].name}: {elapsed:.1f} s")
        for key, value in out.items():
            print(f"  {key} = {value!r}")

    assert out["snid_success"] is True
    for key in EXPECTED_KEYS:
        assert key in out, key
    z = out["snid_consensus_redshift"]
    assert isinstance(z, float)
    assert -0.01 <= z <= 0.7
