"""Tests for the SNID-SAGE classifier wrapper."""

import dataclasses
from types import SimpleNamespace

import numpy as np
import pytest
from snid_sage.snid.snid import SNIDResult

from bosscat.classifiers import snid as snid_mod
from bosscat.classifiers.snid import make_snid, snid_result_to_fields
from bosscat.config import SnidConfig


def stub_result(**overrides):
    """Build a stub with exactly the real SNIDResult attribute names."""
    values = {f.name: None for f in dataclasses.fields(SNIDResult)}
    values.update(overrides)
    return SimpleNamespace(**values)


class Recorder:
    """Records SNID calls and returns a canned result."""

    def __init__(self, result=None, raises=None):
        self.result = result or stub_result(success=True)
        self.raises = raises
        self.pre_inputs = []
        self.run_calls = []

    def api(self):
        def pre(input_spectrum):
            self.pre_inputs.append(input_spectrum)
            return {"processed": True}, {}

        def run(processed, templates_dir, **kwargs):
            if self.raises:
                raise self.raises
            self.run_calls.append((processed, templates_dir, kwargs))
            return self.result, {}

        return pre, run


def spectrum(n=20):
    wave = 10 ** (3.55 + 1e-4 * np.arange(n))
    return wave, np.ones(n), np.full(n, 4.0)


def build(monkeypatch, tmp_path, rec, **cfg_kwargs):
    monkeypatch.setattr(snid_mod, "_snid_api", rec.api)
    return make_snid(SnidConfig(templates_dir=str(tmp_path), **cfg_kwargs))


def test_flattener_uses_consensus_age_and_both_redshifts():
    res = stub_result(
        success=True, redshift=0.1, redshift_error=0.01,
        consensus_redshift=0.2, consensus_redshift_error=0.02,
        consensus_age=5.0, consensus_age_error=1.5,
        consensus_type="Ia", best_subtype="Ia-norm",
        template_name="sn2004", template_type="Ia", template_subtype="norm",
        lap=0.8, r=9.0, hsigma_lap_ccc=3.0,
    )
    out = snid_result_to_fields(res)
    assert out["snid_redshift"] == 0.1
    assert out["snid_redshift_err"] == 0.01
    assert out["snid_consensus_redshift"] == 0.2
    assert out["snid_consensus_redshift_err"] == 0.02
    assert out["snid_consensus_age"] == 5.0
    assert out["snid_consensus_age_err"] == 1.5
    assert out["snid_consensus_type"] == "Ia"
    assert out["snid_hsigma_lap_ccc"] == 3.0


def test_type_fraction_keys_sanitized():
    res = stub_result(type_fractions={"Ia-norm/91T": 0.6, "II P": 0.4})
    out = snid_result_to_fields(res)
    assert out["snid_frac_Ia-norm_91T"] == 0.6
    assert out["snid_frac_II_P"] == 0.4


def test_numpy_scalars_become_python():
    res = stub_result(
        success=np.bool_(True), redshift=np.float64(0.1),
        type_fractions={"Ia": np.float32(0.5)},
    )
    out = snid_result_to_fields(res)
    assert type(out["snid_success"]) is bool
    assert type(out["snid_redshift"]) is float
    assert type(out["snid_frac_Ia"]) is float
    assert out["snid_consensus_age"] is None


def test_calls_snid_blind_with_plots_off(monkeypatch, tmp_path):
    rec = Recorder()
    clf = build(monkeypatch, tmp_path, rec)
    clf.run(*spectrum(), {"z": 0.5})
    _, templates_dir, kwargs = rec.run_calls[0]
    assert templates_dir == str(tmp_path)
    assert kwargs["show_plots"] is False
    assert kwargs["zmax"] == 0.7
    assert kwargs["zmin"] == -0.01
    assert "forced_redshift" not in kwargs


def test_only_valid_pixels_reach_snid(monkeypatch, tmp_path):
    rec = Recorder()
    clf = build(monkeypatch, tmp_path, rec)
    wave, flux, ivar = spectrum(10)
    flux[2] = np.nan  # ivar > 0
    flux[4] = np.nan
    ivar[4] = 0.0
    ivar[6] = 0.0
    out = clf.run(wave, flux, ivar, {})
    w, f = rec.pre_inputs[0]
    assert len(w) == len(f) == 7
    assert np.all(np.isfinite(f))
    assert out["snid_snr_median"] == pytest.approx(2.0)


def test_rebin_factor_applied(monkeypatch, tmp_path):
    rec = Recorder()
    clf = build(monkeypatch, tmp_path, rec, rebin_factor=2)
    clf.run(*spectrum(20), {})
    assert len(rec.pre_inputs[0][0]) == 10


def test_all_masked_returns_error_without_calling_snid(monkeypatch, tmp_path):
    rec = Recorder()
    clf = build(monkeypatch, tmp_path, rec)
    wave, flux, _ = spectrum(10)
    out = clf.run(wave, flux, np.zeros(10), {})
    assert out == {"snid_success": False, "snid_error": "no valid pixels"}
    assert not rec.pre_inputs and not rec.run_calls


def test_fewer_pixels_than_factor_returns_error(monkeypatch, tmp_path):
    rec = Recorder()
    clf = build(monkeypatch, tmp_path, rec, rebin_factor=8)
    out = clf.run(*spectrum(5), {})
    assert out == {"snid_success": False, "snid_error": "no valid pixels"}
    assert not rec.pre_inputs and not rec.run_calls


def test_short_spectrum_error_propagates(monkeypatch, tmp_path):
    rec = Recorder(raises=ValueError("too short"))
    clf = build(monkeypatch, tmp_path, rec)
    with pytest.raises(ValueError, match="too short"):
        clf.run(*spectrum(), {})


def test_failed_result_keeps_columns(monkeypatch, tmp_path):
    rec = Recorder(result=stub_result(success=False, consensus_age=None))
    clf = build(monkeypatch, tmp_path, rec)
    out = clf.run(*spectrum(), {})
    assert out["snid_success"] is False
    assert out["snid_consensus_age"] is None


def test_provenance(monkeypatch, tmp_path):
    clf = build(monkeypatch, tmp_path, Recorder())
    prov = clf.provenance()
    assert clf.name == "snid"
    assert prov["templates_dir"] == str(tmp_path)
    assert prov["config"]["zmax"] == 0.7
    assert prov["snid_sage_version"]
