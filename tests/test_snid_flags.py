"""Tests for SNID quality flags."""

import pytest

from bosscat.classifiers import snid as snid_mod
from bosscat.classifiers.snid import make_snid, snid_flags
from bosscat.config import SnidConfig

from test_snid_classifier import Recorder, spectrum, stub_result

CFG = SnidConfig()
KEYS = {
    "snid_flag_low_quality",
    "snid_flag_low_snr",
    "snid_flag_nonsn_template",
    "snid_flag_zedge",
}


def good(**overrides):
    values = dict(
        success=True, lap=0.8, hsigma_lap_ccc=3.0, consensus_redshift=0.3,
        consensus_type="Ia", template_type="Ia",
    )
    values.update(overrides)
    return stub_result(**values)


def test_keys_and_bool_types():
    flags = snid_flags(good(), 10.0, CFG)
    assert set(flags) == KEYS
    assert all(type(v) is bool for v in flags.values())
    assert not any(flags.values())


@pytest.mark.parametrize(
    "overrides, expected",
    [
        ({"lap": 0.1}, True),
        ({"lap": 0.3}, False),
        ({"hsigma_lap_ccc": 1.0}, True),
        ({"hsigma_lap_ccc": 1.5}, False),
    ],
)
def test_low_quality(overrides, expected):
    assert snid_flags(good(**overrides), 10.0, CFG)["snid_flag_low_quality"] is expected


@pytest.mark.parametrize("snr, expected", [(2.9, True), (3.0, False), (50.0, False)])
def test_low_snr(snr, expected):
    assert snid_flags(good(), snr, CFG)["snid_flag_low_snr"] is expected


@pytest.mark.parametrize(
    "overrides, expected",
    [
        ({"consensus_type": "AGN"}, True),
        ({"consensus_type": "Galaxy"}, True),
        ({"template_type": "Star"}, True),
        ({"template_type": "CV"}, True),
        ({"consensus_type": "Ia", "template_type": "II"}, False),
        ({"consensus_type": None, "template_type": None}, False),
    ],
)
def test_nonsn_template(overrides, expected):
    assert snid_flags(good(**overrides), 10.0, CFG)["snid_flag_nonsn_template"] is expected


@pytest.mark.parametrize(
    "z, expected",
    [(0.69, True), (0.68, True), (0.6799, False), (0.3, False)],
)
def test_zedge(z, expected):
    flags = snid_flags(good(consensus_redshift=z), 10.0, CFG)
    assert flags["snid_flag_zedge"] is expected


def test_none_fields_flag_true():
    res = good(lap=None, hsigma_lap_ccc=None, consensus_redshift=None)
    flags = snid_flags(res, 10.0, CFG)
    assert flags["snid_flag_low_quality"] is True
    assert flags["snid_flag_zedge"] is True
    res = good(hsigma_lap_ccc=None)
    assert snid_flags(res, 10.0, CFG)["snid_flag_low_quality"] is True


def run_clf(monkeypatch, tmp_path, result):
    rec = Recorder(result=result)
    monkeypatch.setattr(snid_mod, "_snid_api", rec.api)
    clf = make_snid(SnidConfig(templates_dir=str(tmp_path)))
    return clf.run(*spectrum(), {})


def test_run_merges_flags(monkeypatch, tmp_path):
    out = run_clf(monkeypatch, tmp_path, good())
    assert KEYS <= set(out)
    assert out["snid_flag_low_snr"] is True  # median S/N of spectrum() is 2.0
    assert out["snid_flag_low_quality"] is False
    assert out["snid_consensus_type"] == "Ia"


def test_error_path_has_no_flags(monkeypatch, tmp_path):
    rec = Recorder()
    monkeypatch.setattr(snid_mod, "_snid_api", rec.api)
    clf = make_snid(SnidConfig(templates_dir=str(tmp_path)))
    wave, flux, _ = spectrum(10)
    out = clf.run(wave, flux, wave * 0, {})
    assert out == {"snid_success": False, "snid_error": "no valid pixels"}


def test_flags_do_not_change_type(monkeypatch, tmp_path):
    quiet = run_clf(monkeypatch, tmp_path, good(consensus_type="Ia"))
    loud = run_clf(
        monkeypatch, tmp_path,
        good(consensus_type="Ia", lap=None, hsigma_lap_ccc=None,
             consensus_redshift=0.7, template_type="AGN"),
    )
    assert not quiet["snid_flag_nonsn_template"]
    assert all(loud[k] for k in ("snid_flag_low_quality", "snid_flag_zedge",
                                 "snid_flag_nonsn_template"))
    assert quiet["snid_consensus_type"] == loud["snid_consensus_type"] == "Ia"
