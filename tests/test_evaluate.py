"""Tests for the SNID vs NGSF/TNS evaluation module (synthetic CSVs only)."""

import csv
import json
import sys

import pytest

from bosscat import evaluate as ev

MANIFEST_COLS = ["name", "file_name", "tns_type", "sn_family", "z_fit"]


def write_csv(path, rows, cols=None):
    cols = cols or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(rows)
    return path


def man(file_name, tns="", fam="Ia-norm", z="0.1"):
    return {"name": "AT" + file_name, "file_name": file_name,
            "tns_type": tns, "sn_family": fam, "z_fit": z}


def res(file_name, success="True", cz="0.1", z="0.1", typ="Ia", **extra):
    row = {"filename": file_name + ".fits", "snid_success": success,
           "snid_consensus_redshift": cz, "snid_redshift": z,
           "snid_consensus_type": typ}
    row.update(extra)
    return row


def test_broad_type_mapping():
    assert ev.broad_type("SN Ia-91T-like") == "Ia"
    assert ev.broad_type("IIn") == "II"
    assert ev.broad_type("SLSN-I") == "SLSN"
    assert ev.broad_type("") is None
    assert ev.broad_type(None) is None
    assert ev.broad_type("Ia 99aa-like") == "Ia"
    assert ev.broad_type("Ib-norm") == "Ib"
    assert ev.broad_type("Ic-BL") == "Ic"
    assert ev.broad_type("SN IIb") == "II"
    assert ev.broad_type("IIP") == "II"
    assert ev.broad_type("TDE") is None
    assert ev.broad_type("Galaxy") is None


def test_join_on_filename_suffix(tmp_path):
    m = write_csv(tmp_path / "m.csv", [man("a"), man("b"), man("c")])
    r = write_csv(tmp_path / "r.csv", [res("a"), res("b"), res("zzz")])
    out = ev.evaluate(r, m)
    assert out["n_results"] == 3
    assert out["n_matched"] == 2
    assert out["n_unmatched_results"] == 1
    assert out["n_unmatched_manifest"] == 1
    json.dumps(out)


def test_residual_stats_overall_and_high_z(tmp_path):
    m = write_csv(tmp_path / "m.csv", [
        man("a", z="0.1"), man("b", z="0.2"), man("c", z="0.5"),
        man("d", z="0.6"), man("e", z="0.7")])
    r = write_csv(tmp_path / "r.csv", [
        res("a", cz="0.11", z="0.1"), res("b", cz="0.18", z="0.2"),
        res("c", cz="0.52", z="0.5"), res("d", cz="0.64", z="0.6"),
        res("e", cz="", z="0.7")])
    out = ev.evaluate(r, m, high_z=0.3)
    allz = out["redshift"]["snid_consensus_redshift"]["all"]
    assert allz["n"] == 4
    # residuals: 0.01, -0.02, 0.02, 0.04 -> median 0.015
    assert allz["median"] == pytest.approx(0.015)
    # |dev|: 0.005, 0.035, 0.005, 0.025 -> median 0.015
    assert allz["mad"] == pytest.approx(0.015)
    hi = out["redshift"]["snid_consensus_redshift"]["high_z"]
    assert hi["n"] == 2
    assert hi["median"] == pytest.approx(0.03)
    assert hi["mad"] == pytest.approx(0.01)
    best = out["redshift"]["snid_redshift"]
    assert best["all"]["n"] == 5
    assert best["all"]["median"] == pytest.approx(0.0)


def test_agreement_denominators_exclude_blank_tns(tmp_path):
    m = write_csv(tmp_path / "m.csv", [
        man("a", tns="SN Ia", fam="Ia-norm"),
        man("b", tns="", fam="IIn"),
        man("c", tns="SN II", fam="IIP"),
        man("d", tns="Other", fam="Ib")])
    r = write_csv(tmp_path / "r.csv", [
        res("a", typ="Ia"), res("b", typ="II"), res("c", typ="Ia"),
        res("d", typ="Ib")])
    out = ev.evaluate(r, m)
    ngsf = out["type_agreement"]["vs_ngsf"]
    assert (ngsf["n"], ngsf["agree"]) == (4, 3)
    assert ngsf["fraction"] == pytest.approx(0.75)
    tns = out["type_agreement"]["vs_tns"]
    # b blank, d unmappable: only a and c count
    assert (tns["n"], tns["agree"]) == (2, 1)
    assert tns["fraction"] == pytest.approx(0.5)


def test_failed_snid_rows_counted_not_scored(tmp_path):
    m = write_csv(tmp_path / "m.csv", [man("a"), man("b"), man("c"), man("d")])
    r = write_csv(tmp_path / "r.csv", [
        res("a", snid_flag_zedge="True"),
        res("b", success="False", cz="0.9", snid_flag_zedge="False"),
        res("c", success="", cz="", typ=""),
        {"filename": "d.fits"}], )
    out = ev.evaluate(r, m)
    assert out["n_snid_success"] == 1
    assert out["n_matched"] == 4
    assert out["redshift"]["snid_consensus_redshift"]["all"]["n"] == 1
    assert out["type_agreement"]["vs_ngsf"]["n"] == 1
    assert out["n_flagged"] == {"snid_flag_zedge": 1}


def test_empty_stats_are_none(tmp_path):
    m = write_csv(tmp_path / "m.csv", [man("a")])
    r = write_csv(tmp_path / "r.csv", [res("a", success="False")])
    out = ev.evaluate(r, m)
    s = out["redshift"]["snid_consensus_redshift"]["high_z"]
    assert s == {"n": 0, "median": None, "mad": None}
    assert out["type_agreement"]["vs_tns"]["fraction"] is None


def test_cli_text_and_json(tmp_path, monkeypatch, capsys):
    m = write_csv(tmp_path / "m.csv", [man("a", tns="Ia")])
    r = write_csv(tmp_path / "r.csv", [res("a")])
    monkeypatch.setattr(sys, "argv", ["evaluate", str(r), "--manifest", str(m)])
    ev.main()
    assert "n_matched" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv",
                        ["evaluate", str(r), "--manifest", str(m), "--json"])
    ev.main()
    assert json.loads(capsys.readouterr().out)["n_matched"] == 1


def test_snid_type_normalized_to_broad(tmp_path):
    m = write_csv(tmp_path / "m.csv", [
        man("a", fam="Ia-norm"), man("b", fam="II"), man("c", fam="Ib"),
        man("d", fam="Ia"), man("e", fam="Ia")])
    r = write_csv(tmp_path / "r.csv", [
        res("a", typ="Ia-norm"), res("b", typ="IIb"), res("c", typ="Ibn"),
        res("d", typ="AGN"), res("e", typ="")])
    ngsf = ev.evaluate(r, m)["type_agreement"]["vs_ngsf"]
    assert (ngsf["n"], ngsf["agree"]) == (4, 3)


def test_failure_count(tmp_path):
    m = write_csv(tmp_path / "m.csv", [man("a"), man("b")])
    r = write_csv(tmp_path / "r.csv", [res("a"), res("b", success="False")])
    out = ev.evaluate(r, m)
    assert out["n_snid_success"] == 1
    assert out["n_snid_failed"] == 1
    assert "failures" in ev.format_report(out)
