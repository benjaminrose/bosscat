"""Tests for run_pipeline, base_row, and write_provenance."""

import csv
import json
import subprocess
import sys

import pytest

from bosscat.classifiers import Classifier
from bosscat.io.results import base_row, write_provenance, write_results
from bosscat.pipeline import run_pipeline


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fake(name, fields=None, exc=None, prov=None, record=None):
    def run(*args):
        if record is not None:
            record.append(args)
        if exc is not None:
            raise exc
        return dict(fields or {})

    return Classifier(name, run, lambda: dict(prov or {}))


def test_classifier_columns_merged(make_fits, tmp_path):
    out = tmp_path / "r.csv"
    run_pipeline([str(make_fits())], out, [fake("fake", {"fake_a": 1})])
    row = read_csv(out)[0]
    assert row["fake_a"] == "1"
    assert row["boss_class"] == "GALAXY"


def test_one_classifier_fails_other_survives(make_fits, tmp_path):
    out = tmp_path / "r.csv"
    clfs = [fake("bad", exc=RuntimeError("kaboom")), fake("ok", {"ok_x": 7})]
    run_pipeline([str(make_fits("a.fits")), str(make_fits("b.fits"))], out, clfs)
    rows = read_csv(out)
    assert len(rows) == 2
    for row in rows:
        assert row["bad_error"] == "kaboom"
        assert row["bad_success"] == "False"
        assert row["ok_x"] == "7"
        assert row["boss_class"] == "GALAXY"


def test_read_error_row_kept(make_fits, tmp_path):
    out = tmp_path / "r.csv"
    bad = make_fits("bad.fits", with_spall=False)
    record = []
    run_pipeline([str(bad), str(make_fits("good.fits"))], out,
                 [fake("fake", {"fake_a": 1}, record=record)])
    rows = {r["filepath"]: r for r in read_csv(out)}
    assert "SPALL" in rows[str(bad)]["error"]
    assert rows[str(bad)]["fake_a"] == ""
    assert len(record) == 1  # classifier not run on the unreadable file


def test_classifier_receives_unmasked_arrays_and_meta(make_fits, tmp_path):
    record = []
    run_pipeline([str(make_fits())], tmp_path / "r.csv", [fake("fake", record=record)])
    (args,) = record
    assert len(args) == 4
    wave, flux, ivar, meta = args
    assert len(wave) == len(flux) == len(ivar) == 5  # ivar == 0 pixel kept
    assert meta["z"] == pytest.approx(0.5)


def test_provenance_called_once_before_files_and_error_propagates(make_fits, tmp_path):
    def bad_prov():
        raise RuntimeError("no templates")

    run_calls = []
    clf = Classifier("x", lambda *a: run_calls.append(a) or {}, bad_prov)
    out = tmp_path / "r.csv"
    with pytest.raises(RuntimeError, match="no templates"):
        run_pipeline([str(make_fits())], out, [clf])
    assert run_calls == []
    assert not out.exists()


def test_provenance_called_once(make_fits, tmp_path):
    count = []
    clf = Classifier("x", lambda *a: {}, lambda: count.append(1) or {})
    run_pipeline([str(make_fits("a.fits")), str(make_fits("b.fits"))],
                 tmp_path / "r.csv", [clf])
    assert len(count) == 1


def test_sidecar_json_written(make_fits, tmp_path):
    out = tmp_path / "res" / "r.csv"
    f = make_fits()
    run_pipeline([str(f)], out, [fake("fake", prov={"v": "1.2"})])
    meta = json.loads((tmp_path / "res" / "r.json").read_text())
    assert set(meta) == {"bosscat_version", "python_version", "created",
                         "input_files", "classifiers"}
    assert meta["input_files"] == [str(f)]
    assert meta["classifiers"] == {"fake": {"v": "1.2"}}
    assert meta["created"].endswith(("Z", "+00:00"))


def test_write_provenance_returns_path(tmp_path):
    path = write_provenance(tmp_path / "x.csv", [], ["a"])
    assert path == tmp_path / "x.json"
    assert json.loads(path.read_text())["classifiers"] == {}


def test_write_results_unions_differing_keys(tmp_path):
    out = tmp_path / "r.csv"
    write_results([{"a": 1, "snid_frac_Ia": 0.5},
                   {"a": 2, "snid_success": False, "snid_error": "x"}], out)
    rows = read_csv(out)
    assert set(rows[0]) == {"a", "snid_frac_Ia", "snid_success", "snid_error"}
    assert rows[0]["snid_error"] == "" and rows[1]["snid_frac_Ia"] == ""


def test_base_row_has_no_snid_columns():
    row = base_row({"filename": "f", "z": 0.1, "class": "QSO"})
    assert row["boss_z"] == 0.1 and row["boss_class"] == "QSO"
    assert not any(k.startswith("snid_") for k in row)


def test_importing_pipeline_and_cli_does_not_import_snid_sage():
    code = ("import sys, bosscat.pipeline, bosscat.__main__; "
            "sys.exit('snid_sage' in sys.modules)")
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0
