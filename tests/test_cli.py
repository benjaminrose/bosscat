"""Tests for the bosscat command-line interface and its pipeline wiring."""

import csv
import sys
import pytest

from bosscat import __main__ as cli
from bosscat.classifiers import Classifier
from bosscat.config import SnidConfig


def run_cli(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["bosscat", *argv])
    cli.main()


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# --- argument parsing ------------------------------------------------------

def test_defaults(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits"])
    args = cli.parse_args()
    assert args.input_files == ["a.fits"]
    assert args.classifier == ["snid"]
    assert args.no_classifier is False
    assert args.zmax is None
    assert args.rebin_factor is None
    assert args.templates_dir is None
    assert args.output == "results.csv"


def test_multiple_inputs_and_output_flag(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "b.fits", "-o", "out/x.csv"])
    args = cli.parse_args()
    assert args.input_files == ["a.fits", "b.fits"]
    assert args.output == "out/x.csv"


def test_long_output_flag(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "--output", "y.csv"])
    assert cli.parse_args().output == "y.csv"


def test_classifier_flag_does_not_swallow_positionals(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "--classifier", "snid", "a.fits", "b.fits"])
    args = cli.parse_args()
    assert args.classifier == ["snid"]
    assert args.input_files == ["a.fits", "b.fits"]


def test_unknown_classifier_is_an_error(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "--classifier", "nope"])
    with pytest.raises(SystemExit) as exc:
        cli.parse_args()
    assert exc.value.code == 2


def test_numeric_options_parsed(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "--zmax", "0.5",
                                      "--rebin-factor", "2", "--templates-dir", "/t"])
    args = cli.parse_args()
    assert args.zmax == 0.5
    assert args.rebin_factor == 2
    assert args.templates_dir == "/t"


def test_no_input_files_is_an_error(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["bosscat"])
    with pytest.raises(SystemExit) as exc:
        cli.parse_args()
    assert exc.value.code == 2
    assert "input_files" in capsys.readouterr().err


def test_unknown_option_is_an_error(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "--bogus"])
    with pytest.raises(SystemExit) as exc:
        cli.parse_args()
    assert exc.value.code == 2


# --- main() wiring ---------------------------------------------------------

def _capture_pipeline(monkeypatch):
    calls = []
    monkeypatch.setattr(
        cli, "run_pipeline",
        lambda files, out, classifiers=(): calls.append((files, out, list(classifiers))))
    return calls


def test_main_no_classifier_runs_none(monkeypatch):
    calls = _capture_pipeline(monkeypatch)
    run_cli(monkeypatch, "a.fits", "b.fits", "-o", "o.csv", "--no-classifier")
    assert calls == [(["a.fits", "b.fits"], "o.csv", [])]


def test_main_builds_snid_from_options(monkeypatch):
    calls = _capture_pipeline(monkeypatch)
    seen = []

    def fake_make(cfg):
        seen.append(cfg)
        return Classifier("snid", lambda *a: {}, lambda: {})

    monkeypatch.setitem(cli.CLASSIFIER_FACTORIES, "snid", fake_make)
    run_cli(monkeypatch, "a.fits", "--zmax", "0.3", "--rebin-factor", "2",
            "--templates-dir", "/t")
    assert seen == [SnidConfig(templates_dir="/t", zmax=0.3, rebin_factor=2)]
    assert [c.name for c in calls[0][2]] == ["snid"]


def test_main_unspecified_options_keep_config_defaults(monkeypatch):
    _capture_pipeline(monkeypatch)
    seen = []
    monkeypatch.setitem(
        cli.CLASSIFIER_FACTORIES, "snid",
        lambda cfg: seen.append(cfg) or Classifier("snid", lambda *a: {}, lambda: {}))
    run_cli(monkeypatch, "a.fits")
    assert seen == [SnidConfig()]


def test_main_setup_error_exits_2(monkeypatch, make_fits, tmp_path, capsys):
    monkeypatch.delenv("BOSSCAT_SNID_TEMPLATES", raising=False)
    out = tmp_path / "r.csv"
    with pytest.raises(SystemExit) as exc:
        run_cli(monkeypatch, str(make_fits()), "-o", str(out))
    assert exc.value.code == 2
    assert "templates directory" in capsys.readouterr().err
    assert not out.exists()


# --- end to end (no classifier) --------------------------------------------

def test_end_to_end_no_classifier_writes_expected_csv(monkeypatch, make_fits, tmp_path):
    f1 = make_fits("one.fits", spall={"Z_ERR": -1.0, "CLASS": "STAR"})
    f2 = make_fits("two.fits", drop=["Z", "CLASS"])
    out = tmp_path / "sub" / "results.csv"  # parent dir does not exist yet

    run_cli(monkeypatch, str(f1), str(f2), "-o", str(out), "--no-classifier")

    rows = read_csv(out)
    assert [r["filename"] for r in rows] == ["one.fits", "two.fits"]
    assert rows[0]["boss_class"] == "STAR"
    assert rows[0]["boss_z_err"] == ""       # negative sentinel -> empty cell
    assert float(rows[0]["ra"]) == pytest.approx(150.25)
    assert rows[0]["field"] == "100050"
    assert rows[0]["programname"] == "TOO_TNS"
    assert rows[1]["boss_z"] == ""           # missing column -> empty cell
    assert rows[1]["boss_class"] == ""
    assert not any(k.startswith("snid_") for k in rows[0])
    assert "plate" not in rows[0]


def test_unreadable_file_is_recorded_and_batch_continues(monkeypatch, make_fits,
                                                         tmp_path, capsys):
    good = make_fits("good.fits")
    bad = make_fits("bad.fits", with_spall=False)
    missing = tmp_path / "missing.fits"
    out = tmp_path / "r.csv"

    run_cli(monkeypatch, str(bad), str(missing), str(good), "-o", str(out),
            "--no-classifier")

    rows = {r["filepath"]: r for r in read_csv(out)}
    assert len(rows) == 3
    assert "SPALL" in rows[str(bad)]["error"]
    assert rows[str(missing)]["error"]
    assert rows[str(good)]["filename"] == "good.fits"
    assert rows[str(good)]["error"] == ""
    assert "[ERROR]" in capsys.readouterr().out
