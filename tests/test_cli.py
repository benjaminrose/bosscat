"""Tests for the bosscat command-line interface and its pipeline wiring."""

import csv
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from bosscat import __main__ as cli
from bosscat import pipeline


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
    assert args.snid is True
    assert args.output == "results.csv"


def test_multiple_inputs_and_output_flag(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "b.fits", "-o", "out/x.csv"])
    args = cli.parse_args()
    assert args.input_files == ["a.fits", "b.fits"]
    assert args.output == "out/x.csv"


def test_long_output_flag(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", "--output", "y.csv"])
    assert cli.parse_args().output == "y.csv"


@pytest.mark.parametrize("flag, expected", [("--snid", True), ("--no-snid", False)])
def test_snid_flag(monkeypatch, flag, expected):
    monkeypatch.setattr(sys, "argv", ["bosscat", "a.fits", flag])
    assert cli.parse_args().snid is expected


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

def test_main_passes_arguments_to_run_pipeline(monkeypatch):
    calls = []
    monkeypatch.setattr(cli, "run_pipeline",
                        lambda files, out, run_snid: calls.append((files, out, run_snid)))
    run_cli(monkeypatch, "a.fits", "b.fits", "-o", "o.csv", "--no-snid")
    assert calls == [(["a.fits", "b.fits"], "o.csv", False)]


# --- end to end (SNID disabled or stubbed) ---------------------------------

def test_end_to_end_no_snid_writes_expected_csv(monkeypatch, make_fits, tmp_path):
    f1 = make_fits("one.fits", spall={"Z_ERR": -1.0, "CLASS": "STAR"})
    f2 = make_fits("two.fits", drop=["Z", "CLASS"])
    out = tmp_path / "sub" / "results.csv"  # parent dir does not exist yet

    run_cli(monkeypatch, str(f1), str(f2), "-o", str(out), "--no-snid")

    rows = read_csv(out)
    assert [r["filename"] for r in rows] == ["one.fits", "two.fits"]
    assert rows[0]["boss_class"] == "STAR"
    assert rows[0]["boss_z_err"] == ""       # negative sentinel -> empty cell
    assert float(rows[0]["ra"]) == pytest.approx(150.25)
    assert rows[0]["field"] == "100050"
    assert rows[0]["programname"] == "TOO_TNS"
    assert rows[1]["boss_z"] == ""           # missing column -> empty cell
    assert rows[1]["boss_class"] == ""
    assert rows[0]["snid_success"] == "False"  # classifier did not run
    assert "plate" not in rows[0]


def test_unreadable_file_is_recorded_and_batch_continues(monkeypatch, make_fits,
                                                         tmp_path, capsys):
    good = make_fits("good.fits")
    bad = make_fits("bad.fits", with_spall=False)
    missing = tmp_path / "missing.fits"
    out = tmp_path / "r.csv"

    run_cli(monkeypatch, str(bad), str(missing), str(good), "-o", str(out), "--no-snid")

    rows = {r["filepath"]: r for r in read_csv(out)}
    assert len(rows) == 3
    assert "SPALL" in rows[str(bad)]["error"]
    assert rows[str(missing)]["error"]
    assert rows[str(good)]["filename"] == "good.fits"
    assert rows[str(good)]["error"] == ""
    assert "[ERROR]" in capsys.readouterr().out


def test_masked_pixels_are_dropped_before_classification(monkeypatch, make_fits, tmp_path):
    seen = {}

    def fake_classify(wave, flux):
        seen["wave"], seen["flux"] = wave, flux
        result = SimpleNamespace(
            success=True, consensus_type="Ia", redshift=0.1, redshift_error=0.01,
            age=3.0, age_error=1.0, lap=5.0, type_fractions={"Ia/norm": 0.9},
        )
        return result, {}

    monkeypatch.setattr(pipeline, "classify_spectrum", fake_classify)
    out = tmp_path / "r.csv"
    run_cli(monkeypatch, str(make_fits()), "-o", str(out))

    # default fixture has ivar == [1, 2, 0, 4, 5]: pixel 2 must be removed
    np.testing.assert_allclose(seen["flux"], [1, 2, 4, 5])
    assert len(seen["wave"]) == 4
    row = read_csv(out)[0]
    assert row["snid_success"] == "True"
    assert row["snid_consensus_type"] == "Ia"
    assert row["snid_frac_Ia_norm"] == "0.9"


def test_classifier_failure_is_recorded_not_raised(monkeypatch, make_fits, tmp_path):
    def boom(wave, flux):
        raise RuntimeError("snid exploded")

    monkeypatch.setattr(pipeline, "classify_spectrum", boom)
    out = tmp_path / "r.csv"
    run_cli(monkeypatch, str(make_fits()), "-o", str(out))

    row = read_csv(out)[0]
    assert row["snid_error"] == "snid exploded"
    assert row["snid_success"] == "False"
    assert row["boss_class"] == "GALAXY"  # metadata still written
