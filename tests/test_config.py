"""Tests for bosscat.config.SnidConfig."""

import pytest

from bosscat.config import SnidConfig


def test_templates_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("BOSSCAT_SNID_TEMPLATES", str(tmp_path))
    assert SnidConfig().resolved_templates_dir() == str(tmp_path)


def test_explicit_overrides_env(monkeypatch, tmp_path):
    env_dir = tmp_path / "env"
    explicit = tmp_path / "explicit"
    env_dir.mkdir()
    explicit.mkdir()
    monkeypatch.setenv("BOSSCAT_SNID_TEMPLATES", str(env_dir))
    cfg = SnidConfig(templates_dir=str(explicit))
    assert cfg.resolved_templates_dir() == str(explicit)


def test_missing_templates_raises(monkeypatch):
    monkeypatch.delenv("BOSSCAT_SNID_TEMPLATES", raising=False)
    with pytest.raises(RuntimeError) as exc:
        SnidConfig(templates_dir=None).resolved_templates_dir()
    msg = str(exc.value)
    assert "BOSSCAT_SNID_TEMPLATES" in msg
    assert "SnidConfig(templates_dir=" in msg
    assert "--templates-dir" in msg


def test_nonexistent_dir_raises(tmp_path):
    cfg = SnidConfig(templates_dir=str(tmp_path / "nope"))
    with pytest.raises(RuntimeError):
        cfg.resolved_templates_dir()


def test_defaults():
    cfg = SnidConfig()
    assert cfg.zmax == 0.7
    assert cfg.zmin == -0.01
    assert cfg.rebin_factor == 1
    assert cfg.as_dict()["min_lap"] == 0.3
