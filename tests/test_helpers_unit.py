"""Opt-in in-process unit tests for the skill's pure helpers.

This is the ONLY module that imports ``dnaerys`` (transitively, by importing the
script in-process). It is collected only when ``--run-unit`` is passed, so the
default black-box suite never imports the client.

Run with:
    uv run --with pytest --with dnaerys pytest --run-unit tests/test_helpers_unit.py
"""

from __future__ import annotations

import importlib.util
import types
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

# Load the script as a module (executes its top-level `from dnaerys import ...`).
_SCRIPT = (Path(__file__).parent.parent / "skills" / "onekgpd" / "scripts" / "onekgpd_api.py").resolve()
_spec = importlib.util.spec_from_file_location("onekgpd_api", _SCRIPT)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)

from dnaerys import DnaerysConnectionError, DnaerysInvalidRequestError  # noqa: E402


# --- _call_with_retry -----------------------------------------------------


def test_retry_success_runs_once():
    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        return "ok"

    assert mod._call_with_retry(fn) == "ok"
    assert calls["n"] == 1


def test_retry_retries_retryable_then_succeeds(monkeypatch):
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)
    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        if calls["n"] < mod.MAX_RETRIES:
            raise DnaerysConnectionError("transient")  # is_retryable = True
        return "ok"

    assert mod._call_with_retry(fn) == "ok"
    assert calls["n"] == mod.MAX_RETRIES


def test_retry_nonretryable_propagates_immediately():
    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        raise DnaerysInvalidRequestError("bad")  # is_retryable = False

    with pytest.raises(DnaerysInvalidRequestError):
        mod._call_with_retry(fn)
    assert calls["n"] == 1


def test_retry_exhausts_attempts(monkeypatch):
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)
    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        raise DnaerysConnectionError("always")

    with pytest.raises(DnaerysConnectionError):
        mod._call_with_retry(fn)
    assert calls["n"] == mod.MAX_RETRIES


# --- pure helpers ---------------------------------------------------------


def test_zygosity_mapping():
    assert mod._zygosity(types.SimpleNamespace(het_only=False, hom_only=False)) == (True, True)
    assert mod._zygosity(types.SimpleNamespace(het_only=True, hom_only=False)) == (False, True)
    assert mod._zygosity(types.SimpleNamespace(het_only=False, hom_only=True)) == (True, False)


def test_split_csv():
    assert mod._split_csv("A, B ,C") == ["A", "B", "C"]
    assert mod._split_csv("") == []
    assert mod._split_csv(" , ,X") == ["X"]


def test_parse_region_and_chr_display():
    r = mod._parse_region_str("chrX:100-200")
    assert mod._chr_to_str(r.chr) == "chrX"
    assert (r.start, r.end) == (100, 200)
    assert mod._chr_to_str(mod._parse_region_str("MT:1-2").chr) == "chrMT"


def test_parse_region_bad_format_raises():
    with pytest.raises(ValueError):
        mod._parse_region_str("chr1_1_2")


def test_build_annotation_filter_none_when_empty():
    assert mod._build_annotation_filter(types.SimpleNamespace()) is None


def test_build_annotation_filter_resolves_csv():
    af = mod._build_annotation_filter(types.SimpleNamespace(consequence="MISSENSE_VARIANT"))
    assert af is not None
    assert af.consequence  # resolved to a non-empty tuple of enum members


def test_build_annotation_filter_am_conflict_exits():
    args = types.SimpleNamespace(
        alpha_missense_class="AM_AMBIGUOUS",
        alpha_missense_score_lt=0.5,
    )
    with pytest.raises(SystemExit):
        mod._build_annotation_filter(args)
