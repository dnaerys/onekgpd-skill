"""Shared fixtures and helpers for the OneKGPd skill test suite.

Tests drive the skill as a **black box**: each invokes
``uv run <scripts/onekgpd_api.py> <args>`` as a subprocess and asserts on the
exit code, stdout, stderr, and the JSON output file. The suite depends only on
the standard library + pytest and never imports ``dnaerys``.

The sole exception is the opt-in unit tier in ``test_helpers_unit.py`` (which
imports the script in-process and therefore imports ``dnaerys``). It is *not
collected* unless ``--run-unit`` is passed, so the default run preserves the
"tests never import the client" boundary.

Live tests request the session-scoped ``live_server`` fixture, which probes the
public endpoint once: it SKIPS the live tier on a connection-class failure and
FAILS on wrong data (``samples_total != 3202``). An all-skipped live run is made
loud by ``pytest_terminal_summary``.
"""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

# --- Locations ------------------------------------------------------------

SCRIPT = (
    Path(__file__).parent.parent / "skills" / "onekgpd" / "scripts" / "onekgpd_api.py"
).resolve()

# --- Stable test inputs ---------------------------------------------------

# BRCA1 window on GRCh38 (matches the script's own example); re-verify against
# Ensembl if the dataset ever changes.
BRCA1_CHROM = "chr17"
BRCA1_START = 43044292
BRCA1_END = 43170245

# Documented 1000 Genomes Yoruba trio (family Y117); both relationships
# verified live before encoding: child<->parent FIRST_DEGREE, founders UNRELATED.
TRIO_MOTHER = "NA19238"
TRIO_FATHER = "NA19239"
TRIO_CHILD = "NA19240"

EXPECTED_SAMPLES_TOTAL = 3202

# Exact key set of a serialized variant (no cadd_*; gnomAD keys plural).
VARIANT_KEYS = {
    "chr", "start", "end", "ref", "alt", "af", "ac", "an",
    "homc", "hetc", "misc", "homfc", "hetfc", "misfc",
    "gnomad_exomes_af", "gnomad_genomes_af", "am_score", "amino_acids", "biallelic",
}

_SAVED_RE = re.compile(r"\[\*\] Full JSON saved to (.+)$", re.MULTILINE)

# Set by the live_server probe when the endpoint is unreachable, read by the
# terminal-summary hook to make an all-skipped live tier loud.
_LIVE_SKIP = {"reason": None}


# --- Subprocess helper ----------------------------------------------------


@dataclass
class CliResult:
    code: int
    stdout: str
    stderr: str
    json: dict | None
    saved_path: str | None


def run_cli(args, *, output=None, timeout=120) -> CliResult:
    """Run the skill as a subprocess and capture the result.

    *args* is the list of CLI tokens after the script name. If *output* is given
    it is appended as ``--output <path>``; otherwise the JSON file is located by
    parsing the ``[*] Full JSON saved to <path>`` stdout line.
    """
    cmd = ["uv", "run", str(SCRIPT), *args]
    if output is not None:
        cmd += ["--output", str(output)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        pytest.fail(f"`{' '.join(str(a) for a in args)}` timed out after {timeout}s")

    saved_path = None
    m = _SAVED_RE.search(proc.stdout)
    if m:
        saved_path = m.group(1).strip()

    path = str(output) if output is not None else saved_path
    data = None
    if path and Path(path).is_file():
        try:
            data = json.loads(Path(path).read_text())
        except json.JSONDecodeError:
            data = None
    return CliResult(proc.returncode, proc.stdout, proc.stderr, data, saved_path)


@pytest.fixture
def cli(tmp_path):
    """Run the skill with an isolated per-test ``--output`` file by default."""
    counter = {"n": 0}

    def _run(args, *, use_output=True, timeout=120):
        out = None
        if use_output:
            counter["n"] += 1
            out = tmp_path / f"out_{counter['n']}.json"
        return run_cli(args, output=out, timeout=timeout)

    return _run


@pytest.fixture(scope="session")
def live_server():
    """Probe the live endpoint once and gate the live tier.

    - Connection-class failure (non-zero exit) -> SKIP the whole live tier.
    - Connected but wrong data (samples_total != 3202) -> FAIL.
    Returns the parsed ``dataset-info`` JSON for reuse by live tests.
    """
    res = run_cli(["dataset-info"], output=None, timeout=120)
    if res.code != 0:
        last = (res.stderr.strip().splitlines() or ["unknown connection error"])[-1]
        _LIVE_SKIP["reason"] = last
        pytest.skip(f"live 1000 Genomes endpoint unavailable: {last}")
    if res.json is None:
        _LIVE_SKIP["reason"] = "dataset-info returned no parseable JSON"
        pytest.skip(_LIVE_SKIP["reason"])
    total = res.json.get("samples_total")
    if total != EXPECTED_SAMPLES_TOTAL:
        pytest.fail(
            f"live data mismatch: samples_total={total!r}, expected "
            f"{EXPECTED_SAMPLES_TOTAL} — refusing to run against an unexpected dataset"
        )
    return res.json


# --- Opt-in unit tier gating ----------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--run-unit",
        action="store_true",
        default=False,
        help="Also run the opt-in in-process unit tier (imports dnaerys).",
    )


def pytest_ignore_collect(collection_path, config):
    """Never collect (hence never import) the unit tier unless --run-unit is set.

    This keeps the default run free of any ``dnaerys`` import.
    """
    if collection_path.name == "test_helpers_unit.py" and not config.getoption("--run-unit"):
        return True
    return None


# --- Loud skip summary ----------------------------------------------------


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    if _LIVE_SKIP["reason"]:
        skipped = len(terminalreporter.stats.get("skipped", []))
        terminalreporter.write_sep("!", "LIVE TIER SKIPPED", red=True, bold=True)
        terminalreporter.write_line(
            f"1000 Genomes endpoint unreachable: {_LIVE_SKIP['reason']}"
        )
        terminalreporter.write_line(
            f"{skipped} live test(s) skipped; offline contract tests still ran."
        )
