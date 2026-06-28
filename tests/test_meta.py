"""Offline tests for the sample/population metadata script (onekgpd_meta.py).

These are fully deterministic and need no network — they read the bundled
``assets/kgpe.json`` — so they run unconditionally (no live gate). Each test runs
``uv run scripts/onekgpd_meta.py …`` as a subprocess and asserts on the JSON
output file, stdout summary, exit code, and stderr.

Anchors are stable facts of the bundled cohort (3,202 samples, 5 superpopulations,
26 populations) and the documented Yoruba trio (NA19238/39/40). The 1604 F /
1598 M total also guards the HG02300 sex correction.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

META_SCRIPT = (
    Path(__file__).parent.parent / "skills" / "onekgpd" / "scripts" / "onekgpd_meta.py"
).resolve()


@dataclass
class Result:
    code: int
    stdout: str
    stderr: str
    json: dict | None


@pytest.fixture
def meta(tmp_path):
    counter = {"n": 0}

    def _run(args, *, use_output=True, timeout=60):
        cmd = ["uv", "run", str(META_SCRIPT), *args]
        out = None
        if use_output:
            counter["n"] += 1
            out = tmp_path / f"out_{counter['n']}.json"
            cmd += ["--output", str(out)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        data = None
        if out is not None and out.is_file():
            data = json.loads(out.read_text())
        return Result(proc.returncode, proc.stdout, proc.stderr, data)

    return _run


# --- list commands --------------------------------------------------------


def test_list_superpopulations(meta):
    r = meta(["list-superpopulations"])
    assert r.code == 0
    sp = {s["superpopulation_code"]: s for s in r.json["superpopulations"]}
    assert set(sp) == {"AFR", "AMR", "EAS", "EUR", "SAS"}
    assert sp["AFR"]["sample_count"] == 893 and len(sp["AFR"]["populations"]) == 7
    assert sp["AMR"]["sample_count"] == 490 and len(sp["AMR"]["populations"]) == 4
    assert sp["EAS"]["sample_count"] == 585 and len(sp["EAS"]["populations"]) == 5
    assert sp["EUR"]["sample_count"] == 633 and len(sp["EUR"]["populations"]) == 5
    assert sp["SAS"]["sample_count"] == 601 and len(sp["SAS"]["populations"]) == 5
    # ordered by superpopulation code
    assert [s["superpopulation_code"] for s in r.json["superpopulations"]] == \
        ["AFR", "AMR", "EAS", "EUR", "SAS"]


def test_list_populations(meta):
    r = meta(["list-populations"])
    assert r.code == 0
    pops = r.json["populations"]
    assert len(pops) == 26
    assert sum(p["sample_count"] for p in pops) == 3202
    # ordered by (superpopulation, population)
    keys = [(p["superpopulation_code"], p["population_code"]) for p in pops]
    assert keys == sorted(keys)
    assert {p["population_code"] for p in pops} >= {"YRI", "CHS", "GBR", "BEB"}


# --- sample-metadata ------------------------------------------------------


def test_sample_metadata_trio(meta):
    r = meta(["sample-metadata", "--samples", "NA19240,NA19238"])
    assert r.code == 0
    by_id = {s["sample_id"]: s for s in r.json["samples"]}
    # ordered by sample_id
    assert [s["sample_id"] for s in r.json["samples"]] == ["NA19238", "NA19240"]

    child = by_id["NA19240"]
    assert child["paternal_id"] == "NA19239"
    assert child["maternal_id"] == "NA19238"
    assert child["relationship"] == "child"
    assert child["children"] == []
    assert child["population_code"] == "YRI" and child["superpopulation_code"] == "AFR"
    assert child["phase3"] == "FALSE"

    mother = by_id["NA19238"]
    assert mother["paternal_id"] is None and mother["maternal_id"] is None
    assert mother["relationship"] == "mother"
    assert mother["children"] == ["NA19240"]
    assert mother["family_id"] == "Y117"


def test_sample_metadata_hg02300_is_female(meta):
    # Regression guard for the corrected HG02300 sex.
    r = meta(["sample-metadata", "--samples", "HG02300"])
    assert r.code == 0
    assert r.json["samples"][0]["gender"] == "female"


# --- population-stats -----------------------------------------------------


def test_population_stats_yri(meta):
    r = meta(["population-stats", "--populations", "YRI"])
    assert r.code == 0
    p = r.json["populations"][0]
    assert p["population_code"] == "YRI"
    assert p["sample_count"] == 178
    assert p["male_count"] == 97 and p["female_count"] == 81
    assert p["phase3_count"] == 108
    assert p["trio_count"] == 56


def test_population_stats_full_name_with_comma(meta):
    # Full name contains a comma; repeatable flag (not CSV) must keep it intact.
    r = meta(["population-stats", "--populations", "Yoruba in Ibadan, Nigeria"])
    assert r.code == 0
    assert len(r.json["populations"]) == 1
    assert r.json["populations"][0]["population_code"] == "YRI"


def test_population_stats_case_insensitive_and_dedup(meta):
    # lowercase code + full name of the same population collapse to one group.
    r = meta([
        "population-stats",
        "--populations", "yri",
        "--populations", "Yoruba in Ibadan, Nigeria",
    ])
    assert r.code == 0
    assert len(r.json["populations"]) == 1
    assert r.json["populations"][0]["sample_count"] == 178


# --- superpopulation-summary ----------------------------------------------


def test_superpopulation_summary_eas(meta):
    r = meta(["superpopulation-summary", "--superpopulations", "EAS"])
    assert r.code == 0
    s = r.json["superpopulations"][0]
    assert s["superpopulation_code"] == "EAS"
    assert s["sample_count"] == 585
    assert s["male_count"] == 292 and s["female_count"] == 293
    assert s["phase3_count"] == 504 and s["trio_count"] == 72
    assert len(s["populations"]) == 5
    # summary totals equal the sum of the nested per-population breakdown
    assert s["sample_count"] == sum(p["sample_count"] for p in s["populations"])
    assert s["male_count"] == sum(p["male_count"] for p in s["populations"])


def test_superpopulation_summary_totals_match_cohort(meta):
    # All 5 superpopulations must sum to the corrected 3202 / 1598 M / 1604 F.
    r = meta([
        "superpopulation-summary",
        "--superpopulations", "AFR", "--superpopulations", "AMR",
        "--superpopulations", "EAS", "--superpopulations", "EUR",
        "--superpopulations", "SAS",
    ])
    assert r.code == 0
    sp = r.json["superpopulations"]
    assert sum(s["sample_count"] for s in sp) == 3202
    assert sum(s["male_count"] for s in sp) == 1598
    assert sum(s["female_count"] for s in sp) == 1604


# --- select-samples-by-population -----------------------------------------


def test_select_by_population(meta):
    r = meta(["select-samples-by-population", "--population", "YRI", "--limit", "5"])
    assert r.code == 0
    assert r.json["count"] == len(r.json["samples"]) == 5
    assert r.json["samples"] == sorted(r.json["samples"])  # ordered ascending
    assert r.json["request"] == {"population": "YRI", "superpopulation": None,
                                 "skip": 0, "limit": 5}


def test_select_skip_offsets(meta):
    a = meta(["select-samples-by-population", "--population", "YRI", "--limit", "10"])
    b = meta(["select-samples-by-population", "--population", "YRI", "--skip", "5", "--limit", "5"])
    assert a.code == 0 and b.code == 0
    # deterministic ordering: skip=5 starts at a's 6th element
    assert b.json["samples"] == a.json["samples"][5:10]


def test_select_intersection(meta):
    # population YRI is within AFR -> intersect is non-empty; within EUR -> empty.
    inside = meta(["select-samples-by-population", "--population", "YRI",
                   "--superpopulation", "AFR", "--limit", "3202"])
    outside = meta(["select-samples-by-population", "--population", "YRI",
                    "--superpopulation", "EUR", "--limit", "3202"])
    assert inside.json["count"] == 178
    assert outside.json["count"] == 0


def test_select_case_insensitive_full_name(meta):
    r = meta(["select-samples-by-population",
              "--population", "southern han chinese, china", "--limit", "3"])
    assert r.code == 0
    assert r.json["count"] == 3


# --- error paths ----------------------------------------------------------


def test_err_unknown_sample(meta):
    r = meta(["sample-metadata", "--samples", "NA19240,BOGUS_ID"], use_output=False)
    assert r.code == 1
    assert "Unknown sample IDs: [BOGUS_ID]" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_sample_id_case_sensitive(meta):
    r = meta(["sample-metadata", "--samples", "na19240"], use_output=False)
    assert r.code == 1
    assert "Unknown sample IDs: [na19240]" in r.stderr


def test_err_unknown_population(meta):
    r = meta(["population-stats", "--populations", "NOPE"], use_output=False)
    assert r.code == 1
    assert "Unrecognised population values: [NOPE]" in r.stderr


def test_err_neither_pop_nor_superpop(meta):
    r = meta(["select-samples-by-population"], use_output=False)
    assert r.code == 1
    assert "At least one parameter" in r.stderr


def test_err_skip_negative(meta):
    r = meta(["select-samples-by-population", "--population", "YRI", "--skip", "-1"], use_output=False)
    assert r.code == 1
    assert "'skip' must be >= 0" in r.stderr


def test_err_limit_out_of_range(meta):
    r = meta(["select-samples-by-population", "--population", "YRI", "--limit", "9999"], use_output=False)
    assert r.code == 1
    assert "'limit' must be between 1 and 3202" in r.stderr


def test_err_missing_required_flag(meta):
    r = meta(["population-stats"], use_output=False)
    assert r.code == 2  # argparse: --populations required
