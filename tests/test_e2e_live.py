"""Live end-to-end tests against the public 1000 Genomes endpoint.

Every test requests the session-scoped ``live_server`` fixture, so the whole
module is skipped (loudly) when the endpoint is unreachable. Assertions favor
**stable anchors** (samples_total=3202, a fixed trio pedigree) and
**invariants/monotonicity** over volatile exact counts, so the suite stays green
across dataset reloads.

Anchors verified live before encoding:
- mother<->child / father<->child -> FIRST_DEGREE; founders mother<->father -> UNRELATED
- count-variants == select-variants full-walk count on the small window
- hom-ref: present SNV pos -> count != -1; uncovered pos -> count == -1
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import (
    BRCA1_CHROM,
    BRCA1_START,
    EXPECTED_SAMPLES_TOTAL,
    TRIO_CHILD,
    TRIO_FATHER,
    TRIO_MOTHER,
    VARIANT_KEYS,
    run_cli,
)

pytestmark = pytest.mark.live

# A small, fast sub-window of BRCA1 used for most variant/sample tests.
WIN_START = BRCA1_START
WIN_END = BRCA1_START + 2000
REGION = ["--chrom", BRCA1_CHROM, "--start", str(WIN_START), "--end", str(WIN_END)]


# --- dataset-info ---------------------------------------------------------


def test_dataset_info(live_server):
    d = live_server  # fetched and gated by the probe
    assert d["command"] == "dataset-info"
    assert d["samples_total"] == EXPECTED_SAMPLES_TOTAL
    assert d["females_total"] + d["males_total"] == d["samples_total"]
    assert d["assembly"] == "GRCh38"
    assert isinstance(d["cohorts"], list) and len(d["cohorts"]) >= 1


def test_default_output_path(cli, live_server):
    r = cli(["dataset-info"], use_output=False)
    assert r.code == 0
    assert r.saved_path and Path(r.saved_path).is_file()
    assert r.json is not None and r.json["samples_total"] == EXPECTED_SAMPLES_TOTAL
    assert "1000 Genomes Project -" in r.stdout
    assert "[*] Full JSON saved to" in r.stdout


def test_explicit_output(live_server, tmp_path):
    out = tmp_path / "explicit.json"
    r = run_cli(["count-samples", *REGION], output=out)
    assert r.code == 0
    assert out.is_file()
    assert f"Full JSON saved to {out}" in r.stdout


# --- count / select variants ---------------------------------------------


def test_count_variants(cli, live_server):
    r = cli(["count-variants", *REGION])
    assert r.code == 0
    assert r.json["command"] == "count-variants"
    assert isinstance(r.json["count"], int) and r.json["count"] >= 0
    assert isinstance(r.json["result_incomplete"], bool)
    assert "variants match in" in r.stdout


def test_select_variants_schema(cli, live_server):
    r = cli(["select-variants", *REGION, "--limit", "5"])
    assert r.code == 0
    assert r.json["command"] == "select-variants"
    variants = r.json["variants"]
    assert r.json["count_returned"] == len(variants)
    assert len(variants) <= 5
    for v in variants:
        assert set(v.keys()) == VARIANT_KEYS          # exactly the 19 in-scope keys
        assert "cadd_raw" not in v and "cadd_phred" not in v
        assert WIN_START <= v["start"] <= WIN_END
        assert 0.0 <= v["af"] <= 1.0
        assert v["ref"] and v["alt"]


def test_select_variants_limit_truncation(cli, live_server):
    r = cli(["select-variants", *REGION, "--limit", "5"])
    assert r.code == 0
    # The window has many variants, so a cap of 5 should truncate.
    if r.json["count_returned"] == 5:
        assert r.json["truncated"] is True
        assert "Truncated at --limit" in r.stdout


def test_select_variants_page_size_full_walk(cli, live_server):
    total = cli(["count-variants", *REGION]).json["count"]
    r = cli(["select-variants", *REGION, "--page-size", "500"], timeout=180)
    assert r.code == 0
    assert r.json["truncated"] is False
    # Verified invariant: full walk returns exactly the counted set.
    assert r.json["count_returned"] == total


def test_stdout_preview_cap(cli, live_server):
    r = cli(["select-variants", *REGION, "--limit", "50"])
    assert r.code == 0
    preview = [ln for ln in r.stdout.splitlines() if ln.startswith("  chr")]
    assert len(preview) <= 10
    if r.json["count_returned"] > 10:
        assert len(preview) == 10  # PREVIEW_ROWS cap; full set only in the file


# --- in-samples variants --------------------------------------------------


def test_count_variants_in_samples_subset(cli, live_server):
    cohort = cli(["count-variants", *REGION])
    sub = cli(["count-variants-in-samples", *REGION, "--samples", TRIO_CHILD])
    assert sub.code == 0
    assert sub.json["command"] == "count-variants-in-samples"
    assert sub.json["request"]["samples"] == [TRIO_CHILD]
    assert sub.json["count"] <= cohort.json["count"]  # one individual is a subset


def test_select_variants_in_samples(cli, live_server):
    r = cli(["select-variants-in-samples", *REGION,
             "--samples", f"{TRIO_CHILD},{TRIO_MOTHER}", "--limit", "5"])
    assert r.code == 0
    assert r.json["request"]["samples"] == [TRIO_CHILD, TRIO_MOTHER]
    for v in r.json["variants"]:
        assert set(v.keys()) == VARIANT_KEYS


# --- count / select samples ----------------------------------------------


def test_count_samples(cli, live_server):
    r = cli(["count-samples", *REGION])
    assert r.code == 0
    assert isinstance(r.json["count"], int) and r.json["count"] >= 0
    assert "individuals carry a matching variant" in r.stdout


def test_select_samples(cli, live_server):
    r = cli(["select-samples", *REGION, "--limit", "5"])
    assert r.code == 0
    assert r.json["count"] == len(r.json["samples"])
    assert len(r.json["samples"]) <= 5
    for name in r.json["samples"]:
        assert isinstance(name, str) and name


def test_select_samples_skip_limit(cli, live_server):
    r = cli(["select-samples", *REGION, "--skip", "5", "--limit", "10"])
    assert r.code == 0
    assert r.json["request"]["skip"] == 5
    assert r.json["request"]["limit"] == 10
    assert len(r.json["samples"]) <= 10


# --- zygosity & filter invariants ----------------------------------------


def test_zygosity_monotonic(cli, live_server):
    both = cli(["count-samples", *REGION])
    het = cli(["count-samples", *REGION, "--het-only"])
    hom = cli(["count-samples", *REGION, "--hom-only"])
    assert het.json["count"] <= both.json["count"]
    assert hom.json["count"] <= both.json["count"]
    assert both.json["request"]["zygosity"] == "hom+het"
    assert het.json["request"]["zygosity"] == "het only"
    assert hom.json["request"]["zygosity"] == "hom only"


def test_filter_monotonic(cli, live_server):
    base = cli(["count-variants", *REGION])
    filt = cli(["count-variants", *REGION, "--consequence", "MISSENSE_VARIANT"])
    assert filt.code == 0
    assert filt.json["count"] <= base.json["count"]  # narrowing never increases
    assert filt.json["request"]["filters"]["consequence"] == ["MISSENSE_VARIANT"]


def test_multi_region(cli, live_server):
    mid = (WIN_START + WIN_END) // 2
    a = cli(["count-variants", "--chrom", BRCA1_CHROM, "--start", str(WIN_START), "--end", str(mid)])
    b = cli(["count-variants", "--chrom", BRCA1_CHROM, "--start", str(mid + 1), "--end", str(WIN_END)])
    multi = cli([
        "count-variants",
        "--region", f"{BRCA1_CHROM}:{WIN_START}-{mid}",
        "--region", f"{BRCA1_CHROM}:{mid + 1}-{WIN_END}",
    ])
    assert multi.code == 0
    assert multi.json["request"]["regions"] == [
        f"{BRCA1_CHROM}:{WIN_START}-{mid}", f"{BRCA1_CHROM}:{mid + 1}-{WIN_END}",
    ]
    # Disjoint union count is at least each part.
    assert multi.json["count"] >= a.json["count"]
    assert multi.json["count"] >= b.json["count"]


# --- homozygous-reference sentinels --------------------------------------


def _derive_positions(cli):
    """Return (present_snv_pos, uncovered_pos) derived from real variant data."""
    sv = cli(["select-variants", *REGION, "--page-size", "500"], timeout=180)
    assert sv.code == 0 and sv.json["variants"], "expected variants to derive positions"
    occupied = set()
    snv = []
    for v in sv.json["variants"]:
        for p in range(v["start"], v["end"] + 1):
            occupied.add(p)
        if v["start"] == v["end"]:
            snv.append(v["start"])
    present = min(snv) if snv else min(v["start"] for v in sv.json["variants"])
    uncovered = next((p for p in range(WIN_START, WIN_END + 1) if p not in occupied), None)
    assert uncovered is not None
    return present, uncovered


def test_hom_ref_count_sentinels(cli, live_server):
    present, uncovered = _derive_positions(cli)

    # Present position: a variant exists -> count != -1.
    pr = cli(["count-samples-hom-ref", "--chrom", BRCA1_CHROM, "--position", str(present)])
    assert pr.code == 0
    assert pr.json["variant_present"] is True
    assert pr.json["count"] != -1 and pr.json["count"] >= 0

    # Uncovered position: no variant -> count == -1 (sentinel).
    ab = cli(["count-samples-hom-ref", "--chrom", BRCA1_CHROM, "--position", str(uncovered)])
    assert ab.code == 0
    assert ab.json["count"] == -1
    assert ab.json["variant_present"] is False
    assert "homozygous-reference count is undefined here" in ab.stdout


def test_select_samples_hom_ref(cli, live_server):
    present, _ = _derive_positions(cli)
    r = cli(["select-samples-hom-ref", "--chrom", BRCA1_CHROM, "--position", str(present)])
    assert r.code == 0
    assert r.json["count"] == len(r.json["samples"])
    assert "homozygous reference at" in r.stdout


# --- kinship --------------------------------------------------------------


def test_kinship_related_first_degree(cli, live_server):
    r = cli(["kinship", "--sample1", TRIO_MOTHER, "--sample2", TRIO_CHILD])
    assert r.code == 0
    assert r.json["command"] == "kinship"
    assert r.json["degree"] == "FIRST_DEGREE"      # verified live
    assert isinstance(r.json["phi_bwf"], float)
    assert "KING kinship coefficient phi" in r.stdout


def test_kinship_unrelated_founders(cli, live_server):
    r = cli(["kinship", "--sample1", TRIO_MOTHER, "--sample2", TRIO_FATHER])
    assert r.code == 0
    assert r.json["degree"] == "UNRELATED"          # verified live (founders)


def test_kinship_unknown_sample_errors(cli, live_server):
    """An unknown kinship sample yields a clean error, not a result.

    The server returns an empty pair list for an unknown sample (verified live
    2026-06-26). The wrapper's empty-``pairs`` guard turns this into a clean
    exit 1 with no JSON file and no traceback.

    Historical note: an earlier server bug instead dropped the unknown name and
    returned a degenerate self-comparison (e.g.
    ``NA19240 <-> NA19240: TWINS_MONOZYGOTIC``, exit 0). That bug is fixed; this
    test gates the corrected behavior.
    """
    r = cli(["kinship", "--sample1", "NOT_A_SAMPLE_XYZ", "--sample2", TRIO_CHILD])
    assert r.code == 1
    assert "no relatedness result returned" in r.stderr
    assert "Traceback" not in r.stderr
    assert r.json is None  # nothing written on the error path
