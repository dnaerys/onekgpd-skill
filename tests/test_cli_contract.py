"""Offline CLI-contract tests.

Every test here triggers an error in argparse or in the in-process input
builders (``_build_regions`` / ``_build_annotation_filter``), which run *before*
any network call. They need no live server and always run. Each asserts the exit
code, a specific stderr message, and that no traceback leaks.

Exit-code convention:
- 1  -> handled error (clean ``Error: ...`` to stderr)
- 2  -> argparse usage error (mutually-exclusive / missing-required)
"""

from __future__ import annotations

# A valid region prefix so that region construction succeeds and a *filter*
# error is the one that surfaces.
VALID_REGION = ["--chrom", "chr17", "--start", "1", "--end", "2"]


def test_help_lists_all_commands(cli):
    r = cli(["--help"], use_output=False)
    assert r.code == 0
    for cmd in (
        "dataset-info", "count-variants", "select-variants",
        "count-variants-in-samples", "select-variants-in-samples",
        "count-samples", "select-samples", "count-samples-hom-ref",
        "select-samples-hom-ref", "kinship",
    ):
        assert cmd in r.stdout
    assert "Traceback" not in r.stderr


def test_subcommand_help(cli):
    r = cli(["select-variants", "--help"], use_output=False)
    assert r.code == 0
    assert "--limit" in r.stdout and "--page-size" in r.stdout
    assert "Traceback" not in r.stderr


def test_err_missing_coords(cli):
    r = cli(["count-variants", "--chrom", "chr17"], use_output=False)
    assert r.code == 1
    assert "Error: --chrom requires --start and --end" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_region_and_chrom(cli):
    r = cli(["count-variants", *VALID_REGION, "--region", "chr1:1-2"], use_output=False)
    assert r.code == 1
    assert "use either --chrom/--start/--end or --region" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_bad_region_format(cli):
    r = cli(["count-variants", "--region", "chr1_1_2"], use_output=False)
    assert r.code == 1
    assert "invalid --region" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_ref_with_multiregion(cli):
    r = cli(["count-variants", "--region", "chr1:1-2", "--ref", "A"], use_output=False)
    assert r.code == 1
    assert "--ref/--alt apply only to a single" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_unknown_enum_lists_valid(cli):
    r = cli(["count-variants", *VALID_REGION, "--consequence", "NOT_A_TERM"], use_output=False)
    assert r.code == 1
    assert "NOT_A_TERM" in r.stderr
    assert "MISSENSE_VARIANT" in r.stderr  # message enumerates valid values
    assert "Traceback" not in r.stderr


def test_err_am_class_vs_score_conflict(cli):
    r = cli(
        ["count-variants", *VALID_REGION,
         "--alpha-missense-class", "AM_AMBIGUOUS",
         "--alpha-missense-score-lt", "0.5"],
        use_output=False,
    )
    assert r.code == 1
    assert "alpha-missense-class cannot be combined" in r.stderr
    assert "Traceback" not in r.stderr


def test_err_zygosity_mutually_exclusive(cli):
    r = cli(["count-variants", *VALID_REGION, "--het-only", "--hom-only"], use_output=False)
    assert r.code == 2  # argparse
    assert "not allowed with argument" in r.stderr


def test_err_biallelic_mutually_exclusive(cli):
    r = cli(["count-variants", *VALID_REGION, "--biallelic-only", "--multiallelic-only"], use_output=False)
    assert r.code == 2


def test_err_exclude_sex_mutually_exclusive(cli):
    r = cli(["count-variants", *VALID_REGION, "--exclude-males", "--exclude-females"], use_output=False)
    assert r.code == 2


def test_err_missing_required_samples(cli):
    r = cli(["count-variants-in-samples", *VALID_REGION], use_output=False)
    assert r.code == 2  # --samples is required


def test_err_missing_required_sample1(cli):
    r = cli(["kinship", "--sample2", "NA19240"], use_output=False)
    assert r.code == 2  # --sample1 is required


def test_err_no_subcommand(cli):
    r = cli([], use_output=False)
    assert r.code == 2  # a subcommand is required
