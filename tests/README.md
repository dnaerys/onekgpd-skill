# OneKGPd test suite

Black-box end-to-end tests for the `onekgpd` skill (which lives at
`skills/onekgpd/`). Tests invoke `uv run skills/onekgpd/scripts/onekgpd_api.py …`
as a subprocess and assert on exit code, stdout, stderr, and the JSON output
file. They depend only on the standard library + pytest and never import
`dnaerys` (the opt-in unit tier is the sole exception and is not collected by
default).

## Tiers

- **Offline contract** (`test_cli_contract.py`) — argparse + input-validation
  error paths. No live endpoint needed; always runs.
- **Live E2E** (`test_e2e_live.py`) — the ten commands and edge cases against the
  public endpoint. The session `live_server` fixture probes once: it **skips**
  the live tier on a connection failure and **fails** on wrong data
  (`samples_total != 3202`). An all-skipped live run prints a loud banner.
- **Opt-in unit** (`test_helpers_unit.py`) — in-process tests of `_call_with_retry`
  and pure helpers. Imports `dnaerys`; collected only with `--run-unit`.

## Run

```bash
# Default suite (offline always; live auto-skips if the endpoint is down):
uv run --with pytest pytest tests -v

# Offline contract tests only:
uv run --with pytest pytest tests/test_cli_contract.py -v

# Opt-in unit tier (needs dnaerys importable):
uv run --with pytest --with dnaerys pytest --run-unit tests/test_helpers_unit.py -v
```

Both PyPI (for `uv` to provision `dnaerys`) and the gRPC endpoint must be
reachable for the live tier.
