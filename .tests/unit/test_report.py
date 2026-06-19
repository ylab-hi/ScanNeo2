"""Tests for workflow/scripts/report.py (issue #94).

Each test materializes a synthetic ``results/<sample>/prioritization/`` tree,
a synthetic ``.snakemake/log/*.snakemake.log``, and an optional
``config/config.yaml`` under ``tmp_path``, then runs ``report.py`` from that
working directory and asserts on stdout.

The script tolerates a missing PyYAML; class-aware completion checks expect
PyYAML to be importable in the test env (the same env runs snakemake itself).
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow/scripts/report.py"


def write_log(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / ".snakemake" / "log" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body)
    return p


def write_results(
    tmp_path: Path,
    sample: str,
    *,
    marker: bool = False,
    mhc_i: bool = False,
    mhc_ii: bool = False,
    mhc_i_text: str = "record\n",
    mhc_ii_text: str = "record\n",
) -> Path:
    pri = tmp_path / "results" / sample / "prioritization"
    pri.mkdir(parents=True)
    if marker:
        (pri / ".snakemake_timestamp").touch()
    if mhc_i:
        (pri / "mhc-I_neoepitopes_all.txt").write_text(mhc_i_text)
    if mhc_ii:
        (pri / "mhc-II_neoepitopes_all.txt").write_text(mhc_ii_text)
    return pri


def write_config(tmp_path: Path, prio_class: str) -> Path:
    cfg = tmp_path / "config" / "config.yaml"
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text(
        "prioritization:\n"
        f"  class: {prio_class}\n"
        "  lengths:\n"
        "    MHC-I: '8,9'\n"
        "    MHC-II: '15'\n"
    )
    return cfg


def run_report(tmp_path: Path, *extra_args: str) -> "subprocess.CompletedProcess[str]":
    return subprocess.run(
        [sys.executable, str(SCRIPT), *extra_args],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


# ---------------------------------------------------------------------------
# Log-block templates
# ---------------------------------------------------------------------------

def rule_block(rule: str, jobid: int, sample: str, log_path: str) -> str:
    return (
        "[Mon Jun 16 12:00:00 2026]\n"
        f"rule {rule}:\n"
        f"    log: {log_path}\n"
        f"    jobid: {jobid}\n"
        f"    wildcards: sample={sample}\n"
    )


def finished_block(jobid: int) -> str:
    return f"[Mon Jun 16 12:00:30 2026]\nFinished job {jobid}.\n"


def error_block(rule: str, jobid: int, log_path: str) -> str:
    return (
        "[Mon Jun 16 12:01:00 2026]\n"
        f"Error in rule {rule}:\n"
        f"    jobid: {jobid}\n"
        f"    log: {log_path} (check log file(s) for error details)\n"
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_all_complete(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)
    write_results(tmp_path, "S2", marker=True, mhc_i=True)
    body = (
        rule_block("prioritization", 1, "S1", "logs/S1/prioritization/prio.log")
        + finished_block(1)
        + rule_block("prioritization", 2, "S2", "logs/S2/prioritization/prio.log")
        + finished_block(2)
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "2 complete" in r.stdout
    assert "S1" in r.stdout and "S2" in r.stdout


def test_one_rule_errored(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1")  # results dir exists, no completion artifacts
    rule_log = tmp_path / "logs" / "S1" / "calling" / "var_call.log"
    rule_log.parent.mkdir(parents=True, exist_ok=True)
    rule_log.write_text("\n".join(f"diag line {i}" for i in range(30)) + "\n")

    body = (
        rule_block("var_call", 7, "S1", str(rule_log))
        + error_block("var_call", 7, str(rule_log))
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "error" in r.stdout
    assert "var_call" in r.stdout
    assert "diag line 29" in r.stdout
    assert "diag line 14" not in r.stdout  # outside the 15-line tail


def test_incomplete_run(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1")
    body = (
        rule_block("step_a", 1, "S1", "logs/S1/a.log")
        + finished_block(1)
        + rule_block("step_b", 2, "S1", "logs/S1/b.log")
        + rule_block("step_c", 3, "S1", "logs/S1/c.log")
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "incomplete" in r.stdout
    assert "1/3 jobs finished" in r.stdout


def test_not_started_via_samples_filter(tmp_path):
    write_config(tmp_path, "I")
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", "")

    r = run_report(tmp_path, "--samples", "ghost")

    assert r.returncode == 0, r.stderr + r.stdout
    assert "ghost" in r.stdout
    assert "not_started" in r.stdout


def test_multi_sample_mixed_states(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "ok", marker=True, mhc_i=True)
    write_results(tmp_path, "broken")
    write_results(tmp_path, "running")
    broken_log = tmp_path / "logs" / "broken" / "step.log"
    broken_log.parent.mkdir(parents=True, exist_ok=True)
    broken_log.write_text("boom\n")
    body = (
        rule_block("done_rule", 1, "ok", "logs/ok/done.log")
        + finished_block(1)
        + rule_block("step", 2, "broken", str(broken_log))
        + error_block("step", 2, str(broken_log))
        + rule_block("ongoing", 3, "running", "logs/running/x.log")
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "1 complete, 1 errored, 1 incomplete" in r.stdout


def test_malformed_log_lines_skipped(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)
    body = (
        "garbage line\n"
        "@@@ random noise\n"
        + rule_block("done", 1, "S1", "logs/S1/done.log")
        + "    something_weird: value\n"
        + finished_block(1)
        + "more garbage\n"
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "complete" in r.stdout
    assert "S1" in r.stdout


def test_missing_per_rule_log_for_error(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1")
    body = (
        rule_block("var_call", 9, "S1", "logs/S1/nope.log")
        + error_block("var_call", 9, "logs/S1/nope.log")
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "log file not found" in r.stdout


def test_shared_rule_without_wildcards_ignored(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)
    body = (
        "[Mon Jun 16 12:00:00 2026]\n"
        "localrule bwa_index:\n"
        "    log: logs/ref/bwa_index.log\n"
        "    jobid: 42\n"
        "[Mon Jun 16 12:00:30 2026]\n"
        "Finished job 42.\n"
        + rule_block("done", 1, "S1", "logs/S1/done.log")
        + finished_block(1)
    )
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "complete" in r.stdout
    assert "S1" in r.stdout
    # The shared rule must not appear as a row in its own right.
    assert "\n| bwa_index " not in r.stdout


def test_markdown_output(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)
    body = rule_block("done", 1, "S1", "logs/S1/done.log") + finished_block(1)
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path, "--markdown")

    assert r.returncode == 0, r.stderr + r.stdout
    assert "| --- |" in r.stdout
    assert "| Sample |" in r.stdout


def test_picks_latest_master_log(tmp_path):
    write_config(tmp_path, "I")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)
    # The older log lacks any rule for S1 (so S1 would be classified by
    # results state alone) — same effect, but the report header must name the
    # newer log either way.
    write_log(tmp_path, "2026-01-01T000000.0.snakemake.log", "")
    write_log(
        tmp_path,
        "2026-06-16T120000.0.snakemake.log",
        rule_block("done", 1, "S1", "logs/S1/done.log") + finished_block(1),
    )

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "2026-06-16T120000.0.snakemake.log" in r.stdout
    assert "2026-01-01" not in r.stdout


def test_complete_class_II_only(tmp_path):
    write_config(tmp_path, "II")
    write_results(tmp_path, "S1", marker=True, mhc_ii=True)
    body = rule_block("done", 1, "S1", "logs/S1/done.log") + finished_block(1)
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "complete" in r.stdout
    assert "Prioritization class: II" in r.stdout


def test_class_BOTH_requires_both_files(tmp_path):
    """BOTH config + only class-I combined output present → incomplete, not complete."""
    write_config(tmp_path, "BOTH")
    write_results(tmp_path, "S1", marker=True, mhc_i=True)  # mhc-II missing
    body = rule_block("done", 1, "S1", "logs/S1/done.log") + finished_block(1)
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    # Sample must NOT classify as complete; it has one finished job, so
    # classification falls through to incomplete (or remains incomplete via
    # no-jobs-for-sample path if log carries no S1 rows — both acceptable).
    # We assert the strong claim: the summary shows 0 complete.
    assert "0 complete" in r.stdout


def test_fallback_when_config_missing(tmp_path):
    """No config/config.yaml — script warns, falls back to any-one-combined."""
    write_results(tmp_path, "S1", marker=True, mhc_i=True)  # only mhc-I present
    body = rule_block("done", 1, "S1", "logs/S1/done.log") + finished_block(1)
    write_log(tmp_path, "2026-06-16T120000.0.snakemake.log", body)

    r = run_report(tmp_path)

    assert r.returncode == 0, r.stderr + r.stdout
    assert "WARNING" in r.stderr
    assert "1 complete" in r.stdout
