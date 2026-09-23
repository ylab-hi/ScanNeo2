"""Regression tests for the IEDB prediction subprocesses' stdin (PR #198).

IEDB's `predict_binding.py` and `mhc_II_binding.py` both contain:

    if not sys.stdin.isatty():
        stdin = sys.stdin.readline().strip()

to support piped input. If we invoke them without redirecting stdin they
inherit ours, and under the SLURM executor every job step runs via `srun`,
which supplies an open stdin pipe that is never written and never closed.
That `readline()` then blocks until `PREDICTION_TIMEOUT_SEC` fires and the
batch is discarded -- silently, for every batch, yielding header-only
neoepitope tables while the run still reports progress.

Passing `stdin=subprocess.DEVNULL` gives them an immediate EOF, which is the
input shape they already handle. These tests pin that, because the failure
mode raises nothing and exits zero: there is no other in-band signal.
"""

import ast
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PRIORITIZATION = REPO_ROOT / "workflow/scripts/prioritization"

sys.path.insert(0, str(PRIORITIZATION))
import prediction  # noqa: E402

# One header line (skipped by the parser) plus one binder row. Field order
# follows predict_binding.py's netmhcpan output: index 1 seq_num, 2 start,
# 3 end, 5 peptide, 8 ic50, 9 rank.
STUB_TOOL = '''\
import sys
if not sys.stdin.isatty():
    sys.stdin.readline()
print("allele\\tseq_num\\tstart\\tend\\tlength\\tpeptide\\tmethod\\tscore\\tic50\\trank")
print("HLA-A*02:01\\t1\\t1\\t9\\t9\\tSIINFEKLL\\tnetmhcpan\\t0.5\\t123.4\\t0.35")
'''

EXPECTED = {1: {"SIINFEKLL": ("HLA-A*02:01", 0, 8, 123.4, 0.35)}}


def test_run_prediction_passes_devnull_as_stdin(monkeypatch):
    """The subprocess must be given an explicit DEVNULL stdin.

    Guards the keyword itself: inheriting stdin is what caused the hang.
    """
    captured = {}

    def fake_run(call, **kwargs):
        captured.update(kwargs)
        return subprocess.CompletedProcess(call, 0, stdout="header\n", stderr="")

    monkeypatch.setattr(prediction.subprocess, "run", fake_run)
    prediction.BindingAffinities._run_prediction(
        ["true"], "batch.fa", "mt", "mhc-I"
    )

    assert "stdin" in captured, "_run_prediction must not let the child inherit stdin"
    assert captured["stdin"] is subprocess.DEVNULL


def test_run_prediction_survives_an_open_stdin_pipe(tmp_path, monkeypatch):
    """End-to-end against a tool that mimics IEDB's stdin idiom.

    fd 0 is replaced with the read end of a pipe that is never written, which
    is what srun hands a job step. Without `stdin=DEVNULL` the child inherits
    it, blocks in readline(), and _run_prediction returns None.
    """
    stub = tmp_path / "stub_iedb_tool.py"
    stub.write_text(STUB_TOOL)

    # Keep the failure fast: a regression would otherwise wait the full hour.
    monkeypatch.setattr(prediction, "PREDICTION_TIMEOUT_SEC", 20)

    read_fd, write_fd = os.pipe()
    saved_stdin = os.dup(0)
    try:
        os.dup2(read_fd, 0)  # an open pipe nobody will ever write to
        result = prediction.BindingAffinities._run_prediction(
            [sys.executable, str(stub)], str(stub), "mt", "mhc-I"
        )
    finally:
        os.dup2(saved_stdin, 0)
        for fd in (saved_stdin, read_fd, write_fd):
            os.close(fd)

    assert result is not None, (
        "prediction batch was dropped -- the child blocked on the inherited "
        "stdin pipe instead of receiving EOF"
    )
    assert result == EXPECTED


def _subprocess_run_calls(path):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "run"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "subprocess"
        ):
            yield node


def _passes_devnull_stdin(call_node):
    for kw in call_node.keywords:
        if kw.arg == "stdin":
            return (
                isinstance(kw.value, ast.Attribute)
                and kw.value.attr == "DEVNULL"
                and isinstance(kw.value.value, ast.Name)
                and kw.value.value.id == "subprocess"
            )
    return False


def test_immunogenicity_call_closes_stdin():
    """filtering.py shells out to a third IEDB tool and needs the same guard.

    Checked via AST rather than by importing: filtering.py needs pandas and
    blosum, which the CI test env does not install. This call has no timeout,
    so an inherited stdin would hang unbounded rather than for an hour.
    """
    filtering = PRIORITIZATION / "filtering.py"
    matches = [
        node
        for node in _subprocess_run_calls(filtering)
        if any(
            isinstance(arg, ast.Constant)
            and isinstance(arg.value, str)
            and "predict_immunogenicity.py" in arg.value
            for child in node.args
            for arg in ast.walk(child)
        )
    ]

    assert matches, "no subprocess.run invoking predict_immunogenicity.py found"
    for node in matches:
        assert _passes_devnull_stdin(node), (
            f"{filtering.name}:{node.lineno} invokes an IEDB tool without "
            "stdin=subprocess.DEVNULL"
        )
