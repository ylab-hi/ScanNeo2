"""Post-run status report for ScanNeo2.

Scans the most recent Snakemake master log (``.snakemake/log/*.snakemake.log``)
together with the on-disk ``results/<sample>/prioritization/`` artifacts and
emits a per-sample status table (complete / error / incomplete / not started).
For errored samples it appends the tail of the per-rule log file so the cause
is one read away. See issue #94.

The completion check follows ``config['prioritization']['class']``: a run
configured for class I requires ``mhc-I_neoepitopes_all.txt`` non-empty, class
II requires ``mhc-II_neoepitopes_all.txt`` non-empty, and BOTH requires both.
If the config can't be loaded (no PyYAML, missing path, malformed) the script
falls back to "marker plus any-one combined output non-empty" and warns on
stderr.

Fragility: this script parses ``.snakemake/log/*.snakemake.log``, which is
**human-facing console output, not a stable Snakemake API**. Upstream changes
to Snakemake's log format (job header style, ``Finished job N.`` wording, the
``Error in rule X:`` block) will break this parser. If the report misclassifies
after a Snakemake upgrade, re-validate the regexes at the top of this file
against a fresh log file. Tested against the log format emitted by
ScanNeo2 v0.4.x running Snakemake >= 8.0.

Usage:
    python workflow/scripts/report.py
    python workflow/scripts/report.py --markdown --output status.md
"""

import argparse
import re
import sys
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

try:
    import yaml  # type: ignore
    _HAVE_YAML = True
except ImportError:
    _HAVE_YAML = False


RE_TS = re.compile(
    r"^\[(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) "
    r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) "
    r"\d{1,2} \d{2}:\d{2}:\d{2} \d{4}\]$"
)
RE_RULE_OPEN = re.compile(r"^(rule|localrule|checkpoint) ([A-Za-z_]\w*):$")
RE_COMPACT = re.compile(r"^Job (\d+): (.+)$")
RE_FINISHED = re.compile(r"^Finished job (\d+)\.$")
RE_ERROR_OPEN = re.compile(r"^Error in rule ([A-Za-z_]\w*):$")
RE_JOBID = re.compile(r"^    jobid: (\d+)$")
RE_LOG_FIELD = re.compile(r"^    log: (\S+)(?: \(.*\))?$")
RE_WILDCARDS = re.compile(r"^    wildcards: (.+)$")
RE_SAMPLE_INLINE = re.compile(r"\bsample:([A-Za-z0-9_.-]+)")

CLASS_TO_COMBINED = {
    "I": ["mhc-I_neoepitopes_all.txt"],
    "II": ["mhc-II_neoepitopes_all.txt"],
    "BOTH": ["mhc-I_neoepitopes_all.txt", "mhc-II_neoepitopes_all.txt"],
}


@dataclass
class JobBlock:
    jobid: int
    rule: str = "?"
    sample: Optional[str] = None
    log_path: Optional[str] = None
    finished: bool = False
    errored: bool = False
    error_log_path: Optional[str] = None
    description: Optional[str] = None


@dataclass
class SampleStatus:
    sample: str
    status: str  # one of: complete, error, incomplete, not_started
    failed_rule: Optional[str] = None
    failed_log: Optional[str] = None
    finished_jobs: int = 0
    total_jobs: int = 0
    excerpt: List[str] = field(default_factory=list)


@dataclass
class ResultsInfo:
    marker_present: bool = False
    combined_files: Dict[str, bool] = field(default_factory=dict)


class _Pending:
    __slots__ = ("kind", "rule", "jobid", "log_path", "sample")

    def __init__(self, kind: str, rule: Optional[str] = None) -> None:
        self.kind = kind
        self.rule = rule
        self.jobid: Optional[int] = None
        self.log_path: Optional[str] = None
        self.sample: Optional[str] = None


def _parse_sample_from_wildcards(wc: str) -> Optional[str]:
    for kv in wc.split(", "):
        if "=" in kv:
            k, v = kv.split("=", 1)
            if k == "sample":
                return v
    return None


def parse_master_log(path: Path) -> List[JobBlock]:
    """Parse a single Snakemake master log into a list of ``JobBlock``s."""
    jobs: Dict[int, JobBlock] = {}
    pending: Optional[_Pending] = None

    def flush(p: Optional[_Pending]) -> None:
        if p is None or p.jobid is None:
            return
        jb = jobs.get(p.jobid)
        if jb is None:
            jb = JobBlock(jobid=p.jobid)
            jobs[p.jobid] = jb
        if p.rule and (jb.rule == "?" or p.kind == "error"):
            jb.rule = p.rule
        if p.sample and jb.sample is None:
            jb.sample = p.sample
        if p.kind == "rule" and p.log_path and jb.log_path is None:
            jb.log_path = p.log_path
        if p.kind == "error":
            jb.errored = True
            if p.log_path and jb.error_log_path is None:
                jb.error_log_path = p.log_path

    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if RE_TS.match(line):
                flush(pending)
                pending = None
                continue

            m = RE_RULE_OPEN.match(line)
            if m:
                flush(pending)
                pending = _Pending(kind="rule", rule=m.group(2))
                continue

            m = RE_ERROR_OPEN.match(line)
            if m:
                flush(pending)
                pending = _Pending(kind="error", rule=m.group(1))
                continue

            m = RE_COMPACT.match(line)
            if m:
                flush(pending)
                pending = None
                jobid = int(m.group(1))
                desc = m.group(2)
                ms = RE_SAMPLE_INLINE.search(desc)
                sample = ms.group(1) if ms else None
                existing = jobs.get(jobid)
                if existing is None:
                    jobs[jobid] = JobBlock(
                        jobid=jobid, rule="?", sample=sample, description=desc
                    )
                else:
                    if existing.sample is None and sample:
                        existing.sample = sample
                    if existing.description is None:
                        existing.description = desc
                continue

            m = RE_FINISHED.match(line)
            if m:
                flush(pending)
                pending = None
                jobid = int(m.group(1))
                if jobid in jobs:
                    jobs[jobid].finished = True
                continue

            if pending is None:
                continue

            m = RE_JOBID.match(line)
            if m:
                pending.jobid = int(m.group(1))
                continue

            m = RE_LOG_FIELD.match(line)
            if m:
                pending.log_path = m.group(1)
                continue

            m = RE_WILDCARDS.match(line)
            if m:
                pending.sample = _parse_sample_from_wildcards(m.group(1))
                continue

        flush(pending)

    return list(jobs.values())


def find_latest_master_log(log_dir: Path) -> Optional[Path]:
    """Return the lex-greatest ``*.snakemake.log`` under ``log_dir`` or None.

    Snakemake writes one log per invocation, named with an ISO-8601-ish
    timestamp; lexicographic max coincides with chronological max.
    """
    if not log_dir.is_dir():
        return None
    candidates = sorted(log_dir.glob("*.snakemake.log"))
    return candidates[-1] if candidates else None


def load_prioritization_class(
    config_path: Optional[Path],
) -> Optional[str]:
    """Return ``config['prioritization']['class']`` or ``None`` on any failure.

    Failures (missing path, no PyYAML, malformed YAML, missing keys) are
    silently absorbed and surfaced as ``None`` so the caller can warn once
    and fall back to the looser completion check.
    """
    if config_path is None or not config_path.exists() or not _HAVE_YAML:
        return None
    try:
        with open(config_path, encoding="utf-8") as fh:
            cfg = yaml.safe_load(fh)
    except (OSError, yaml.YAMLError):
        return None
    if not isinstance(cfg, dict):
        return None
    cls = cfg.get("prioritization", {}).get("class") if isinstance(
        cfg.get("prioritization"), dict
    ) else None
    if cls in ("I", "II", "BOTH"):
        return cls
    return None


def scan_results_dir(
    results_dir: Path, required_combined: List[str]
) -> Dict[str, ResultsInfo]:
    """Inspect ``results/<sample>/prioritization/`` for each subdir.

    ``required_combined`` is the union of combined-output filenames the report
    may need to check across all completion regimes; the per-class subset is
    applied at classification time.
    """
    info: Dict[str, ResultsInfo] = {}
    if not results_dir.is_dir():
        return info
    for sub in sorted(results_dir.iterdir()):
        if not sub.is_dir():
            continue
        pri = sub / "prioritization"
        marker = (pri / ".snakemake_timestamp").exists()
        combined = {}
        for fname in required_combined:
            p = pri / fname
            combined[fname] = p.exists() and p.stat().st_size > 0
        info[sub.name] = ResultsInfo(marker_present=marker, combined_files=combined)
    return info


def _is_complete(info: Optional[ResultsInfo], required: List[str]) -> bool:
    if info is None or not info.marker_present:
        return False
    if not required:
        return any(info.combined_files.values())
    return all(info.combined_files.get(f, False) for f in required)


def classify(
    sample: str,
    results_info: Optional[ResultsInfo],
    jobs_for_sample: List[JobBlock],
    required_combined: List[str],
) -> SampleStatus:
    """Apply the locked decision tree (#94)."""
    if _is_complete(results_info, required_combined):
        return SampleStatus(sample=sample, status="complete")

    errored = [j for j in jobs_for_sample if j.errored]
    if errored:
        first = errored[0]
        return SampleStatus(
            sample=sample,
            status="error",
            failed_rule=first.rule,
            failed_log=first.error_log_path or first.log_path,
        )

    if not jobs_for_sample and results_info is None:
        return SampleStatus(sample=sample, status="not_started")

    total = len(jobs_for_sample)
    finished = sum(1 for j in jobs_for_sample if j.finished)
    return SampleStatus(
        sample=sample,
        status="incomplete",
        finished_jobs=finished,
        total_jobs=total,
    )


def tail_log(path: Optional[Path], n: int) -> List[str]:
    if path is None:
        return ["<no log path recorded>"]
    if not path.exists():
        return [f"<log file not found: {path}>"]
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            tail = deque(fh, maxlen=n)
    except OSError as exc:
        return [f"<could not read log: {exc}>"]
    lines = [t.rstrip("\n") for t in tail]
    return lines if lines else ["<log file empty>"]


def _ascii_table(rows: List[List[str]]) -> str:
    if not rows:
        return ""
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    out = []
    for r_idx, row in enumerate(rows):
        cells = " | ".join(c.ljust(widths[i]) for i, c in enumerate(row))
        out.append(f"| {cells} |")
        if r_idx == 0:
            sep = "-+-".join("-" * widths[i] for i in range(len(row)))
            out.append(f"+-{sep}-+")
    return "\n".join(out)


def _md_table(rows: List[List[str]]) -> str:
    if not rows:
        return ""
    out = ["| " + " | ".join(rows[0]) + " |"]
    out.append("| " + " | ".join("---" for _ in rows[0]) + " |")
    for row in rows[1:]:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def _summary_line(statuses: List[SampleStatus]) -> str:
    n_complete = sum(1 for s in statuses if s.status == "complete")
    n_error = sum(1 for s in statuses if s.status == "error")
    n_incomplete = sum(1 for s in statuses if s.status == "incomplete")
    n_total = len(statuses)
    return (
        f"{n_complete} complete, {n_error} errored, "
        f"{n_incomplete} incomplete of {n_total} samples"
    )


def _notes(s: SampleStatus) -> str:
    if s.status == "incomplete" and s.total_jobs:
        return f"{s.finished_jobs}/{s.total_jobs} jobs finished"
    return ""


def render(
    statuses: List[SampleStatus],
    master_log: Optional[Path],
    prioritization_class: Optional[str],
    show_excerpts: bool,
    markdown: bool,
) -> str:
    header_rows = [["Sample", "Status", "Failed rule", "Log", "Notes"]]
    for s in statuses:
        header_rows.append(
            [
                s.sample,
                s.status,
                s.failed_rule or "",
                s.failed_log or "",
                _notes(s),
            ]
        )

    parts: List[str] = []
    title = "ScanNeo2 status report"
    if markdown:
        parts.append(f"# {title}")
    else:
        parts.append(title)
        parts.append("=" * len(title))

    parts.append(
        f"Master log: {master_log}" if master_log else "Master log: <none found>"
    )
    parts.append(
        f"Prioritization class: {prioritization_class}"
        if prioritization_class
        else "Prioritization class: <unknown - fell back to any-one combined output>"
    )
    parts.append(f"Summary: {_summary_line(statuses)}")
    parts.append("")
    parts.append(_md_table(header_rows) if markdown else _ascii_table(header_rows))

    if show_excerpts:
        errored = [s for s in statuses if s.status == "error" and s.excerpt]
        if errored:
            parts.append("")
            tail_title = "Errored samples - log tails"
            if markdown:
                parts.append(f"## {tail_title}")
            else:
                parts.append(tail_title)
                parts.append("-" * len(tail_title))
            for s in errored:
                parts.append("")
                ref = f"{s.sample} - {s.failed_rule} ({s.failed_log})"
                if markdown:
                    parts.append(f"### {ref}")
                    parts.append("```")
                    parts.extend(s.excerpt)
                    parts.append("```")
                else:
                    parts.append(f"{ref}:")
                    for line in s.excerpt:
                        parts.append(f"    {line}")

    return "\n".join(parts) + "\n"


def _expected_samples(
    results_info: Dict[str, ResultsInfo],
    jobs: List[JobBlock],
    explicit: Optional[List[str]],
) -> List[str]:
    discovered = set(results_info.keys())
    for j in jobs:
        if j.sample:
            discovered.add(j.sample)
    if explicit:
        return sorted(set(explicit) | discovered)
    return sorted(discovered)


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Per-sample status report for a ScanNeo2 run (issue #94)."
    )
    p.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="root of per-sample results trees (default: results/)",
    )
    p.add_argument(
        "--log-dir",
        type=Path,
        default=Path(".snakemake/log"),
        help="directory holding snakemake master logs (default: .snakemake/log/)",
    )
    p.add_argument(
        "--master-log",
        type=Path,
        default=None,
        help="explicit master log to parse (default: latest under --log-dir)",
    )
    p.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="ScanNeo2 config to read prioritization.class from "
        "(default: config/config.yaml)",
    )
    p.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="restrict / extend the report to these sample names",
    )
    p.add_argument(
        "--markdown",
        action="store_true",
        help="emit GitHub-flavoured markdown instead of ASCII",
    )
    p.add_argument(
        "--excerpt-lines",
        type=int,
        default=15,
        help="lines to print from each errored per-rule log (default: 15)",
    )
    p.add_argument(
        "--no-excerpts",
        action="store_true",
        help="omit per-rule log tails for errored samples",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="write report here instead of stdout",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)

    master_log = args.master_log
    if master_log is None:
        master_log = find_latest_master_log(args.log_dir)
    elif not master_log.exists():
        print(f"ERROR: --master-log not found: {master_log}", file=sys.stderr)
        return 1

    prio_class = load_prioritization_class(args.config)
    if prio_class is None:
        print(
            f"WARNING: could not derive prioritization class from "
            f"{args.config}; falling back to any-one combined output check.",
            file=sys.stderr,
        )
        required_combined: List[str] = []
        scan_targets = CLASS_TO_COMBINED["BOTH"]
    else:
        required_combined = CLASS_TO_COMBINED[prio_class]
        scan_targets = required_combined

    jobs: List[JobBlock] = []
    if master_log is not None:
        jobs = parse_master_log(master_log)

    results_info = scan_results_dir(args.results_dir, scan_targets)

    samples = _expected_samples(results_info, jobs, args.samples)
    if args.samples:
        samples = [s for s in samples if s in set(args.samples)]

    statuses: List[SampleStatus] = []
    for sample in samples:
        sample_jobs = [j for j in jobs if j.sample == sample]
        st = classify(
            sample, results_info.get(sample), sample_jobs, required_combined
        )
        if st.status == "error" and not args.no_excerpts:
            st.excerpt = tail_log(
                Path(st.failed_log) if st.failed_log else None,
                args.excerpt_lines,
            )
        statuses.append(st)

    out = render(
        statuses,
        master_log,
        prio_class,
        show_excerpts=not args.no_excerpts,
        markdown=args.markdown,
    )

    if args.output:
        args.output.write_text(out)
    else:
        sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
