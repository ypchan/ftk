from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from rich_argparse import RichHelpFormatter

from ftk.core.stats import (
    Metric,
    ALL_METRICS,
    compute_assembly_stats,
    format_tsv_rows,
    per_contig_tsv_rows,
)
from ftk.core.gc_peaks import GcPeaksMethod


def add_stat_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "stat",
        help="Compute genome FASTA assembly statistics.",
        description=(
            "Compute genome assembly statistics from FASTA (plain or .gz). Default: output all metrics.\n\n"
            "Input:\n"
            "  - One or more FASTA files: ftk stat a.fa b.fa.gz\n"
            "  - Stdin: cat a.fa | ftk stat -\n"
            "  - If no FASTA args are given, reads from stdin: cat a.fa | ftk stat\n\n"
            "Output:\n"
            "  - TSV to stdout (one row per input)."
        ),
        formatter_class=RichHelpFormatter,
    )

    p.add_argument(
        "fasta",
        nargs="*",
        help="Input genome FASTA file(s). Use '-' for stdin. If omitted, reads from stdin.",
    )

    p.add_argument(
        "--metrics",
        "-m",
        nargs="+",
        choices=[m.value for m in ALL_METRICS],
        default=None,
        help="Choose metrics to output. Default: all.",
    )

    p.add_argument(
        "--threads",
        "-t",
        type=int,
        default=0,
        help="Number of worker processes. 0 means use all CPUs.",
    )

    p.add_argument(
        "--n-run-scaffold-threshold",
        type=int,
        default=2,
        help="If any run of consecutive 'N' is >= this threshold, classify as scaffold.",
    )

    p.add_argument(
        "--per-contig",
        type=Path,
        default=None,
        help="Write per-record stats TSV (source, name, length, N_count, unknown_count, gc, contig_or_scaffold).",
    )

    p.add_argument(
        "--gc-peaks-method",
        choices=[m.value for m in GcPeaksMethod],
        default=GcPeaksMethod.kde.value,
        help="Method for gc_peaks: kde (scipy) or gmm (scikit-learn).",
    )

    p.add_argument(
        "--gc-kde-grid",
        type=int,
        default=1024,
        help="Grid size for KDE peak detection.",
    )

    p.add_argument(
        "--gc-kde-min-prominence",
        type=float,
        default=0.01,
        help="Minimum relative KDE peak height (0-1) to count as a peak.",
    )

    p.add_argument(
        "--gc-gmm-max-components",
        type=int,
        default=6,
        help="Max components for GMM (BIC).",
    )


def run_stat(args: argparse.Namespace) -> None:
    inputs = args.fasta if args.fasta else ["-"]

    metrics: Optional[list[Metric]] = None
    if args.metrics is not None:
        metrics = [Metric(m) for m in args.metrics]

    results = []
    per_contig_rows_all: list[dict] = []

    for inp in inputs:
        source = "stdin" if inp == "-" else inp

        res = compute_assembly_stats(
            source=source,
            fasta_path=None if inp == "-" else Path(inp),
            use_stdin=(inp == "-"),
            n_run_scaffold_threshold=args.n_run_scaffold_threshold,
            metrics=metrics,
            gc_peaks_method=GcPeaksMethod(args.gc_peaks_method),
            gc_kde_grid=args.gc_kde_grid,
            gc_kde_min_prominence=args.gc_kde_min_prominence,
            gc_gmm_max_components=args.gc_gmm_max_components,
        )
        results.append(res)

        if args.per_contig is not None and res.per_contig is not None:
            per_contig_rows_all.extend(per_contig_tsv_rows(res))

    # Pure TSV output only (stdout)
    print(format_tsv_rows(results, metrics=metrics), end="")

    if args.per_contig is not None:
        _write_per_contig(args.per_contig, per_contig_rows_all)


def _write_per_contig(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if rows:
            cols = list(rows[0].keys())
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
        else:
            f.write("source\tname\tlength\tN_count\tunknown_count\tgc\tcontig_or_scaffold\n")
