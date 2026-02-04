from enum import Enum
from typing import List, Optional

import typer

from ftk.core.stats import run_stat

app = typer.Typer(help="Compute genome FASTA assembly statistics.",
                  context_settings={"help_option_names": ["-h", "--help"]})


class Metric(str, Enum):
    genome_size = "genome_size"
    contig_or_scaffold = "contig_or_scaffold"
    unknown_count = "unknown_count"
    longest = "longest"
    shortest = "shortest"
    gc = "gc"
    n_stats = "n_stats"
    gc_avg = "gc_avg"
    gc_max = "gc_max"
    gc_min = "gc_min"
    gc_peaks = "gc_peaks"


class GCPeaksMethod(str, Enum):
    kde = "kde"
    gmm = "gmm"


ALL_METRICS = [
    Metric.genome_size,
    Metric.contig_or_scaffold,
    Metric.unknown_count,
    Metric.longest,
    Metric.shortest,
    Metric.gc,
    Metric.n_stats,
    Metric.gc_avg,
    Metric.gc_max,
    Metric.gc_min,
    Metric.gc_peaks,
]


STAT_HELP = """\
Compute genome assembly statistics from FASTA (plain or .gz). Default: output all metrics.

Definitions:
- genome_size:
  - genome_size_with_N: sum of record lengths (includes N and ambiguous/unknown characters).
  - genome_size_ungapped: genome_size_with_N minus total count of 'N'.
- contig_or_scaffold:
  - If any run of consecutive 'N' is >= --n-run-scaffold-threshold => scaffold; else contig.
- unknown_count:
  - Count of non-IUPAC characters (IUPAC set: A C G T N R Y K M S W B D H V).
- gc:
  - Global GC = (G+C)/(A+C+G+T), aggregated across all records.
- n_stats:
  - N50/L50/N90/L90 based on lengths and genome_size_with_N.
- gc_avg/max/min:
  - Unweighted statistics across per-record GC values.
- gc_peaks:
  - Number of modes in per-record GC distribution via KDE (default) or GMM.
"""


@app.command(help=STAT_HELP)
def main(
    fasta: str = typer.Argument(..., help="Input genome FASTA file (optionally .gz)."),
    metrics: List[Metric] = typer.Option(
        None, "-m", "--metrics",
        help="Choose metrics to output. Default: all."
    ),
    threads: int = typer.Option(
        0, "-t", "--threads",
        help="Number of worker processes. 0 means use all CPUs."
    ),
    n_run_scaffold_threshold: int = typer.Option(
        2, "--n-run-scaffold-threshold",
        help="If any run of consecutive 'N' is >= this threshold, classify as scaffold.",
        min=1,
    ),
    per_contig: Optional[str] = typer.Option(
        None, "--per-contig",
        help="Write per-record stats TSV: name, length, N_count, unknown_count, gc, contig_or_scaffold."
    ),
    gc_peaks_method: GCPeaksMethod = typer.Option(
        GCPeaksMethod.kde, "--gc-peaks-method",
        help="Method for gc_peaks: kde (scipy) or gmm (scikit-learn)."
    ),
    gc_kde_grid: int = typer.Option(
        1024, "--gc-kde-grid", help="Grid size for KDE peak detection.", min=128
    ),
    gc_kde_min_prominence: float = typer.Option(
        0.01, "--gc-kde-min-prominence",
        help="Minimum relative KDE peak height (0-1) to count as a peak.", min=0.0
    ),
    gc_gmm_max_components: int = typer.Option(
        6, "--gc-gmm-max-components",
        help="Max components for GMM (BIC).", min=1
    ),
):
    sel = metrics if metrics else ALL_METRICS
    result = run_stat(
        fasta=fasta,
        metrics=set(sel),
        threads=threads,
        n_run_scaffold_threshold=n_run_scaffold_threshold,
        per_contig_path=per_contig,
        gc_peaks_method=gc_peaks_method.value,
        gc_kde_grid=gc_kde_grid,
        gc_kde_min_prominence=gc_kde_min_prominence,
        gc_gmm_max_components=gc_gmm_max_components,
    )

    typer.echo(f"# ftk stat: {result['input']}")
    for k in sorted(result["summary"].keys()):
        typer.echo(f"{k}\t{result['summary'][k]}")
    if result.get("per_contig_written"):
        typer.echo(f"# per-contig TSV written: {result['per_contig_written']}")
