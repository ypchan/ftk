from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from ftk.core.fasta import iter_fasta, iter_fasta_stdin, FastaRecord
from ftk.core.gc_peaks import count_gc_peaks, GcPeaksMethod

IUPAC_SET: Set[str] = set("ACGTNRYKMSWBDHV")


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


ALL_METRICS: List[Metric] = [
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


@dataclass
class PerContig:
    source: str
    name: str
    length: int
    n_count: int
    unknown_count: int
    gc: Optional[float]
    contig_or_scaffold: str


@dataclass
class AssemblyResult:
    source: str
    record_count: int
    genome_size_with_N: int
    genome_size_ungapped: int
    total_N: int
    total_unknown: int
    longest: int
    shortest: int
    global_gc: Optional[float]
    gc_avg: Optional[float]
    gc_max: Optional[float]
    gc_min: Optional[float]
    n50: Optional[int]
    l50: Optional[int]
    n90: Optional[int]
    l90: Optional[int]
    gc_peaks: Optional[int]
    per_contig: Optional[List[PerContig]]


def compute_assembly_stats(
    *,
    source: str,
    fasta_path: Optional[Path],
    use_stdin: bool,
    n_run_scaffold_threshold: int,
    metrics: Optional[List[Metric]],
    gc_peaks_method: GcPeaksMethod,
    gc_kde_grid: int,
    gc_kde_min_prominence: float,
    gc_gmm_max_components: int,
) -> AssemblyResult:
    wanted = set(metrics) if metrics else set(ALL_METRICS)

    records = _iter_records(fasta_path=fasta_path, use_stdin=use_stdin)

    per_contig_list: Optional[List[PerContig]] = []
    lengths: List[int] = []

    total_len = 0
    total_N = 0
    total_unknown = 0
    longest = 0
    shortest = 0
    record_count = 0

    gc_num = 0
    gc_den = 0
    per_record_gc: List[float] = []

    for rec in records:
        record_count += 1
        seq = rec.seq
        L = len(seq)
        lengths.append(L)

        total_len += L
        longest = max(longest, L)
        shortest = L if shortest == 0 else min(shortest, L)

        n_count = seq.count("N")
        total_N += n_count

        unknown = _count_unknown(seq)
        total_unknown += unknown

        gc_val, g_num, g_den = _gc_stats(seq)
        if g_den > 0:
            gc_num += g_num
            gc_den += g_den
        if gc_val is not None:
            per_record_gc.append(gc_val)

        c_or_s = _classify_contig_or_scaffold(seq, threshold=n_run_scaffold_threshold)

        per_contig_list.append(
            PerContig(
                source=source,
                name=rec.name,
                length=L,
                n_count=n_count,
                unknown_count=unknown,
                gc=gc_val,
                contig_or_scaffold=c_or_s,
            )
        )

    genome_size_with_N = total_len
    genome_size_ungapped = total_len - total_N

    global_gc = (gc_num / gc_den) if gc_den > 0 else None

    if per_record_gc:
        gc_avg = sum(per_record_gc) / len(per_record_gc)
        gc_max = max(per_record_gc)
        gc_min = min(per_record_gc)
    else:
        gc_avg = gc_max = gc_min = None

    n50 = l50 = n90 = l90 = None
    if lengths and Metric.n_stats in wanted:
        n50, l50 = _nx_lx(lengths, genome_size_with_N, x=50)
        n90, l90 = _nx_lx(lengths, genome_size_with_N, x=90)

    gc_peaks_val: Optional[int] = None
    if Metric.gc_peaks in wanted:
        gc_peaks_val = count_gc_peaks(
            per_record_gc,
            gc_peaks_method,
            kde_grid=gc_kde_grid,
            kde_min_prominence=gc_kde_min_prominence,
            gmm_max_components=gc_gmm_max_components,
        )

    return AssemblyResult(
        source=source,
        record_count=record_count,
        genome_size_with_N=genome_size_with_N,
        genome_size_ungapped=genome_size_ungapped,
        total_N=total_N,
        total_unknown=total_unknown,
        longest=longest,
        shortest=shortest,
        global_gc=global_gc,
        gc_avg=gc_avg,
        gc_max=gc_max,
        gc_min=gc_min,
        n50=n50,
        l50=l50,
        n90=n90,
        l90=l90,
        gc_peaks=gc_peaks_val,
        per_contig=per_contig_list,
    )


def _iter_records(*, fasta_path: Optional[Path], use_stdin: bool) -> Iterable[FastaRecord]:
    if use_stdin:
        return iter_fasta_stdin()
    if fasta_path is None:
        raise ValueError("fasta_path is required when use_stdin is False")
    return iter_fasta(fasta_path)


def _count_unknown(seq: str) -> int:
    return sum(1 for ch in seq if ch not in IUPAC_SET)


def _gc_stats(seq: str) -> tuple[Optional[float], int, int]:
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    den = a + c + g + t
    num = g + c
    if den == 0:
        return None, 0, 0
    return (num / den), num, den


def _classify_contig_or_scaffold(seq: str, *, threshold: int) -> str:
    run = 0
    for ch in seq:
        if ch == "N":
            run += 1
            if run >= threshold:
                return "scaffold"
        else:
            run = 0
    return "contig"


def _nx_lx(lengths: List[int], genome_size: int, *, x: int) -> tuple[int, int]:
    target = genome_size * (x / 100.0)
    s = 0
    for i, L in enumerate(sorted(lengths, reverse=True), start=1):
        s += L
        if s >= target:
            return L, i
    return 0, len(lengths)


def format_tsv_rows(results: List[AssemblyResult], *, metrics: Optional[List[Metric]]) -> str:
    wanted = metrics if metrics else ALL_METRICS

    cols: List[str] = ["source", "record_count"]

    for m in wanted:
        if m == Metric.genome_size:
            cols.extend(["genome_size_with_N", "genome_size_ungapped", "total_N"])
        elif m == Metric.contig_or_scaffold:
            cols.extend(["contig_count", "scaffold_count"])
        elif m == Metric.unknown_count:
            cols.append("unknown_count")
        elif m == Metric.longest:
            cols.append("longest")
        elif m == Metric.shortest:
            cols.append("shortest")
        elif m == Metric.gc:
            cols.append("gc")
        elif m == Metric.n_stats:
            cols.extend(["N50", "L50", "N90", "L90"])
        elif m == Metric.gc_avg:
            cols.append("gc_avg")
        elif m == Metric.gc_max:
            cols.append("gc_max")
        elif m == Metric.gc_min:
            cols.append("gc_min")
        elif m == Metric.gc_peaks:
            cols.append("gc_peaks")

    seen = set()
    cols = [c for c in cols if not (c in seen or seen.add(c))]

    lines = ["\t".join(cols)]
    for r in results:
        row = _row_dict(r)
        lines.append("\t".join("" if row.get(c) is None else str(row.get(c)) for c in cols))
    return "\n".join(lines) + "\n"


def _row_dict(r: AssemblyResult) -> Dict[str, object]:
    contig_count = 0
    scaffold_count = 0
    if r.per_contig:
        for pc in r.per_contig:
            if pc.contig_or_scaffold == "scaffold":
                scaffold_count += 1
            else:
                contig_count += 1

    return {
        "source": r.source,
        "record_count": r.record_count,
        "genome_size_with_N": r.genome_size_with_N,
        "genome_size_ungapped": r.genome_size_ungapped,
        "total_N": r.total_N,
        "contig_count": contig_count,
        "scaffold_count": scaffold_count,
        "unknown_count": r.total_unknown,
        "longest": r.longest,
        "shortest": r.shortest,
        "gc": None if r.global_gc is None else round(r.global_gc, 6),
        "N50": r.n50,
        "L50": r.l50,
        "N90": r.n90,
        "L90": r.l90,
        "gc_avg": None if r.gc_avg is None else round(r.gc_avg, 6),
        "gc_max": None if r.gc_max is None else round(r.gc_max, 6),
        "gc_min": None if r.gc_min is None else round(r.gc_min, 6),
        "gc_peaks": r.gc_peaks,
    }


def per_contig_tsv_rows(r: AssemblyResult) -> List[Dict[str, object]]:
    if not r.per_contig:
        return []
    rows: List[Dict[str, object]] = []
    for pc in r.per_contig:
        rows.append(
            {
                "source": pc.source,
                "name": pc.name,
                "length": pc.length,
                "N_count": pc.n_count,
                "unknown_count": pc.unknown_count,
                "gc": None if pc.gc is None else round(pc.gc, 6),
                "contig_or_scaffold": pc.contig_or_scaffold,
            }
        )
    return rows
