import os
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List, Optional, Set, Tuple

from ftk.core.fasta import iter_fasta_records
from ftk.core.gc_peaks import gc_peaks_kde, gc_peaks_gmm


_IUPAC = set("ACGTNRYKMSWBDHV")


@dataclass
class ContigStat:
    name: str
    length: int
    n_count: int
    unknown_count: int
    gc: Optional[float]
    is_scaffold: bool
    a: int
    c: int
    g: int
    t: int


def compute_contig_stat(name: str, seq: str, n_run_scaffold_threshold: int) -> ContigStat:
    length = len(seq)
    n_count = 0
    unknown_count = 0
    a = c = g = t = 0
    max_n_run = 0
    cur_n_run = 0

    for ch in seq:
        if ch == "N":
            n_count += 1
            cur_n_run += 1
            if cur_n_run > max_n_run:
                max_n_run = cur_n_run
        else:
            cur_n_run = 0
            if ch == "A":
                a += 1
            elif ch == "C":
                c += 1
            elif ch == "G":
                g += 1
            elif ch == "T":
                t += 1

        if ch not in _IUPAC:
            unknown_count += 1

    denom = a + c + g + t
    gc = None if denom == 0 else (g + c) / denom

    return ContigStat(
        name=name,
        length=length,
        n_count=n_count,
        unknown_count=unknown_count,
        gc=gc,
        is_scaffold=(max_n_run >= n_run_scaffold_threshold),
        a=a, c=c, g=g, t=t
    )


def _worker(args: Tuple[str, str, int]) -> ContigStat:
    name, seq, thr = args
    return compute_contig_stat(name, seq, thr)


def calc_nxx_lxx(lengths: List[int], total_size: int, xx: int) -> Tuple[int, int]:
    if not lengths or total_size <= 0:
        return 0, 0
    target = total_size * (xx / 100.0)
    cum = 0
    for i, L in enumerate(sorted(lengths, reverse=True), start=1):
        cum += L
        if cum >= target:
            return L, i
    return sorted(lengths, reverse=True)[-1], len(lengths)


def run_stat(
    fasta: str,
    metrics: Set,
    threads: int,
    n_run_scaffold_threshold: int,
    per_contig_path: Optional[str],
    gc_peaks_method: str,
    gc_kde_grid: int,
    gc_kde_min_prominence: float,
    gc_gmm_max_components: int,
) -> Dict:
    if not os.path.exists(fasta):
        raise FileNotFoundError(fasta)

    cpu = os.cpu_count() or 1
    if threads <= 0:
        threads = cpu

    lengths: List[int] = []
    gc_list: List[float] = []

    scaffold_count = 0
    contig_count = 0
    unknown_total = 0
    n_total = 0
    genome_size_with_N = 0

    A = C = G = T = 0

    tsv = None
    if per_contig_path:
        tsv = open(per_contig_path, "w", encoding="utf-8")
        tsv.write("name\tlength\tN_count\tunknown_count\tgc\tcontig_or_scaffold\n")

    record_iter = iter_fasta_records(fasta)

    if threads == 1:
        for name, seq in record_iter:
            cs = compute_contig_stat(name, seq, n_run_scaffold_threshold)
            lengths.append(cs.length)
            genome_size_with_N += cs.length
            n_total += cs.n_count
            unknown_total += cs.unknown_count
            if cs.is_scaffold:
                scaffold_count += 1
            else:
                contig_count += 1
            if cs.gc is not None:
                gc_list.append(cs.gc)
            A += cs.a; C += cs.c; G += cs.g; T += cs.t

            if tsv:
                gc_s = "NA" if cs.gc is None else f"{cs.gc:.6f}"
                tsv.write(f"{cs.name}\t{cs.length}\t{cs.n_count}\t{cs.unknown_count}\t"
                          f"{gc_s}\t{'scaffold' if cs.is_scaffold else 'contig'}\n")
    else:
        with Pool(processes=threads) as pool:
            args_iter = ((name, seq, n_run_scaffold_threshold) for name, seq in record_iter)
            for cs in pool.imap_unordered(_worker, args_iter, chunksize=64):
                lengths.append(cs.length)
                genome_size_with_N += cs.length
                n_total += cs.n_count
                unknown_total += cs.unknown_count
                if cs.is_scaffold:
                    scaffold_count += 1
                else:
                    contig_count += 1
                if cs.gc is not None:
                    gc_list.append(cs.gc)
                A += cs.a; C += cs.c; G += cs.g; T += cs.t

                if tsv:
                    gc_s = "NA" if cs.gc is None else f"{cs.gc:.6f}"
                    tsv.write(f"{cs.name}\t{cs.length}\t{cs.n_count}\t{cs.unknown_count}\t"
                              f"{gc_s}\t{'scaffold' if cs.is_scaffold else 'contig'}\n")

    if tsv:
        tsv.close()

    if not lengths:
        raise ValueError("No sequences found in FASTA.")

    genome_size_ungapped = genome_size_with_N - n_total
    denom = A + C + G + T
    gc_global = None if denom == 0 else (G + C) / denom

    n50, l50 = calc_nxx_lxx(lengths, genome_size_with_N, 50)
    n90, l90 = calc_nxx_lxx(lengths, genome_size_with_N, 90)

    gc_avg = (sum(gc_list) / len(gc_list)) if gc_list else None
    gc_min = min(gc_list) if gc_list else None
    gc_max = max(gc_list) if gc_list else None

    if gc_peaks_method == "kde":
        peaks = gc_peaks_kde(gc_list, grid=gc_kde_grid, min_prominence=gc_kde_min_prominence) if gc_list else 0
    else:
        peaks = gc_peaks_gmm(gc_list, max_components=gc_gmm_max_components) if gc_list else 0

    def fmt(x: Optional[float]) -> str:
        return "NA" if x is None else f"{x:.6f}"

    summary: Dict[str, str] = {}
    # Output keys are stable for downstream parsing.
    if "genome_size" in metrics or getattr(metrics, "__contains__", lambda x: False)("genome_size"):
        summary["genome_size_with_N(bp)"] = str(genome_size_with_N)
        summary["genome_size_ungapped(bp)"] = str(genome_size_ungapped)
        summary["total_N_count"] = str(n_total)

    if "contig_or_scaffold" in metrics:
        summary["contig_count"] = str(contig_count)
        summary["scaffold_count"] = str(scaffold_count)
        summary["scaffold_rule"] = f"N_run >= {n_run_scaffold_threshold}"

    if "unknown_count" in metrics:
        summary["unknown_char_count"] = str(unknown_total)

    if "longest" in metrics:
        summary["longest_record(bp)"] = str(max(lengths))

    if "shortest" in metrics:
        summary["shortest_record(bp)"] = str(min(lengths))

    if "gc" in metrics:
        summary["gc_content_global"] = fmt(gc_global)

    if "n_stats" in metrics:
        summary["N50(bp)"] = str(n50)
        summary["L50"] = str(l50)
        summary["N90(bp)"] = str(n90)
        summary["L90"] = str(l90)

    if "gc_avg" in metrics:
        summary["gc_avg(per_record_mean)"] = fmt(gc_avg)

    if "gc_max" in metrics:
        summary["gc_max(per_record)"] = fmt(gc_max)

    if "gc_min" in metrics:
        summary["gc_min(per_record)"] = fmt(gc_min)

    if "gc_peaks" in metrics:
        summary["gc_peaks"] = str(peaks)
        summary["gc_peaks_method"] = gc_peaks_method

    return {
        "input": os.path.basename(fasta),
        "summary": summary,
        "per_contig_written": per_contig_path if per_contig_path else None,
    }
