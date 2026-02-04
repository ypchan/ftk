from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Set, Tuple, Optional, Dict
import gzip
import os
import shutil
import subprocess
import tempfile

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, fcluster


@dataclass(frozen=True)
class MsaQcConfig:
    msa_path: Path
    out_prefix: Path

    moltype: str = "auto"  # auto|dna|protein
    gap_chars: Set[str] = None  # set of gap characters

    min_overlap_frac: float = 0.20
    min_mean_pid: float = 0.50
    min_nn_pid: float = 0.60
    rz_cutoff: float = 3.5

    cluster_cutoff: str = "auto"  # "auto" or numeric str
    max_heatmap: int = 200
    no_heatmap: bool = False
    no_plots: bool = False

    re_align: str = "off"  # off|fix
    mafft: str = "mafft"
    threads: int = 0


@dataclass
class TaxonRow:
    taxon: str
    moltype: str
    aligned_len: int
    nongap_frac: float
    mean_overlap_frac: float
    mean_pid: float
    nn_pid: float
    mean_dist: float
    in_main_cluster: int
    is_anomalous: int
    flags: str


def run_msaqc_pipeline(cfg: MsaQcConfig) -> None:
    if cfg.gap_chars is None:
        gap_chars = set("-.")
    else:
        gap_chars = set(c.upper() for c in cfg.gap_chars)

    msa_path = cfg.msa_path

    # Optional realignment
    if cfg.re_align == "fix":
        _ensure_binary(cfg.mafft, what="MAFFT (requested by --re-align fix)")
        msa_path = _run_mafft_realignment(cfg, msa_path)

    taxa, seqs = _read_fasta(msa_path)
    if not taxa:
        raise SystemExit(f"ERROR: empty MSA: {msa_path}")

    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise SystemExit("ERROR: MSA sequences must have identical aligned length.")

    moltype = _infer_moltype(cfg.moltype, seqs, gap_chars)

    msa = _encode_msa(seqs, gap_chars)

    mean_overlap, mean_pid, nn_pid, mean_dist, dist_mat = _pairwise_stats(msa)

    nongap = (msa != 0).sum(axis=1)
    nongap_frac = nongap / float(L)

    # Cluster and main cluster detection
    cutoff = _resolve_cluster_cutoff(cfg.cluster_cutoff, nn_pid)
    in_main = _main_cluster_membership(dist_mat, cutoff)

    # Robust z-scores
    rz_overlap = _robust_z(mean_overlap)
    rz_pid = _robust_z(mean_pid)

    rows: List[TaxonRow] = []
    for i, t in enumerate(taxa):
        flags: List[str] = []

        if mean_overlap[i] < cfg.min_overlap_frac:
            flags.append("LOW_OVERLAP_HARD")
        if mean_pid[i] < cfg.min_mean_pid:
            flags.append("LOW_PID_HARD")
        if nn_pid[i] < cfg.min_nn_pid:
            flags.append("LOW_NN_PID_HARD")

        if rz_overlap[i] < -cfg.rz_cutoff:
            flags.append("LOW_OVERLAP_RZ")
        if rz_pid[i] < -cfg.rz_cutoff:
            flags.append("LOW_PID_RZ")

        if in_main[i] == 0:
            flags.append("OUTSIDE_MAIN_CLUSTER")

        # Heuristic "paralogy/contamination": aligns but isolated and low similarity
        if (
            mean_overlap[i] >= max(cfg.min_overlap_frac, 0.15)
            and mean_pid[i] < cfg.min_mean_pid
            and nn_pid[i] < cfg.min_nn_pid
            and in_main[i] == 0
        ):
            flags.append("POSSIBLE_PARALOGY_OR_CONTAM")

        is_anom = 1 if flags else 0

        rows.append(
            TaxonRow(
                taxon=t,
                moltype=moltype,
                aligned_len=L,
                nongap_frac=float(nongap_frac[i]),
                mean_overlap_frac=float(mean_overlap[i]),
                mean_pid=float(mean_pid[i]),
                nn_pid=float(nn_pid[i]),
                mean_dist=float(mean_dist[i]),
                in_main_cluster=int(in_main[i]),
                is_anomalous=is_anom,
                flags=",".join(flags),
            )
        )

    cfg.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(cfg.out_prefix.with_suffix(".tsv"), rows)

    if not cfg.no_plots:
        _plot_scatter(cfg.out_prefix.with_suffix(".scatter.png"), rows)
        if (not cfg.no_heatmap) and len(rows) <= cfg.max_heatmap:
            _plot_heatmap(cfg.out_prefix.with_suffix(".heatmap.png"), taxa, dist_mat)


def _ensure_binary(bin_path: str, *, what: str) -> None:
    if shutil.which(bin_path) is None and not Path(bin_path).exists():
        raise SystemExit(
            f"ERROR: {what} not found: '{bin_path}'.\n"
            "Explanation: You requested functionality that requires this external program.\n"
            "Fix: Install it (e.g., conda/pixi/apt) or provide a valid path via --mafft.\n"
        )


def _threads_value(threads: int) -> int:
    if threads is None or threads <= 0:
        return max(1, os.cpu_count() or 1)
    return threads


def _run_mafft_realignment(cfg: MsaQcConfig, msa_path: Path) -> Path:
    # MAFFT expects unaligned FASTA too; but we allow aligned input (it will realign).
    threads = _threads_value(cfg.threads)

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        out_msa = td_path / "realigned.fasta"

        cmd = [cfg.mafft, "--thread", str(threads), str(msa_path)]
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(
                "ERROR: MAFFT failed.\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDERR:\n{proc.stderr}"
            )
        out_msa.write_text(proc.stdout, encoding="utf-8")
        # Persist to a named output file next to out_prefix for reproducibility
        fixed_path = cfg.out_prefix.with_suffix(".fixed.fasta")
        fixed_path.write_text(proc.stdout, encoding="utf-8")
        return fixed_path


def _read_fasta(path: Path) -> Tuple[List[str], List[str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    taxa: List[str] = []
    seqs: List[str] = []

    name = None
    chunks: List[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(chunks).upper())
                name = line[1:].split()[0]
                taxa.append(name)
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs.append("".join(chunks).upper())

    return taxa, seqs


def _infer_moltype(moltype: str, seqs: Sequence[str], gap_chars: Set[str]) -> str:
    if moltype in ("dna", "protein"):
        return moltype

    # auto: count fraction of DNA-like letters among non-gap letters
    dna_letters = set("ACGTUNRYKMSWBDHV")  # allow U in some data
    total = 0
    dna_like = 0
    for s in seqs:
        for ch in s:
            c = ch.upper()
            if c in gap_chars:
                continue
            total += 1
            if c in dna_letters:
                dna_like += 1
    if total == 0:
        return "dna"
    frac = dna_like / total
    # If overwhelmingly DNA-like, call DNA; else protein
    return "dna" if frac >= 0.90 else "protein"


def _encode_msa(seqs: Sequence[str], gap_chars: Set[str]) -> np.ndarray:
    n = len(seqs)
    L = len(seqs[0])
    mat = np.zeros((n, L), dtype=np.uint8)

    for i, s in enumerate(seqs):
        b = np.frombuffer(s.encode("ascii", errors="replace"), dtype=np.uint8).copy()
        mask = np.zeros(L, dtype=bool)
        for gc in gap_chars:
            mask |= (b == ord(gc))
        b[mask] = 0
        mat[i, :] = b
    return mat


def _pairwise_stats(msa: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n, L = msa.shape
    non_gap = (msa != 0)

    overlap_sum = np.zeros(n, dtype=float)
    pid_sum = np.zeros(n, dtype=float)
    dist_sum = np.zeros(n, dtype=float)
    count = np.zeros(n, dtype=float)
    nn_pid = np.full(n, np.nan, dtype=float)

    dist_mat = np.zeros((n, n), dtype=float)

    for i in tqdm(range(n), desc="Pairwise stats", unit="taxon"):
        best_pid_i = -1.0
        for j in range(i + 1, n):
            ov_mask = non_gap[i] & non_gap[j]
            ov = int(ov_mask.sum())
            if ov == 0:
                pid = 0.0
                dist = 1.0
                ov_frac = 0.0
            else:
                matches = int((msa[i, ov_mask] == msa[j, ov_mask]).sum())
                pid = matches / ov
                dist = 1.0 - pid
                ov_frac = ov / float(L)

            overlap_sum[i] += ov_frac
            overlap_sum[j] += ov_frac
            pid_sum[i] += pid
            pid_sum[j] += pid
            dist_sum[i] += dist
            dist_sum[j] += dist
            count[i] += 1.0
            count[j] += 1.0

            dist_mat[i, j] = dist
            dist_mat[j, i] = dist

            if pid > best_pid_i:
                best_pid_i = pid
            if np.isnan(nn_pid[j]) or pid > nn_pid[j]:
                nn_pid[j] = pid

        if np.isnan(nn_pid[i]):
            nn_pid[i] = best_pid_i if best_pid_i >= 0 else 0.0
        else:
            nn_pid[i] = max(nn_pid[i], best_pid_i if best_pid_i >= 0 else 0.0)

    mean_overlap = np.divide(overlap_sum, count, out=np.zeros_like(overlap_sum), where=count > 0)
    mean_pid = np.divide(pid_sum, count, out=np.zeros_like(pid_sum), where=count > 0)
    mean_dist = np.divide(dist_sum, count, out=np.zeros_like(dist_sum), where=count > 0)

    if n == 1:
        mean_overlap[:] = 1.0
        mean_pid[:] = 1.0
        nn_pid[:] = 1.0
        mean_dist[:] = 0.0

    return mean_overlap, mean_pid, nn_pid, mean_dist, dist_mat


def _robust_z(x: np.ndarray) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.zeros_like(x)
    return (x - med) / (1.4826 * mad)


def _resolve_cluster_cutoff(cluster_cutoff: str, nn_pid: np.ndarray) -> float:
    if cluster_cutoff != "auto":
        try:
            t = float(cluster_cutoff)
        except ValueError:
            raise SystemExit("ERROR: --cluster-cutoff must be 'auto' or a numeric value in [0,1].")
        if t < 0 or t > 1:
            raise SystemExit("ERROR: --cluster-cutoff numeric value must be in [0,1].")
        return t

    nn_dist = 1.0 - nn_pid
    med = np.median(nn_dist)
    mad = np.median(np.abs(nn_dist - med))
    if mad == 0:
        t = med
    else:
        t = med + 3.0 * 1.4826 * mad
    # Clamp for stability
    t = min(max(t, 0.05), 0.50)
    return float(t)


def _main_cluster_membership(dist_mat: np.ndarray, cutoff: float) -> np.ndarray:
    n = dist_mat.shape[0]
    if n == 1:
        return np.array([1], dtype=int)

    # linkage expects condensed distance vector
    triu = dist_mat[np.triu_indices(n, k=1)]
    Z = linkage(triu, method="average")
    labels = fcluster(Z, t=cutoff, criterion="distance")  # cluster labels start at 1

    # main cluster = largest size
    uniq, counts = np.unique(labels, return_counts=True)
    main_label = uniq[np.argmax(counts)]
    in_main = (labels == main_label).astype(int)
    return in_main


def _write_tsv(path: Path, rows: Sequence[TaxonRow]) -> None:
    cols = [
        "taxon",
        "moltype",
        "aligned_len",
        "nongap_frac",
        "mean_overlap_frac",
        "mean_pid",
        "nn_pid",
        "mean_dist",
        "in_main_cluster",
        "is_anomalous",
        "flags",
    ]
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write(
                "\t".join(
                    [
                        r.taxon,
                        r.moltype,
                        str(r.aligned_len),
                        f"{r.nongap_frac:.6f}",
                        f"{r.mean_overlap_frac:.6f}",
                        f"{r.mean_pid:.6f}",
                        f"{r.nn_pid:.6f}",
                        f"{r.mean_dist:.6f}",
                        str(r.in_main_cluster),
                        str(r.is_anomalous),
                        r.flags,
                    ]
                )
                + "\n"
            )


def _plot_scatter(path: Path, rows: Sequence[TaxonRow]) -> None:
    x = np.array([r.mean_overlap_frac for r in rows])
    y = np.array([r.mean_pid for r in rows])
    flagged = np.array([r.is_anomalous == 1 for r in rows], dtype=bool)

    plt.figure()
    plt.scatter(x[~flagged], y[~flagged], s=20)
    plt.scatter(x[flagged], y[flagged], s=35)
    for r in rows:
        if r.is_anomalous == 1:
            plt.text(r.mean_overlap_frac, r.mean_pid, r.taxon, fontsize=8)

    plt.xlabel("mean_overlap_frac")
    plt.ylabel("mean_pid")
    plt.title("MSA QC: overlap vs identity (flagged taxa labeled)")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def _plot_heatmap(path: Path, taxa: Sequence[str], dist: np.ndarray) -> None:
    plt.figure(figsize=(8, 7))
    plt.imshow(dist, aspect="auto")
    plt.colorbar(label="distance (1 - pid on overlapped sites)")
    plt.title("Pairwise distance heatmap")
    if len(taxa) <= 50:
        plt.xticks(range(len(taxa)), taxa, rotation=90, fontsize=6)
        plt.yticks(range(len(taxa)), taxa, fontsize=6)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
