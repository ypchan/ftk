from __future__ import annotations

import argparse
from pathlib import Path

from rich_argparse import RichRawDescriptionHelpFormatter
from ftk.core.msaqc import MsaQcConfig, run_msaqc_pipeline


def add_msaqc_parser(subparsers: argparse._SubParsersAction) -> None:
    desc = (
        "Detect anomalous taxa in a barcode/marker MSA (DNA or protein).\n\n"
        "Input:\n"
        "  --msa PATH\n"
        "    Aligned FASTA/FASTA.gz. All sequences must have the same aligned length.\n"
        "    Gap characters are excluded from overlap/identity calculations.\n\n"
        "Key metrics (per taxon):\n"
        "  nongap_frac\n"
        "    Fraction of aligned positions that are NOT gaps for this taxon.\n"
        "  mean_overlap_frac\n"
        "    Mean fraction of aligned positions that are non-gap in BOTH sequences (averaged over all other taxa).\n"
        "    Low values often indicate partial sequences, trimming, primer mismatch, or alignment artifacts.\n"
        "  mean_pid\n"
        "    Mean pairwise identity on overlapped (non-gap in both) sites.\n"
        "  nn_pid\n"
        "    Nearest-neighbor identity (best match among other taxa) on overlapped sites.\n"
        "  mean_dist\n"
        "    Mean distance, defined as 1 - pid, averaged over other taxa.\n\n"
        "Main-cluster detection (SciPy):\n"
        "  We build a distance matrix (dist = 1 - pid_on_overlap) and perform hierarchical clustering (average linkage).\n"
        "  The main cluster is defined as the largest cluster after cutting the tree at --cluster-cutoff.\n"
        "  Taxa outside the main cluster are flagged as OUTSIDE_MAIN_CLUSTER.\n\n"
        "Cluster cutoff:\n"
        "  --cluster-cutoff auto (default)\n"
        "    Uses robust statistics of nearest-neighbor distances:\n"
        "      nn_dist = 1 - nn_pid\n"
        "      cutoff = median(nn_dist) + 3 * 1.4826 * MAD(nn_dist)\n"
        "    Then clamps cutoff to [0.05, 0.50] for stability.\n\n"
        "Outputs (pure files; Rich affects help only):\n"
        "  --out-prefix PREFIX writes:\n"
        "    PREFIX.tsv          Per-taxon QC table with is_anomalous and flags.\n"
        "    PREFIX.scatter.png  Scatter: mean_overlap_frac vs mean_pid (flagged taxa labeled).\n"
        "    PREFIX.heatmap.png  Distance heatmap (if taxa <= --max-heatmap and not disabled).\n\n"
        "Optional realignment:\n"
        "  --re-align fix --mafft mafft\n"
        "    Re-aligns sequences with MAFFT and re-runs QC. If MAFFT is requested but not found, an error is raised.\n"
        "    Note: default is no realignment.\n"
    )

    p = subparsers.add_parser(
        "msaqc",
        help="QC an MSA: detect anomalous taxa; auto main cluster; plots + TSV.",
        description=desc,
        formatter_class=RawDescriptionRichHelpFormatter,
    )

    p.add_argument("--msa", required=True, help="Input aligned MSA FASTA/FASTA.gz.")
    p.add_argument("--out-prefix", required=True, help="Output prefix path (writes TSV + PNG plots).")

    p.add_argument(
        "--moltype",
        choices=["auto", "dna", "protein"],
        default="auto",
        help="Molecule type. auto infers from characters (default: auto).",
    )
    p.add_argument(
        "--gap-chars",
        default="-.",
        help="Characters treated as gaps (default: '-.').",
    )

    p.add_argument("--min-overlap-frac", type=float, default=0.20, help="Hard threshold: mean_overlap_frac < this => LOW_OVERLAP_HARD.")
    p.add_argument("--min-mean-pid", type=float, default=0.50, help="Hard threshold: mean_pid < this => LOW_PID_HARD.")
    p.add_argument("--min-nn-pid", type=float, default=0.60, help="Hard threshold: nn_pid < this => LOW_NN_PID_HARD.")
    p.add_argument("--rz-cutoff", type=float, default=3.5, help="Robust z-score cutoff using MAD (default: 3.5).")

    p.add_argument(
        "--cluster-cutoff",
        default="auto",
        help="Cluster cutoff for main-cluster detection. Use 'auto' or a numeric distance in [0,1].",
    )

    p.add_argument("--max-heatmap", type=int, default=200, help="Max taxa to render heatmap (default: 200).")
    p.add_argument("--no-heatmap", action="store_true", help="Disable heatmap generation.")
    p.add_argument("--no-plots", action="store_true", help="Disable all plot generation.")

    p.add_argument(
        "--re-align",
        choices=["off", "fix"],
        default="off",
        help="Optional: re-align with MAFFT then re-run QC (default: off).",
    )
    p.add_argument("--mafft", default="mafft", help="Path to MAFFT binary (default: mafft). Only used if --re-align fix.")
    p.add_argument("--threads", type=int, default=0, help="Threads for MAFFT (if used). 0 means use all CPUs.")


def run_msaqc(args: argparse.Namespace) -> None:
    cfg = MsaQcConfig(
        msa_path=Path(args.msa),
        out_prefix=Path(args.out_prefix),
        moltype=args.moltype,
        gap_chars=set(args.gap_chars),
        min_overlap_frac=args.min_overlap_frac,
        min_mean_pid=args.min_mean_pid,
        min_nn_pid=args.min_nn_pid,
        rz_cutoff=args.rz_cutoff,
        cluster_cutoff=args.cluster_cutoff,
        max_heatmap=args.max_heatmap,
        no_heatmap=args.no_heatmap,
        no_plots=args.no_plots,
        re_align=args.re_align,
        mafft=args.mafft,
        threads=args.threads,
    )
    run_msaqc_pipeline(cfg)
