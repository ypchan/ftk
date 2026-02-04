from __future__ import annotations

import argparse
from pathlib import Path

from rich_argparse import RichHelpFormatter

from ftk.core.fishing import FishingConfig, run_fishing_pipeline


def add_fishing_parser(subparsers: argparse._SubParsersAction) -> None:
    desc = (
        "Extract homologous sequences from genome pools using minimap2.\n\n"
        "Baits (TSV):\n"
        "  label<TAB>sequence\n"
        "  - Multiple lines with the same label represent multiple bait sequences for the same hook.\n\n"
        "Pools:\n"
        "  Either provide --pool (one or more FASTA/FASTA.gz files) OR --pool-list (TSV):\n"
        "    pool_id<TAB>path\n"
        "  - Empty lines and lines starting with '#' are ignored.\n"
        "  - Relative paths in --pool-list are resolved relative to the list file.\n\n"
        "Output (pure text, not Rich):\n"
        "  With --out-prefix PREFIX, writes:\n"
        "    PREFIX.tsv           (full results: one row per pool_id × label)\n"
        "    PREFIX.fasta         (PASS sequences only; oriented to bait direction; with flanks)\n"
        "    PREFIX.presence.tsv  (presence-only long table: pool_id, label, present 0/1)\n"
        "  Optionally:\n"
        "    --save-paf           saves raw PAF to PREFIX.paf (or user-specified path)\n"
    )

    p = subparsers.add_parser(
        "fishing",
        help="Fish homologous sequences from genome pools using minimap2.",
        description=desc,
        formatter_class=RawTextRichHelpFormatter,
    )

    p.add_argument(
        "--baits",
        required=True,
        help="Baits TSV file (label<TAB>sequence). Use '-' to read from stdin.",
    )

    mx = p.add_mutually_exclusive_group(required=True)
    mx.add_argument(
        "--pool",
        nargs="+",
        help="One or more genome FASTA/FASTA.gz files; each file is treated as one pool.",
    )
    mx.add_argument(
        "--pool-list",
        help="TSV file with two columns: pool_id<TAB>path (one pool per line).",
    )

    p.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix path. Writes PREFIX.tsv, PREFIX.fasta, PREFIX.presence.tsv by default.",
    )

    p.add_argument(
        "--preset",
        choices=["asm5", "asm10", "asm20"],
        default="asm20",
        help="minimap2 preset for DNA contigs/scaffolds.",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Number of threads for minimap2. 0 means use all CPUs.",
    )

    p.add_argument(
        "--min-identity",
        type=float,
        default=0.50,
        help="Minimum identity (matches/aln_len) to accept a hit.",
    )
    p.add_argument(
        "--min-aln-len",
        type=int,
        default=50,
        help="Minimum alignment block length (bp).",
    )
    p.add_argument(
        "--min-aln-frac",
        type=float,
        default=0.50,
        help="Minimum aligned fraction relative to bait length (aln_len / qlen).",
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum minimap2 MAPQ to accept a hit.",
    )

    p.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Flank length in bp on both bait-forward and bait-reverse sides.",
    )
    p.add_argument(
        "--flank-f",
        type=int,
        default=None,
        help="Flank length in bp on the bait-forward side (5' upstream in bait orientation).",
    )
    p.add_argument(
        "--flank-r",
        type=int,
        default=None,
        help="Flank length in bp on the bait-reverse side (3' downstream in bait orientation).",
    )

    p.add_argument(
        "--save-paf",
        nargs="?",
        const=True,
        default=False,
        help="Save raw PAF. If provided without a value, saves to PREFIX.paf.",
    )

    p.add_argument(
        "--minimap2",
        default="minimap2",
        help="Path to minimap2 binary (default: minimap2).",
    )

    p.add_argument(
        "--tmpdir",
        type=Path,
        default=None,
        help="Temporary directory for intermediate files (default: system temp).",
    )


def run_fishing(args: argparse.Namespace) -> None:
    cfg = FishingConfig(
        baits_path=args.baits,
        pool_paths=args.pool,
        pool_list_path=args.pool_list,
        out_prefix=Path(args.out_prefix),
        preset=args.preset,
        threads=args.threads,
        min_identity=args.min_identity,
        min_aln_len=args.min_aln_len,
        min_aln_frac=args.min_aln_frac,
        min_mapq=args.min_mapq,
        flank=args.flank,
        flank_f=args.flank_f,
        flank_r=args.flank_r,
        save_paf=args.save_paf,
        minimap2=args.minimap2,
        tmpdir=args.tmpdir,
    )
    run_fishing_pipeline(cfg)
