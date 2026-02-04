from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable
import os
import shutil
import subprocess
import tempfile

from tqdm import tqdm

from ftk.core.paf import parse_paf_lines, PafRecord
from ftk.core.fasta_reader import FastaInMemory, reverse_complement


@dataclass(frozen=True)
class FishingConfig:
    baits_path: str  # file path or '-'
    pool_paths: Optional[List[str]]  # from --pool
    pool_list_path: Optional[str]  # from --pool-list (id\tpath)
    out_prefix: Path

    preset: str = "asm20"
    threads: int = 0

    min_identity: float = 0.50
    min_aln_len: int = 50
    min_aln_frac: float = 0.50
    min_mapq: int = 0

    flank: int = 0
    flank_f: Optional[int] = None
    flank_r: Optional[int] = None

    save_paf: object = False  # False | True (const) | str path
    minimap2: str = "minimap2"
    tmpdir: Optional[Path] = None


@dataclass
class Bait:
    label: str
    bait_id: str  # label#idx
    seq: str


@dataclass
class BestHit:
    pool_id: str
    label: str
    bait_count: int
    status: str  # PASS|NO_HIT|FILTERED
    nhits_raw: int
    nhits_pass: int

    # best-hit fields (only valid when PASS)
    best_bait_id: Optional[str] = None
    target: Optional[str] = None
    strand: Optional[str] = None
    qlen: Optional[int] = None
    qstart: Optional[int] = None
    qend: Optional[int] = None
    tlen: Optional[int] = None
    tstart: Optional[int] = None
    tend: Optional[int] = None
    aln_len: Optional[int] = None
    matches: Optional[int] = None
    identity: Optional[float] = None
    mapq: Optional[int] = None

    flank_f: Optional[int] = None
    flank_r: Optional[int] = None
    extract_start: Optional[int] = None
    extract_end: Optional[int] = None
    extract_len: Optional[int] = None


def run_fishing_pipeline(cfg: FishingConfig) -> None:
    _ensure_minimap2(cfg.minimap2)

    baits = _read_baits(cfg.baits_path)
    label_to_baits: Dict[str, List[Bait]] = {}
    for b in baits:
        label_to_baits.setdefault(b.label, []).append(b)

    pools = _resolve_pools(cfg.pool_paths, cfg.pool_list_path)

    out_tsv = cfg.out_prefix.with_suffix(".tsv")
    out_fasta = cfg.out_prefix.with_suffix(".fasta")
    out_presence = cfg.out_prefix.with_suffix(".presence.tsv")

    save_paf_path = _resolve_save_paf(cfg)

    # Prepare tmp dir
    tmp_root = str(cfg.tmpdir) if cfg.tmpdir is not None else None
    with tempfile.TemporaryDirectory(dir=tmp_root) as td:
        td_path = Path(td)

        bait_fa = td_path / "baits.fa"
        _write_baits_fasta(baits, bait_fa)

        # Output writers (pure text)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        fasta_handle = out_fasta.open("w", encoding="utf-8")

        # Write headers
        tsv_cols = [
            "pool",
            "label",
            "bait_count",
            "status",
            "nhits_raw",
            "nhits_pass",
            "best_bait_id",
            "target",
            "strand",
            "qlen",
            "qstart",
            "qend",
            "tlen",
            "tstart",
            "tend",
            "aln_len",
            "matches",
            "identity",
            "mapq",
            "flank_f",
            "flank_r",
            "extract_start",
            "extract_end",
            "extract_len",
        ]
        with out_tsv.open("w", encoding="utf-8") as tsv_out, out_presence.open("w", encoding="utf-8") as pres_out:
            tsv_out.write("\t".join(tsv_cols) + "\n")
            pres_out.write("pool\tlabel\tpresent\n")

            # If saving PAF, open once and append per pool
            paf_handle = None
            if save_paf_path is not None:
                Path(save_paf_path).parent.mkdir(parents=True, exist_ok=True)
                paf_handle = open(save_paf_path, "w", encoding="utf-8")

            try:
                for pool_id, pool_path in tqdm(pools, desc="Fishing pools", unit="pool"):
                    pool_file = Path(pool_path)
                    paf_text = _run_minimap2(
                        minimap2=cfg.minimap2,
                        preset=cfg.preset,
                        threads=_threads_value(cfg.threads),
                        ref=pool_file,
                        query=bait_fa,
                    )

                    if paf_handle is not None:
                        paf_handle.write(f"# pool_id={pool_id}\tpath={pool_path}\n")
                        paf_handle.write(paf_text)
                        if not paf_text.endswith("\n"):
                            paf_handle.write("\n")

                    records = list(parse_paf_lines(paf_text))

                    # Group hits by label
                    # qname is bait_id = label#idx
                    baitid_to_label = {b.bait_id: b.label for b in baits}

                    raw_by_label: Dict[str, List[PafRecord]] = {lab: [] for lab in label_to_baits.keys()}
                    pass_by_label: Dict[str, List[PafRecord]] = {lab: [] for lab in label_to_baits.keys()}

                    for r in records:
                        lab = baitid_to_label.get(r.qname)
                        if lab is None:
                            continue
                        raw_by_label[lab].append(r)
                        if _pass_filters(r, cfg):
                            pass_by_label[lab].append(r)

                    # Load pool sequences once (in-memory template)
                    # NOTE: This loads the full FASTA; replace with indexed reader for huge genomes if needed.
                    pool_fa = FastaInMemory.from_fasta(pool_file)

                    for lab, bait_list in label_to_baits.items():
                        nh_raw = len(raw_by_label.get(lab, []))
                        nh_pass = len(pass_by_label.get(lab, []))

                        if nh_raw == 0:
                            hit = BestHit(
                                pool_id=pool_id,
                                label=lab,
                                bait_count=len(bait_list),
                                status="NO_HIT",
                                nhits_raw=0,
                                nhits_pass=0,
                            )
                        elif nh_pass == 0:
                            hit = BestHit(
                                pool_id=pool_id,
                                label=lab,
                                bait_count=len(bait_list),
                                status="FILTERED",
                                nhits_raw=nh_raw,
                                nhits_pass=0,
                            )
                        else:
                            best = _choose_best(pass_by_label[lab])
                            flank_f, flank_r = _resolve_flanks(cfg)

                            extract_start, extract_end = _compute_extract_interval_bait_oriented(
                                best, flank_f=flank_f, flank_r=flank_r
                            )
                            seq = pool_fa.get(best.tname, extract_start, extract_end)
                            # Orient sequence to bait direction
                            if best.strand == "-":
                                seq = reverse_complement(seq)

                            # Write FASTA (PASS only)
                            header = (
                                f">{pool_id}|{lab}|bait={best.qname}|"
                                f"{best.tname}:{extract_start}-{extract_end}({best.strand})|"
                                f"id={best.identity:.3f}|aln={best.aln_len}|m={best.matches}"
                            )
                            fasta_handle.write(header + "\n")
                            fasta_handle.write(_wrap_fasta(seq) + "\n")

                            hit = BestHit(
                                pool_id=pool_id,
                                label=lab,
                                bait_count=len(bait_list),
                                status="PASS",
                                nhits_raw=nh_raw,
                                nhits_pass=nh_pass,
                                best_bait_id=best.qname,
                                target=best.tname,
                                strand=best.strand,
                                qlen=best.qlen,
                                qstart=best.qstart,
                                qend=best.qend,
                                tlen=best.tlen,
                                tstart=best.tstart,
                                tend=best.tend,
                                aln_len=best.aln_len,
                                matches=best.matches,
                                identity=round(best.identity, 6),
                                mapq=best.mapq,
                                flank_f=flank_f,
                                flank_r=flank_r,
                                extract_start=extract_start,
                                extract_end=extract_end,
                                extract_len=extract_end - extract_start,
                            )

                        # Write TSV row
                        tsv_out.write(_hit_to_tsv(hit, tsv_cols) + "\n")
                        # presence-only
                        pres_out.write(f"{pool_id}\t{lab}\t{1 if hit.status == 'PASS' else 0}\n")

            finally:
                if paf_handle is not None:
                    paf_handle.close()
                fasta_handle.close()


def _ensure_minimap2(bin_path: str) -> None:
    if shutil.which(bin_path) is None and not Path(bin_path).exists():
        raise SystemExit(
            f"ERROR: minimap2 not found: '{bin_path}'. "
            "Please install minimap2 or provide --minimap2 /path/to/minimap2."
        )


def _threads_value(threads: int) -> int:
    if threads is None or threads <= 0:
        return max(1, os.cpu_count() or 1)
    return threads


def _read_baits(path_or_dash: str) -> List[Bait]:
    import sys

    lines: List[str]
    if path_or_dash == "-":
        lines = sys.stdin.read().splitlines()
    else:
        lines = Path(path_or_dash).read_text(encoding="utf-8").splitlines()

    baits: List[Bait] = []
    label_counts: Dict[str, int] = {}

    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) != 2:
            raise SystemExit(f"ERROR: invalid baits line (expect 2 columns label<TAB>sequence): {ln}")
        label, seq = parts[0].strip(), parts[1].strip().upper()
        if not label:
            raise SystemExit(f"ERROR: empty label in baits line: {ln}")
        if not seq:
            raise SystemExit(f"ERROR: empty sequence for label '{label}'")

        label_counts[label] = label_counts.get(label, 0) + 1
        bait_id = f"{label}#{label_counts[label]}"
        baits.append(Bait(label=label, bait_id=bait_id, seq=seq))

    if not baits:
        raise SystemExit("ERROR: no baits loaded.")
    return baits


def _resolve_pools(pool_paths: Optional[List[str]], pool_list_path: Optional[str]) -> List[Tuple[str, str]]:
    if pool_paths:
        # pool_id defaults to basename without suffix
        out: List[Tuple[str, str]] = []
        for p in pool_paths:
            pp = Path(p)
            pool_id = pp.name
            out.append((pool_id, str(pp)))
        return out

    if not pool_list_path:
        raise SystemExit("ERROR: either --pool or --pool-list is required.")

    lp = Path(pool_list_path)
    base_dir = lp.parent
    rows = lp.read_text(encoding="utf-8").splitlines()

    pools: List[Tuple[str, str]] = []
    for ln in rows:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) != 2:
            raise SystemExit(f"ERROR: invalid pool-list line (expect id<TAB>path): {ln}")
        pid = parts[0].strip()
        pth = parts[1].strip()
        if not pid or not pth:
            raise SystemExit(f"ERROR: invalid pool-list line (empty id or path): {ln}")
        path = (base_dir / pth).resolve() if not Path(pth).is_absolute() else Path(pth)
        pools.append((pid, str(path)))
    if not pools:
        raise SystemExit("ERROR: no pools loaded from --pool-list.")
    return pools


def _write_baits_fasta(baits: List[Bait], out_fa: Path) -> None:
    with out_fa.open("w", encoding="utf-8") as f:
        for b in baits:
            f.write(f">{b.bait_id}\n")
            f.write(_wrap_fasta(b.seq) + "\n")


def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width)) if seq else ""


def _run_minimap2(*, minimap2: str, preset: str, threads: int, ref: Path, query: Path) -> str:
    # Use stdout capture; output is PAF by default.
    cmd = [minimap2, "-x", preset, "-t", str(threads), str(ref), str(query)]
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        raise SystemExit(f"ERROR: minimap2 not found: {minimap2}")

    if proc.returncode != 0:
        raise SystemExit(
            "ERROR: minimap2 failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDERR:\n{proc.stderr}"
        )
    # minimap2 writes PAF to stdout; messages to stderr
    return proc.stdout


def _pass_filters(r: PafRecord, cfg: FishingConfig) -> bool:
    if r.aln_len < cfg.min_aln_len:
        return False
    if r.identity < cfg.min_identity:
        return False
    if r.aln_frac < cfg.min_aln_frac:
        return False
    if r.mapq < cfg.min_mapq:
        return False
    return True


def _choose_best(records: List[PafRecord]) -> PafRecord:
    # Stable deterministic ordering:
    # 1) matches desc
    # 2) identity desc
    # 3) aln_len desc
    # 4) mapq desc
    # 5) tname, tstart, tend asc
    def key(r: PafRecord):
        return (-r.matches, -r.identity, -r.aln_len, -r.mapq, r.tname, r.tstart, r.tend)

    return sorted(records, key=key)[0]


def _resolve_flanks(cfg: FishingConfig) -> Tuple[int, int]:
    # Flanks are defined in bait orientation (forward/backward).
    if cfg.flank_f is not None or cfg.flank_r is not None:
        ff = cfg.flank_f if cfg.flank_f is not None else 0
        fr = cfg.flank_r if cfg.flank_r is not None else 0
        return int(ff), int(fr)
    return int(cfg.flank), int(cfg.flank)


def _compute_extract_interval_bait_oriented(r: PafRecord, *, flank_f: int, flank_r: int) -> Tuple[int, int]:
    # Flanks are defined along bait orientation.
    # If strand is '+':
    #   bait-forward flank extends upstream => target left side
    #   bait-reverse flank extends downstream => target right side
    # If strand is '-':
    #   bait-forward flank extends upstream in bait direction => target right side
    #   bait-reverse flank extends downstream in bait direction => target left side
    if r.strand == "+":
        left = flank_f
        right = flank_r
    else:
        left = flank_r
        right = flank_f

    extract_start = max(0, r.tstart - left)
    extract_end = min(r.tlen, r.tend + right)
    if extract_end < extract_start:
        extract_end = extract_start
    return extract_start, extract_end


def _hit_to_tsv(hit: BestHit, cols: List[str]) -> str:
    d = {
        "pool": hit.pool_id,
        "label": hit.label,
        "bait_count": hit.bait_count,
        "status": hit.status,
        "nhits_raw": hit.nhits_raw,
        "nhits_pass": hit.nhits_pass,
        "best_bait_id": hit.best_bait_id,
        "target": hit.target,
        "strand": hit.strand,
        "qlen": hit.qlen,
        "qstart": hit.qstart,
        "qend": hit.qend,
        "tlen": hit.tlen,
        "tstart": hit.tstart,
        "tend": hit.tend,
        "aln_len": hit.aln_len,
        "matches": hit.matches,
        "identity": hit.identity,
        "mapq": hit.mapq,
        "flank_f": hit.flank_f,
        "flank_r": hit.flank_r,
        "extract_start": hit.extract_start,
        "extract_end": hit.extract_end,
        "extract_len": hit.extract_len,
    }
    return "\t".join("" if d.get(c) is None else str(d.get(c)) for c in cols)


def _resolve_save_paf(cfg: FishingConfig) -> Optional[str]:
    if cfg.save_paf is False:
        return None
    if cfg.save_paf is True:
        return str(cfg.out_prefix.with_suffix(".paf"))
    # if user passed a value, argparse sets it to string
    if isinstance(cfg.save_paf, str):
        return cfg.save_paf
    # fallback
    return str(cfg.out_prefix.with_suffix(".paf"))
