from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional


@dataclass(frozen=True)
class PafRecord:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    aln_len: int
    mapq: int
    tags: str = ""

    @property
    def identity(self) -> float:
        if self.aln_len <= 0:
            return 0.0
        return self.matches / self.aln_len

    @property
    def aln_frac(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return self.aln_len / self.qlen


def parse_paf_lines(text: str) -> Iterator[PafRecord]:
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        rec = parse_paf_line(line)
        if rec is not None:
            yield rec


def parse_paf_line(line: str) -> Optional[PafRecord]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 12:
        return None

    qname = parts[0]
    qlen = int(parts[1])
    qstart = int(parts[2])
    qend = int(parts[3])
    strand = parts[4]
    tname = parts[5]
    tlen = int(parts[6])
    tstart = int(parts[7])
    tend = int(parts[8])
    matches = int(parts[9])
    aln_len = int(parts[10])
    mapq = int(parts[11])
    tags = "\t".join(parts[12:]) if len(parts) > 12 else ""
    return PafRecord(
        qname=qname,
        qlen=qlen,
        qstart=qstart,
        qend=qend,
        strand=strand,
        tname=tname,
        tlen=tlen,
        tstart=tstart,
        tend=tend,
        matches=matches,
        aln_len=aln_len,
        mapq=mapq,
        tags=tags,
    )
