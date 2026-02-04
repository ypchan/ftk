from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Tuple
import gzip


@dataclass
class FastaInMemory:
    """Simple in-memory FASTA store for extracting subsequences.

    Notes:
    - This loads the entire FASTA into memory.
    - For very large genomes, consider replacing this with an indexed reader (pysam/pyfaidx).
    """

    seqs: Dict[str, str]

    @classmethod
    def from_fasta(cls, path: Path) -> "FastaInMemory":
        seqs: Dict[str, str] = {}
        name = None
        chunks = []

        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks).upper()
                    name = line[1:].split()[0] if line[1:] else ""
                    chunks = []
                else:
                    chunks.append(line)
            if name is not None:
                seqs[name] = "".join(chunks).upper()

        return cls(seqs=seqs)

    def get(self, name: str, start: int, end: int) -> str:
        seq = self.seqs[name]
        start = max(0, start)
        end = min(len(seq), end)
        if end < start:
            end = start
        return seq[start:end]

    def length(self, name: str) -> int:
        return len(self.seqs[name])


def reverse_complement(seq: str) -> str:
    comp = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "R": "Y",
        "Y": "R",
        "K": "M",
        "M": "K",
        "S": "S",
        "W": "W",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
    }
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))
