from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterator, Optional
import gzip
import sys


@dataclass(frozen=True)
class FastaRecord:
    name: str
    seq: str


def _open_text_auto(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def iter_fasta_from_handle(handle: IO[str]) -> Iterator[FastaRecord]:
    name: Optional[str] = None
    chunks: list[str] = []

    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield FastaRecord(name=name, seq="".join(chunks).upper())
            name = line[1:].split()[0] if line[1:] else ""
            chunks = []
        else:
            chunks.append(line)

    if name is not None:
        yield FastaRecord(name=name, seq="".join(chunks).upper())


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    with _open_text_auto(path) as f:
        yield from iter_fasta_from_handle(f)


def iter_fasta_stdin() -> Iterator[FastaRecord]:
    yield from iter_fasta_from_handle(sys.stdin)
