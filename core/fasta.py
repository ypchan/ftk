import os
import gzip
from typing import Iterator, Tuple, Optional, List

from tqdm import tqdm


class _TqdmReader:
    def __init__(self, raw, pbar: tqdm):
        self.raw = raw
        self.pbar = pbar

    def read(self, n: int = -1) -> bytes:
        b = self.raw.read(n)
        if b:
            self.pbar.update(len(b))
        return b

    def readline(self, n: int = -1) -> bytes:
        b = self.raw.readline(n)
        if b:
            self.pbar.update(len(b))
        return b

    def close(self):
        try:
            self.raw.close()
        finally:
            self.pbar.close()


def iter_fasta_records(path: str, desc: str = "Reading") -> Iterator[Tuple[str, str]]:
    total = os.path.getsize(path)
    pbar = tqdm(total=total, unit="B", unit_scale=True, unit_divisor=1024, desc=desc, leave=False, mininterval=0.2)

    raw = open(path, "rb")
    wrapped = _TqdmReader(raw, pbar)

    if path.endswith(".gz"):
        handle = gzip.GzipFile(fileobj=wrapped, mode="rb")
        lines = (line.decode("utf-8", errors="replace") for line in handle)
    else:
        handle = wrapped
        lines = (line.decode("utf-8", errors="replace") for line in wrapped)

    name: Optional[str] = None
    chunks: List[str] = []

    try:
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks).upper()
    finally:
        handle.close()
