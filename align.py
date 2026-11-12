"""Multiple sequence alignment (FAMSA) and alignment trimming.

FAMSA is used because it ships a native Windows build and scales to large
barcode sets; the pipeline writes plain FASTA, so an alignment produced by
MAFFT or MUSCLE elsewhere can be dropped in at the same step.

Trimming removes poorly occupied columns - the ragged 5'/3' ends and the
indel-rich stretches of ITS1/ITS2 where homology is not credible. The rule is
the equivalent of trimAl's gap-threshold mode: keep columns whose occupancy
(fraction of non-gap residues) is at or above `min_occupancy`, then discard
retained blocks shorter than `min_block_length`.
"""
from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .config import ALIGN_DEFAULTS


def align_records(records: list[SeqRecord], threads: int = 0,
                  guide_tree: str = "upgma") -> list[SeqRecord]:
    """Align with FAMSA. Returns aligned SeqRecords in input order."""
    from pyfamsa import Aligner, Sequence

    if len(records) < 2:
        return list(records)
    seqs = [Sequence(r.id.encode(), str(r.seq).upper().encode()) for r in records]
    msa = Aligner(threads=threads, guide_tree=guide_tree).align(seqs)
    by_id = {s.id.decode(): s.sequence.decode() for s in msa}
    return [SeqRecord(Seq(by_id[r.id]), id=r.id, description="") for r in records]


def column_occupancy(aln: list[SeqRecord]) -> list[float]:
    n = len(aln)
    width = len(aln[0].seq)
    cols = []
    for i in range(width):
        filled = sum(1 for r in aln if r.seq[i] not in "-.~?")
        cols.append(filled / n)
    return cols


def trim_alignment(aln: list[SeqRecord], min_occupancy: float | None = None,
                   min_block_length: int | None = None) -> tuple[list[SeqRecord], dict]:
    """Occupancy-based column trimming. Returns (trimmed, statistics)."""
    occ_t = ALIGN_DEFAULTS["min_occupancy"] if min_occupancy is None else min_occupancy
    blk = ALIGN_DEFAULTS["min_block_length"] if min_block_length is None else min_block_length
    if not aln:
        return [], {}
    occ = column_occupancy(aln)
    keep = [o >= occ_t for o in occ]

    # drop isolated short blocks of kept columns
    i = 0
    while i < len(keep):
        if keep[i]:
            j = i
            while j < len(keep) and keep[j]:
                j += 1
            if (j - i) < blk:
                for k in range(i, j):
                    keep[k] = False
            i = j
        else:
            i += 1

    idx = [i for i, k in enumerate(keep) if k]
    trimmed = [SeqRecord(Seq("".join(str(r.seq)[i] for i in idx)), id=r.id,
                         description="") for r in aln]
    stats = {
        "n_sequences": len(aln),
        "columns_before": len(occ),
        "columns_after": len(idx),
        "columns_removed": len(occ) - len(idx),
        "pct_columns_retained": round(100 * len(idx) / len(occ), 1) if occ else 0.0,
        "mean_occupancy_before": round(sum(occ) / len(occ), 3) if occ else 0.0,
        "mean_occupancy_after": (round(sum(occ[i] for i in idx) / len(idx), 3) if idx else 0.0),
        "min_occupancy_threshold": occ_t,
        "parsimony_informative_sites": parsimony_informative(trimmed),
        "variable_sites": variable_sites(trimmed),
    }
    return trimmed, stats


def variable_sites(aln: list[SeqRecord]) -> int:
    if not aln:
        return 0
    n = 0
    for i in range(len(aln[0].seq)):
        chars = {str(r.seq)[i].upper() for r in aln} - set("-.~?N")
        if len(chars) > 1:
            n += 1
    return n


def parsimony_informative(aln: list[SeqRecord]) -> int:
    if not aln:
        return 0
    n = 0
    for i in range(len(aln[0].seq)):
        counts: dict[str, int] = {}
        for r in aln:
            c = str(r.seq)[i].upper()
            if c in "ACGT":
                counts[c] = counts.get(c, 0) + 1
        if sum(1 for v in counts.values() if v >= 2) >= 2:
            n += 1
    return n


def write_fasta(records: list[SeqRecord], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(path), "fasta")
    return path
