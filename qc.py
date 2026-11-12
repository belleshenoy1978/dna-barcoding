"""Sanger chromatogram reading, quality trimming and read pairing.

Handles the two inputs a culture-collection workflow produces:
  * raw .ab1 chromatograms (one per primer direction), and
  * already-basecalled FASTA.

Trimming uses the modified Mott algorithm, the same procedure Phred/Phrap and
most assembly packages apply: bases are scored (cutoff - error_probability)
and the maximum-scoring contiguous segment is retained.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .config import QC_DEFAULTS

AMBIG = set("NRYSWKMBDHVn")


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def read_trace(path: str | Path) -> SeqRecord:
    """Read one .ab1 chromatogram into a SeqRecord carrying phred_quality."""
    path = Path(path)
    rec = SeqIO.read(str(path), "abi")
    rec.id = path.stem
    rec.name = path.stem
    rec.description = ""
    return rec


def read_any(path: str | Path) -> list[SeqRecord]:
    """Read .ab1, .fasta/.fas/.fa/.fna or .fastq into SeqRecords."""
    path = Path(path)
    suf = path.suffix.lower()
    if suf in {".ab1", ".abi"}:
        return [read_trace(path)]
    fmt = {".fastq": "fastq", ".fq": "fastq"}.get(suf, "fasta")
    recs = list(SeqIO.parse(str(path), fmt))
    for r in recs:
        r.description = r.description or ""
    return recs


def direction_of(name: str) -> str:
    """Guess primer direction from a file/sequence name.

    Recognises the usual fungal/bacterial primer names and generic F/R tags.
    Returns 'F', 'R' or '?'.
    """
    n = name.upper()
    fwd = ["ITS1F", "ITS5", "ITS1", "ITS3", "LR0R", "NS1", "EF1-728", "EF1-983",
           "FRPB2-5F", "BT2A", "27F", "8F", "515F", "_F", "-F", ".F", "_FWD", "_FOR"]
    rev = ["ITS4", "ITS2", "LR5", "LR3", "NS4", "EF-2", "EF1-1567", "FRPB2-7CR",
           "BT2B", "1492R", "1525R", "806R", "_R", "-R", ".R", "_REV"]
    for tag in rev:                      # test reverse first: ITS4 contains "ITS"
        if n.endswith(tag) or tag in n:
            return "R"
    for tag in fwd:
        if n.endswith(tag) or tag in n:
            return "F"
    return "?"


def isolate_of(name: str) -> str:
    """Strip a primer/direction suffix to recover the isolate identifier."""
    stem = re.split(r"[_\-.]", name)
    if len(stem) > 1 and direction_of(stem[-1]) != "?":
        return "_".join(stem[:-1])
    return name


# --------------------------------------------------------------------------
# trimming
# --------------------------------------------------------------------------
def mott_trim(rec: SeqRecord, cutoff: float = 0.05,
              min_length: int = 50) -> tuple[SeqRecord, int, int]:
    """Modified-Mott quality trim. Returns (trimmed_record, start, end).

    `cutoff` is an error-probability threshold: 0.05 ~ Q13, 0.01 ~ Q20.
    If no segment scores positively, or the best segment is shorter than
    `min_length`, an empty record is returned (start == end).
    """
    qual = rec.letter_annotations.get("phred_quality")
    if not qual:                       # FASTA input: nothing to trim on
        return rec, 0, len(rec)

    scores = [cutoff - 10 ** (-q / 10.0) for q in qual]
    best = cur = 0.0
    best_s = best_e = cur_s = 0
    for i, s in enumerate(scores):
        cur += s
        if cur < 0:
            cur, cur_s = 0.0, i + 1
        elif cur > best:
            best, best_s, best_e = cur, cur_s, i + 1
    if best <= 0 or (best_e - best_s) < min_length:
        return rec[0:0], 0, 0
    return rec[best_s:best_e], best_s, best_e


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
@dataclass
class ReadQC:
    read_id: str
    isolate: str
    direction: str
    source_file: str
    raw_length: int
    trimmed_length: int
    trim_start: int
    trim_end: int
    mean_quality: float | None
    pct_q20: float | None
    ambiguities: int
    ambiguity_fraction: float
    passed: bool
    flags: str

    def as_dict(self) -> dict:
        return asdict(self)


def assess_read(rec: SeqRecord, source: str = "", **kw) -> tuple[SeqRecord, ReadQC]:
    """Trim one read and score it against the QC thresholds."""
    opt = {**QC_DEFAULTS, **kw}
    trimmed, s, e = mott_trim(rec, opt["mott_cutoff"], min_length=50)
    q = trimmed.letter_annotations.get("phred_quality") or []
    mean_q = round(sum(q) / len(q), 2) if q else None
    pct20 = round(100 * sum(1 for x in q if x >= 20) / len(q), 1) if q else None
    seq = str(trimmed.seq)
    amb = sum(1 for c in seq if c in AMBIG)
    frac = round(amb / len(seq), 4) if seq else 1.0

    flags = []
    if len(trimmed) < opt["min_trimmed_length"]:
        flags.append(f"short(<{opt['min_trimmed_length']}bp)")
    if mean_q is not None and mean_q < opt["min_mean_quality"]:
        flags.append(f"low_mean_Q({mean_q})")
    if frac > opt["max_ambiguity_fraction"]:
        flags.append(f"ambiguous({frac:.1%})")

    qc = ReadQC(
        read_id=rec.id, isolate=isolate_of(rec.id), direction=direction_of(rec.id),
        source_file=source, raw_length=len(rec), trimmed_length=len(trimmed),
        trim_start=s, trim_end=e, mean_quality=mean_q, pct_q20=pct20,
        ambiguities=amb, ambiguity_fraction=frac,
        passed=not flags, flags=";".join(flags),
    )
    trimmed.id = rec.id
    trimmed.description = ""
    return trimmed, qc


# --------------------------------------------------------------------------
# pairing forward and reverse reads
# --------------------------------------------------------------------------
def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 5
    a.mismatch_score = -4
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    return a


def merge_pair(fwd: SeqRecord, rev: SeqRecord, min_overlap: int = 40,
               min_overlap_identity: float = 0.90) -> tuple[SeqRecord | None, dict]:
    """Reverse-complement `rev`, overlap it with `fwd` and build a consensus.

    In the overlap the higher-quality base is kept; disagreements between
    equally good bases become the IUPAC ambiguity code. Returns (record, info);
    record is None when the overlap is too short or too divergent, which is the
    signal to keep the two directions as separate sequences.
    """
    rc = rev.reverse_complement()
    rc.id, rc.description = rev.id + "_rc", ""
    if "phred_quality" in rev.letter_annotations:
        rc.letter_annotations["phred_quality"] = rev.letter_annotations["phred_quality"][::-1]

    aln = _aligner().align(str(fwd.seq), str(rc.seq))
    if len(aln) == 0:
        return None, {"merged": False, "reason": "no_alignment"}
    best = aln[0]
    fa, ra = best[0], best[1]
    matches = sum(1 for x, y in zip(fa, ra) if x == y and x != "-")
    overlap = len(fa)
    ident = matches / overlap if overlap else 0.0
    info = {"overlap_length": overlap, "overlap_identity": round(ident, 4)}
    if overlap < min_overlap or ident < min_overlap_identity:
        return None, {**info, "merged": False, "reason": "weak_overlap"}

    fs, fe = best.aligned[0][0][0], best.aligned[0][-1][-1]
    rs, re_ = best.aligned[1][0][0], best.aligned[1][-1][-1]
    fq = fwd.letter_annotations.get("phred_quality") or [30] * len(fwd)
    rq = rc.letter_annotations.get("phred_quality") or [30] * len(rc)

    cons, consq = [], []
    cons += list(str(fwd.seq)[:fs]); consq += list(fq[:fs])
    fi, ri = fs, rs
    for x, y in zip(fa, ra):
        if x == "-":
            cons.append(y); consq.append(rq[ri]); ri += 1
        elif y == "-":
            cons.append(x); consq.append(fq[fi]); fi += 1
        else:
            if x == y:
                cons.append(x); consq.append(min(60, fq[fi] + rq[ri]))
            elif fq[fi] >= rq[ri] + 10:
                cons.append(x); consq.append(fq[fi] - rq[ri])
            elif rq[ri] >= fq[fi] + 10:
                cons.append(y); consq.append(rq[ri] - fq[fi])
            else:
                cons.append(_iupac(x, y)); consq.append(min(fq[fi], rq[ri]))
            fi += 1; ri += 1
    cons += list(str(rc.seq)[re_:]); consq += list(rq[re_:])

    rec = SeqRecord(Seq("".join(cons)), id=isolate_of(fwd.id), description="")
    rec.letter_annotations["phred_quality"] = [int(c) for c in consq]
    info.update({"merged": True, "consensus_length": len(rec),
                 "mismatches_in_overlap": overlap - matches})
    return rec, info


_IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
          frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}


def _iupac(a: str, b: str) -> str:
    return _IUPAC.get(frozenset((a.upper(), b.upper())), "N")
