"""Barcode identification: NCBI BLAST search, hit parsing and name assignment.

Two search back-ends:
  * `blast_ncbi`  - remote megablast/blastn against an NCBI database (nt by
    default), optionally restricted with an Entrez query. Results are cached on
    disk per sequence so a re-run never repeats a search.
  * `search_local_references` - alignment against a curated in-house FASTA
    (type strains, your own vouchers) with Biopython's PairwiseAligner. Use it
    when you want identification against a trusted reference set rather than
    the whole of GenBank.

Name assignment is threshold-based and deliberately cautious: it reports the
rank the data support (species / genus / family-or-higher) and flags cases
where the top hits disagree.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

from .config import LocusConfig

UNCULTURED = re.compile(r"uncultured|unidentified|environmental|metagenome|"
                        r"\bfungal sp\b|\bbacterium\b$", re.I)
_GREEK = {"alpha", "beta", "gamma", "delta", "var", "f", "cf", "aff"}


# --------------------------------------------------------------------------
# remote BLAST
# --------------------------------------------------------------------------
def blast_ncbi(rec: SeqRecord, cfg: LocusConfig, cache_dir: str | Path = "blast_cache",
               hitlist_size: int = 25, expect: float = 1e-10,
               email: str | None = None, use_entrez_filter: bool = True,
               force: bool = False, retries: int = 2) -> list[dict]:
    """Run one remote BLAST search and return parsed hits (cached on disk)."""
    from Bio import Blast

    if email:
        Blast.email = email
    cache_dir = Path(cache_dir); cache_dir.mkdir(parents=True, exist_ok=True)
    eq = cfg.entrez_query if use_entrez_filter else None
    key = hashlib.sha1(
        f"{str(rec.seq).upper()}|{cfg.blast_database}|{cfg.megablast}|{eq}|{hitlist_size}"
        .encode()).hexdigest()[:16]
    cache = cache_dir / f"{rec.id}_{key}.json"
    if cache.exists() and not force:
        return json.loads(cache.read_text(encoding="utf-8"))["hits"]

    last = None
    for attempt in range(retries + 1):
        try:
            stream = Blast.qblast(
                "blastn", cfg.blast_database, f">{rec.id}\n{rec.seq}",
                megablast=cfg.megablast, hitlist_size=hitlist_size, expect=expect,
                entrez_query=eq or "(none)", format_type="XML",
            )
            result = Blast.read(stream)
            break
        except Exception as exc:                      # network / NCBI throttling
            last = exc
            if attempt == retries:
                raise RuntimeError(f"BLAST failed for {rec.id}: {exc}") from exc
            time.sleep(20 * (attempt + 1))

    hits = _parse_blast(result, len(rec))
    cache.write_text(json.dumps({"query": rec.id, "query_length": len(rec),
                                 "database": cfg.blast_database,
                                 "entrez_query": eq, "hits": hits}, indent=1),
                     encoding="utf-8")
    return hits


def _parse_blast(result, query_length: int) -> list[dict]:
    hits = []
    for hit in result:
        best = None
        for hsp in hit:
            ann = dict(getattr(hsp, "annotations", {}) or {})
            length = int(getattr(hsp, "length", 0) or len(hsp[0]))
            ident = int(ann.get("identity", 0) or 0)
            coords = getattr(hsp, "coordinates", None)
            if coords is not None and len(coords):
                qs, qe = int(coords[0].min()), int(coords[0].max())
            else:
                qs, qe = 0, length
            row = {
                "pct_identity": round(100 * ident / length, 2) if length else 0.0,
                "alignment_length": length,
                "identities": ident,
                "gaps": int(ann.get("gaps", 0) or 0),
                "evalue": float(ann.get("evalue", 1.0)),
                "bit_score": float(ann.get("bit score", ann.get("bit_score", 0)) or 0),
                "query_start": qs + 1, "query_end": qe,
                "query_coverage": round(100 * (qe - qs) / query_length, 1) if query_length else 0.0,
            }
            if best is None or row["bit_score"] > best["bit_score"]:
                best = row
        if best is None:
            continue
        target = getattr(hit, "target", None)
        acc = getattr(hit, "accession", None) or getattr(target, "id", "") or ""
        acc = acc.split("|")[-2] if acc.count("|") >= 2 else acc
        desc = (getattr(target, "description", "") or "").strip()
        best.update({"accession": acc, "title": desc,
                     "organism": organism_from_title(desc),
                     "uncultured": bool(UNCULTURED.search(desc))})
        hits.append(best)
    hits.sort(key=lambda h: (-h["bit_score"], h["evalue"]))
    return hits


def organism_from_title(title: str) -> str:
    """Best-effort binomial from a GenBank definition line."""
    t = re.sub(r"^\S+\s+", "", title) if title.startswith(("gi|", "ref|")) else title
    t = t.strip()
    if not t:
        return ""
    words = t.split()
    if not words:
        return ""
    genus = words[0].strip(",;")
    if len(words) > 1:
        sp = words[1].strip(",;")
        if re.fullmatch(r"[a-z-]+", sp) and sp not in _GREEK:
            return f"{genus} {sp}"
        if sp in {"sp.", "cf.", "aff.", "spp."}:
            return f"{genus} {sp}"
    return genus


def genus_of(organism: str) -> str:
    return organism.split()[0] if organism else ""


def is_named_species(organism: str) -> bool:
    """True only for a proper binomial - 'Aspergillus niger', not 'Aspergillus sp.'."""
    parts = organism.split()
    if len(parts) < 2:
        return False
    return bool(re.fullmatch(r"[a-z][a-z-]+", parts[1]))


# --------------------------------------------------------------------------
# local reference search
# --------------------------------------------------------------------------
def search_local_references(rec: SeqRecord, reference_fasta: str | Path,
                            top_n: int = 10) -> list[dict]:
    """Align a query to every sequence in a curated FASTA; return ranked hits."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score, aligner.mismatch_score = 2, -3
    aligner.open_gap_score, aligner.extend_gap_score = -5, -2

    out = []
    for ref in SeqIO.parse(str(reference_fasta), "fasta"):
        aln = aligner.align(str(rec.seq).upper(), str(ref.seq).upper())
        if len(aln) == 0:
            continue
        best = aln[0]
        a, b = best[0], best[1]
        length = len(a)
        ident = sum(1 for x, y in zip(a, b) if x == y and x != "-")
        qs, qe = int(best.coordinates[0].min()), int(best.coordinates[0].max())
        out.append({
            "accession": ref.id, "title": ref.description,
            "organism": organism_from_title(ref.description.replace(ref.id, "").strip()),
            "pct_identity": round(100 * ident / length, 2) if length else 0.0,
            "alignment_length": length, "identities": ident,
            "gaps": sum(1 for x, y in zip(a, b) if "-" in (x, y)),
            "evalue": float("nan"), "bit_score": float(best.score),
            "query_start": qs + 1, "query_end": qe,
            "query_coverage": round(100 * (qe - qs) / len(rec), 1),
            "uncultured": False,
        })
    out.sort(key=lambda h: (-h["pct_identity"], -h["bit_score"]))
    return out[:top_n]


# --------------------------------------------------------------------------
# assignment
# --------------------------------------------------------------------------
def assign_identity(hits: list[dict], cfg: LocusConfig,
                    exclude_uncultured: bool = True,
                    ambiguity_window: float = 0.5) -> dict:
    """Turn a ranked hit list into one identification with an explicit rank.

    `ambiguity_window` is the percent-identity band below the top hit inside
    which competing names are considered equally good; if those names disagree,
    the call is flagged ambiguous and reported at the shared rank.
    """
    usable = [h for h in hits if h["query_coverage"] >= cfg.min_query_coverage]
    if exclude_uncultured:
        named = [h for h in usable if not h["uncultured"] and h["organism"]]
        usable = named or usable
    if not usable:
        return {"identification": "no acceptable hit", "rank": "none",
                "confidence": "none", "top_identity": None, "top_coverage": None,
                "top_accession": None, "top_organism": None, "ambiguous": False,
                "competing_names": "", "locus": cfg.name,
                "note": f"no hit reached {cfg.min_query_coverage}% query coverage"}

    top = max(usable, key=lambda h: (h["pct_identity"], h["bit_score"]))
    band = [h for h in usable if h["pct_identity"] >= top["pct_identity"] - ambiguity_window]
    # a hit identified only to genus ("Aspergillus sp.") is not a competing name,
    # so prefer a fully named record of equal identity as the representative hit
    named_band = [h for h in band if is_named_species(h["organism"])]
    if named_band and not is_named_species(top["organism"]):
        top = max(named_band, key=lambda h: (h["pct_identity"], h["bit_score"]))
    species_set = {h["organism"] for h in band if is_named_species(h["organism"])}
    genus_set = {genus_of(h["organism"]) for h in band if h["organism"]}

    pid = top["pct_identity"]
    if pid >= cfg.species_identity and species_set:
        rank, name = "species", top["organism"]
        conf = "high" if len(species_set) == 1 else "moderate"
        if len(species_set) > 1:
            rank, name = "genus", f"{genus_of(top['organism'])} sp."
            conf = "moderate" if len(genus_set) == 1 else "low"
    elif pid >= cfg.species_identity:
        rank, name, conf = "genus", f"{genus_of(top['organism'])} sp.", "moderate"
    elif pid >= cfg.genus_identity:
        rank, name = "genus", f"{genus_of(top['organism'])} sp."
        conf = "moderate" if len(genus_set) == 1 else "low"
    else:
        rank, name, conf = "above genus", f"related to {genus_of(top['organism'])}", "low"

    return {
        "locus": cfg.name,
        "identification": name,
        "rank": rank,
        "confidence": conf,
        "top_identity": pid,
        "top_coverage": top["query_coverage"],
        "top_accession": top["accession"],
        "top_organism": top["organism"],
        "top_evalue": top["evalue"],
        "ambiguous": len(species_set) > 1,
        "competing_names": "; ".join(sorted(species_set)[:5]),
        "n_hits_in_band": len(band),
        "note": (f"thresholds: species >= {cfg.species_identity}%, "
                 f"genus >= {cfg.genus_identity}%, coverage >= {cfg.min_query_coverage}%"),
    }
