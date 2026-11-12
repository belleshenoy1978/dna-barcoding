"""Fetch reference sequences from NCBI by accession.

Used to add type strains / ex-type sequences and an outgroup to a tree so the
isolates are placed against named anchors rather than against each other.
"""
from __future__ import annotations

from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_accessions(accessions: list[str], email: str, api_key: str | None = None,
                     batch_size: int = 50, cache: str | Path | None = None) -> list[SeqRecord]:
    """Download nucleotide records for a list of accessions (cached if given)."""
    accessions = [a.strip() for a in accessions if a and a.strip()]
    if not accessions:
        return []
    cache = Path(cache) if cache else None
    if cache and cache.exists():
        have = {r.id.split(".")[0]: r for r in SeqIO.parse(str(cache), "fasta")}
        if all(a.split(".")[0] in have for a in accessions):
            return [have[a.split(".")[0]] for a in accessions]

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    out: list[SeqRecord] = []
    for i in range(0, len(accessions), batch_size):
        chunk = accessions[i:i + batch_size]
        with Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="fasta",
                           retmode="text") as h:
            out.extend(SeqIO.parse(h, "fasta"))
    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(out, str(cache), "fasta")
    return out


def label_reference(rec: SeqRecord, max_words: int = 3) -> str:
    """Tree-friendly label: 'Genus_species_ACCESSION'."""
    desc = rec.description.replace(rec.id, "").strip()
    words = [w for w in desc.split() if w][:max_words]
    name = "_".join(w.strip(",;") for w in words[:2]) if words else "reference"
    acc = rec.id.split("|")[-1]
    return f"{name}_{acc}".replace("__", "_")
