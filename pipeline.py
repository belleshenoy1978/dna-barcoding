"""End-to-end orchestration: traces -> QC -> identification -> alignment -> tree.

Run from the command line (see run_barcode.py) or call `run()` directly.
Every stage writes its own numbered folder so intermediate products can be
inspected, reused, or replaced by files prepared elsewhere.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import align as align_mod
from . import identify as id_mod
from . import plot as plot_mod
from . import qc as qc_mod
from . import references as ref_mod
from .config import LOCI, LocusConfig

SEQ_EXT = {".ab1", ".abi", ".fasta", ".fas", ".fa", ".fna", ".seq", ".txt", ".fastq", ".fq"}


def infer_locus(name: str, default: str | None = None) -> str | None:
    """Guess the locus from a file or sequence name."""
    n = name.upper()
    patterns = [
        ("ITS", r"ITS|ITS1F|ITS4|ITS5"), ("LSU", r"LSU|28S|LR0R|LR5|D1D2"),
        ("SSU", r"\bSSU\b|18S|NS1|NS4"), ("TEF1", r"TEF|EF1|EF-1|TEF1"),
        ("RPB2", r"RPB2"), ("BTUB", r"BTUB|TUB2|BENA|BT2A|BT2B"),
        ("RPOB", r"RPOB"), ("16S", r"16S|27F|1492R|1525R|515F|806R"),
    ]
    for locus, pat in patterns:
        if re.search(pat, n):
            return locus
    return default


# --------------------------------------------------------------------------
def collect_inputs(input_path: str | Path) -> list[Path]:
    p = Path(input_path)
    if p.is_file():
        return [p]
    return sorted(f for f in p.rglob("*") if f.suffix.lower() in SEQ_EXT and f.is_file())


def stage_qc(files: list[Path], outdir: Path, default_locus: str | None,
             merge_pairs: bool = True, qc_overrides: dict | None = None) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Read, trim and (optionally) merge reads; group clean sequences by locus."""
    d = outdir / "01_sequences"; d.mkdir(parents=True, exist_ok=True)
    qc_rows, trimmed_by_key = [], defaultdict(list)

    for f in files:
        try:
            recs = qc_mod.read_any(f)
        except Exception as exc:
            qc_rows.append({"read_id": f.stem, "source_file": str(f), "passed": False,
                            "flags": f"unreadable: {exc}"})
            continue
        for rec in recs:
            locus = infer_locus(rec.id, None) or infer_locus(f.name, default_locus)
            tr, q = qc_mod.assess_read(rec, source=str(f), **(qc_overrides or {}))
            row = q.as_dict(); row["locus"] = locus or ""
            qc_rows.append(row)
            if len(tr) > 0:
                trimmed_by_key[(locus, q.isolate)].append((q.direction, tr, q))

    merge_rows: list[dict] = []
    final: dict[str, list[SeqRecord]] = defaultdict(list)
    for (locus, isolate), items in sorted(trimmed_by_key.items(), key=lambda kv: str(kv[0])):
        dirs = {d_: (r, q) for d_, r, q in items}
        rec = None
        if merge_pairs and "F" in dirs and "R" in dirs:
            rec, info = qc_mod.merge_pair(dirs["F"][0], dirs["R"][0])
            merge_rows.append({"isolate": isolate, "locus": locus or "",
                               "forward_read": dirs["F"][0].id, "reverse_read": dirs["R"][0].id,
                               **info})
        if rec is None:                     # keep the best single read
            best = max(items, key=lambda t: (t[2].mean_quality or 0) * len(t[1]))
            rec = best[1]
            rec = rec.reverse_complement(id=isolate, description="") if best[0] == "R" else rec
            rec.id, rec.description = isolate, ""
        final[locus or "unassigned"].append(rec)

    paths = {}
    for locus, recs in final.items():
        p = d / f"{locus}_queries.fasta"
        SeqIO.write(recs, str(p), "fasta")
        paths[locus] = p
    qc_df = pd.DataFrame(qc_rows)
    merge_df = pd.DataFrame(merge_rows)
    qc_df.to_csv(d / "qc_report.csv", index=False)
    if len(merge_df):
        merge_df.to_csv(d / "read_merging_report.csv", index=False)
    return {"sequences": paths, "records": final}, qc_df, merge_df


def stage_identify(records_by_locus: dict[str, list[SeqRecord]], outdir: Path,
                   email: str | None = None, local_reference_fasta: str | Path | None = None,
                   skip_blast: bool = False, hitlist_size: int = 25,
                   cache_dir: str | Path = "blast_cache",
                   use_entrez_filter: bool = True) -> pd.DataFrame:
    """BLAST (or local-reference) identification for every query sequence."""
    d = outdir / "02_identification"; d.mkdir(parents=True, exist_ok=True)
    id_rows, hit_rows = [], []
    for locus, recs in records_by_locus.items():
        cfg = LOCI.get(locus)
        if cfg is None:
            continue
        for rec in recs:
            if local_reference_fasta:
                hits = id_mod.search_local_references(rec, local_reference_fasta,
                                                      top_n=hitlist_size)
                source = f"local:{Path(local_reference_fasta).name}"
            elif skip_blast:
                continue
            else:
                hits = id_mod.blast_ncbi(rec, cfg, cache_dir=cache_dir, email=email,
                                         hitlist_size=hitlist_size,
                                         use_entrez_filter=use_entrez_filter)
                source = f"NCBI:{cfg.blast_database}"
            for rank, h in enumerate(hits[:10], 1):
                hit_rows.append({"isolate": rec.id, "locus": locus, "hit_rank": rank,
                                 "search": source, **h})
            call = id_mod.assign_identity(hits, cfg)
            id_rows.append({"isolate": rec.id, "query_length": len(rec),
                            "search": source, **call})
    id_df = pd.DataFrame(id_rows)
    hit_df = pd.DataFrame(hit_rows)
    if len(id_df):
        id_df.to_csv(d / "identification_summary.csv", index=False)
    if len(hit_df):
        hit_df.to_csv(d / "blast_top_hits.csv", index=False)
    return id_df


def stage_align(records_by_locus: dict[str, list[SeqRecord]], outdir: Path,
                reference_records: dict[str, list[SeqRecord]] | None = None,
                threads: int = 0, min_occupancy: float | None = None) -> dict:
    """Align queries (plus references) per locus and trim poorly occupied columns."""
    d = outdir / "03_alignment"; d.mkdir(parents=True, exist_ok=True)
    out = {}
    for locus, recs in records_by_locus.items():
        pool = list(recs) + list((reference_records or {}).get(locus, []))
        if len(pool) < 3:
            out[locus] = {"skipped": f"only {len(pool)} sequences"}
            continue
        aligned = align_mod.align_records(pool, threads=threads)
        align_mod.write_fasta(aligned, d / f"{locus}_aligned.fasta")
        trimmed, stats = align_mod.trim_alignment(aligned, min_occupancy=min_occupancy)
        align_mod.write_fasta(trimmed, d / f"{locus}_trimmed.fasta")
        stats.update({"aligned_fasta": str(d / f"{locus}_aligned.fasta"),
                      "trimmed_fasta": str(d / f"{locus}_trimmed.fasta"),
                      "query_labels": [r.id for r in recs]})
        out[locus] = stats
    (d / "alignment_stats.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


def stage_tree(align_results: dict, outdir: Path, model: str | None = None,
               bootstrap: int | None = None, threads: int = 1,
               outgroup: list[str] | None = None, seed: int = 1) -> dict:
    """ModelFinder + ML tree + UFBoot per locus, then a labelled figure."""
    from . import tree as tree_mod
    d = outdir / "04_tree"; d.mkdir(parents=True, exist_ok=True)
    results = {}
    for locus, stats in align_results.items():
        if "trimmed_fasta" not in stats or stats.get("n_sequences", 0) < 4:
            results[locus] = {"skipped": stats.get("skipped", "fewer than 4 sequences")}
            continue
        res = tree_mod.infer_tree(stats["trimmed_fasta"], model=model,
                                  bootstrap_replicates=bootstrap, threads=threads,
                                  seed=seed, model_set=LOCI[locus].substitution_models
                                  if locus in LOCI else None)
        rooted = (tree_mod.root_with_outgroup(res["newick"], outgroup) if outgroup
                  else tree_mod.midpoint_root_newick(res["newick"]))
        tree_mod.write_newick(res["newick"], d / f"{locus}_ML_unrooted.treefile")
        tree_mod.write_newick(rooted, d / f"{locus}_ML_rooted.treefile")
        fig = plot_mod.plot_tree(
            rooted, d / f"{locus}_tree.png",
            query_labels=set(stats.get("query_labels", [])),
            title=f"{locus} maximum-likelihood phylogeny  |  {res['model']}"
                  + (f"  |  UFBoot {res['bootstrap_replicates']}" if res["bootstrap_replicates"] else ""),
            support_label="UFBoot" if res["bootstrap_replicates"] else "support",
            support_threshold=70.0)
        res.update({"rooted_newick": str(d / f"{locus}_ML_rooted.treefile"),
                    "figure": str(fig), "rooting": "outgroup" if outgroup else "midpoint"})
        res.pop("newick", None)
        results[locus] = res
    (d / "tree_stats.json").write_text(json.dumps(results, indent=2, default=str),
                                       encoding="utf-8")
    return results


# --------------------------------------------------------------------------
def run(input_path: str | Path, outdir: str | Path = "results",
        locus: str | None = None, email: str | None = None,
        references: str | Path | None = None, outgroup: list[str] | None = None,
        skip_blast: bool = False, local_reference_fasta: str | Path | None = None,
        merge_pairs: bool = True, threads: int = 1, bootstrap: int | None = None,
        model: str | None = None, min_occupancy: float | None = None,
        use_entrez_filter: bool = True, seed: int = 1) -> dict:
    """Run the whole pipeline and return a manifest of what was produced."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    files = collect_inputs(input_path)
    if not files:
        raise FileNotFoundError(f"no sequence files under {input_path}")

    seqs, qc_df, merge_df = stage_qc(files, outdir, locus, merge_pairs=merge_pairs)
    records = seqs["records"]

    ref_records: dict[str, list[SeqRecord]] = {}
    if references:
        ref_records = _references_by_locus(references, list(records), email=email,
                                           outdir=outdir)

    id_df = stage_identify(records, outdir, email=email, skip_blast=skip_blast,
                           local_reference_fasta=local_reference_fasta,
                           use_entrez_filter=use_entrez_filter)
    aln = stage_align(records, outdir, reference_records=ref_records, threads=threads,
                      min_occupancy=min_occupancy)
    trees = stage_tree(aln, outdir, model=model, bootstrap=bootstrap, threads=threads,
                       outgroup=outgroup, seed=seed)

    manifest = {
        "input": str(input_path),
        "files_read": len(files),
        "reads_assessed": int(len(qc_df)),
        "reads_passing_qc": int(qc_df["passed"].sum()) if "passed" in qc_df else 0,
        "pairs_merged": int(merge_df["merged"].sum()) if len(merge_df) and "merged" in merge_df else 0,
        "loci": {k: len(v) for k, v in records.items()},
        "identifications": int(len(id_df)),
        "alignments": {k: {kk: vv for kk, vv in v.items() if kk != "query_labels"}
                       for k, v in aln.items()},
        "trees": trees,
    }
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str),
                                              encoding="utf-8")
    _write_summary(outdir, qc_df, id_df, aln, trees)
    return manifest


def _references_by_locus(references: str, loci: list[str], email: str | None,
                         outdir: Path) -> dict[str, list[SeqRecord]]:
    """Resolve the --references argument into reference sequences per locus.

    Accepts either a single path (FASTA or a list of GenBank accessions) or a
    comma-separated, locus-tagged list such as
    ``ITS=its_refs.fasta,16S=16s_refs.txt``. An untagged path is used for every
    locus only when the run contains a single locus; with several loci it is
    assigned by the locus named in each reference's own description, because
    mixing markers in one alignment is never correct.
    """
    out: dict[str, list[SeqRecord]] = {}
    tagged = [chunk for chunk in str(references).split(",") if chunk.strip()]
    explicit = all("=" in c for c in tagged) and len(tagged) >= 1

    if explicit:
        for chunk in tagged:
            loc, path = chunk.split("=", 1)
            out[loc.strip().upper()] = _load_references(path.strip(), email, outdir)
        return out

    refs = _load_references(str(references).strip(), email, outdir)
    real_loci = [l for l in loci if l != "unassigned"]
    if len(real_loci) <= 1:
        return {(real_loci[0] if real_loci else "unassigned"): refs}

    for r in refs:
        loc = infer_locus(f"{r.id} {r.description}", None)
        if loc:
            out.setdefault(loc, []).append(r)
    return out


def _load_references(references: str | Path, email: str | None, outdir: Path) -> list[SeqRecord]:
    p = Path(references)
    if p.suffix.lower() in {".fasta", ".fas", ".fa", ".fna"}:
        return list(SeqIO.parse(str(p), "fasta"))
    accs = [ln.strip() for ln in p.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.startswith("#")]
    if not email:
        raise ValueError("an email address is required for NCBI Entrez downloads")
    recs = ref_mod.fetch_accessions(accs, email=email,
                                    cache=outdir / "01_sequences" / "references.fasta")
    for r in recs:
        r.id = ref_mod.label_reference(r)
        r.description = ""
    return recs


def _markdown_table(df) -> str:
    """Minimal markdown table (avoids a hard dependency on tabulate)."""
    cols = list(df.columns)
    rows = ["| " + " | ".join(str(c) for c in cols) + " |",
            "| " + " | ".join("---" for _ in cols) + " |"]
    for _, r in df.iterrows():
        rows.append("| " + " | ".join("" if v is None else str(v) for v in r.tolist()) + " |")
    return "\n".join(rows)


def _write_summary(outdir: Path, qc_df, id_df, aln, trees) -> None:
    lines = ["# Barcode run summary", ""]
    if len(qc_df):
        ok = int(qc_df["passed"].sum()) if "passed" in qc_df else 0
        lines += [f"- reads assessed: **{len(qc_df)}**, passing QC: **{ok}**", ""]
    if len(id_df):
        cols = ["isolate", "locus", "identification", "rank", "confidence",
                "top_identity", "top_coverage", "top_accession"]
        lines += ["## Identifications", "", _markdown_table(id_df[cols]), ""]
    for locus, s in aln.items():
        if "columns_after" in s:
            lines.append(f"- **{locus}** alignment: {s['n_sequences']} seqs, "
                         f"{s['columns_before']} -> {s['columns_after']} columns "
                         f"({s['pct_columns_retained']}% retained), "
                         f"{s['parsimony_informative_sites']} parsimony-informative sites")
    lines.append("")
    for locus, t in trees.items():
        if "model" in t:
            lines.append(f"- **{locus}** ML tree: model {t['model']}, "
                         f"{t['n_tips']} tips, support {t['support_type']}")
    (outdir / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
