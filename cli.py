"""Command-line interface for the DNA-barcode identification/phylogeny pipeline.

Installed as the `barcode-pipeline` command; `run_barcode.py` in the repository
root calls the same entry point for people working from a clone.
"""
import argparse

from .pipeline import run


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help=".ab1/FASTA file or a folder of them")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--locus", default=None,
                    help="ITS, LSU, SSU, TEF1, RPB2, BTUB, 16S, RPOB "
                         "(default: inferred from file names)")
    ap.add_argument("--email", default=None, help="contact address for NCBI services")
    ap.add_argument("--references", default=None,
                    help="FASTA of reference sequences or a text file of GenBank "
                         "accessions; for a multi-locus run tag them per locus, "
                         "e.g. ITS=its_refs.fasta,16S=16s_refs.txt")
    ap.add_argument("--outgroup", nargs="*", default=None,
                    help="tip labels to root on (default: midpoint rooting)")
    ap.add_argument("--local-db", dest="local_db", default=None,
                    help="identify against this curated FASTA instead of NCBI BLAST")
    ap.add_argument("--skip-blast", action="store_true")
    ap.add_argument("--no-merge", dest="merge", action="store_false",
                    help="keep forward and reverse reads separate")
    ap.add_argument("--no-entrez-filter", dest="entrez_filter", action="store_false",
                    help="search all of the database instead of the locus taxon filter")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--bootstrap", type=int, default=1000,
                    help="ultrafast bootstrap replicates (>=1000; 0 disables)")
    ap.add_argument("--model", default=None, help="fix the model, e.g. GTR+G (default: ModelFinder)")
    ap.add_argument("--min-occupancy", type=float, default=None,
                    help="alignment column occupancy threshold (default 0.5)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    manifest = run(a.input, outdir=a.outdir, locus=a.locus, email=a.email,
                   references=a.references, outgroup=a.outgroup,
                   skip_blast=a.skip_blast, local_reference_fasta=a.local_db,
                   merge_pairs=a.merge, threads=a.threads,
                   bootstrap=(a.bootstrap or None), model=a.model,
                   min_occupancy=a.min_occupancy, use_entrez_filter=a.entrez_filter,
                   seed=a.seed)
    print(f"done: {manifest['identifications']} identifications, "
          f"{len(manifest['trees'])} tree(s) -> {a.outdir}")


if __name__ == "__main__":
    main()
