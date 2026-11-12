# barcode-pipeline

DNA-barcode identification and phylogeny for fungal and bacterial isolates:
a reproducible route from Sanger chromatograms to an identification table and a
labelled maximum-likelihood tree, with a browser interface for day-to-day use.

Loci supported out of the box: **ITS, LSU, SSU, TEF1, RPB2, BTUB** (fungi) and
**16S, rpoB** (bacteria).

---

## Install

```bash
git clone https://github.com/belleshenoy1978/dna-barcoding.git
cd dna-barcoding
pip install -e ".[web]"
```

Python 3.12 or newer. All external tools are Python wheels - there is nothing
to compile and no separate BLAST, MAFFT or IQ-TREE installation to manage.

## Web interface

```bash
python -m webapp          # then open http://127.0.0.1:8000
```

Drop `.ab1` chromatograms or FASTA files onto the page, choose the options and
run. The page shows progress stage by stage (the BLAST search is the slow one),
then the identification table with confidence badges and GenBank links, the
read-QC table, the BLAST hit list, the tree figure, and a download of every
output file as a zip.

Uploaded reads and all intermediate files stay on the machine running the
server. The only outbound request is the NCBI search, and ticking **skip BLAST**
removes even that.

## Command line

```bash
# folder of .ab1 traces; locus taken from the file names
barcode-pipeline --input traces/ --outdir results --email you@institute.org

# one locus, GenBank reference accessions, outgroup rooting
barcode-pipeline --input ITS_queries.fasta --locus ITS \
    --references its_refs.txt --outgroup Saccharomyces_cerevisiae_NR_111007 \
    --email you@institute.org --threads 4

# several loci, references tagged per locus
barcode-pipeline --input traces/ \
    --references ITS=its_refs.fasta,16S=16s_refs.fasta --email you@institute.org

# offline: align and build trees from sequences you already have
barcode-pipeline --input seqs/ --locus 16S --skip-blast
```

From a clone without installing, `python run_barcode.py ...` is equivalent.

Try it on the bundled example data:

```bash
barcode-pipeline --input examples/queries --outdir demo \
  --references ITS=examples/references_ITS.fasta,16S=examples/references_16S.fasta \
  --email you@institute.org
```

## What each stage does

1. **Read handling** (`barcode_pipeline/qc.py`) - reads `.ab1` or FASTA. Quality
   trimming uses the modified Mott algorithm (the Phred/Phrap procedure): bases
   are scored `cutoff - P(error)` and the maximum-scoring segment is kept.
   Forward and reverse reads of an isolate are recognised from primer names, the
   reverse read is reverse-complemented, and the pair is merged into a consensus
   where the higher-quality base wins and unresolved conflicts become IUPAC
   ambiguity codes. Reads are scored on length, mean Q, %Q20 and ambiguity.
2. **Identification** (`identify.py`) - taxon-restricted megablast against NCBI,
   or alignment against a curated in-house FASTA. Every search is cached on
   disk, so re-runs never repeat a query.
3. **Alignment** (`align.py`) - FAMSA, then removal of columns whose occupancy
   falls below a threshold (the equivalent of trimAl's gap-threshold mode), with
   short retained blocks discarded. Reports variable and parsimony-informative
   sites.
4. **Phylogeny** (`tree.py`) - IQ-TREE 2 through the `piqtree` bindings:
   ModelFinder selects the model by BIC (AIC/AICc also reported), then an ML
   tree with 1000 ultrafast bootstrap replicates. Rooting on a named outgroup,
   or midpoint.
5. **Figure** (`plot.py`) - phylogram with query isolates in bold blue, UFBoot
   values on branches, and a substitutions-per-site scale bar.

Locus and primer direction are read from sequence and file names (`ITS1F`,
`ITS4`, `LR0R`, `28S`, `TEF`, `RPB2`, `benA`, `16S`, `27F`, `1492R`, ...), then
from `--locus`.

## Output layout

```
results/
  01_sequences/       <LOCUS>_queries.fasta, qc_report.csv, read_merging_report.csv
  02_identification/  identification_summary.csv, blast_top_hits.csv
  03_alignment/       <LOCUS>_aligned.fasta, <LOCUS>_trimmed.fasta, alignment_stats.json
  04_tree/            <LOCUS>_ML_unrooted.treefile, <LOCUS>_ML_rooted.treefile,
                      <LOCUS>_tree.png, tree_stats.json
  run_manifest.json, SUMMARY.md
```

## Identification thresholds

Set in `barcode_pipeline/config.py` and applied per locus - edit them to match
your own criteria.

| Locus | Species | Genus | Min. query coverage |
|-------|---------|-------|---------------------|
| ITS   | 98.5 %  | 95 %  | 85 % |
| LSU   | 99.5 %  | 97 %  | 85 % |
| SSU   | 99.8 %  | 98 %  | 85 % |
| TEF1 / RPB2 / BTUB | 99 % | 95 % | 80 % |
| 16S   | 98.7 %  | 94.5 %| 85 % |
| rpoB  | 97.7 %  | 94 %  | 80 % |

The assignment is deliberately conservative:

* A hit is usable only if it reaches the locus coverage threshold.
* Records flagged `uncultured` / `unidentified` / `environmental` are set aside
  when named records are available.
* Hits within 0.5 % identity of the best hit form an ambiguity band. If that
  band holds more than one binomial the call drops to genus and `ambiguous` is
  set, with the competing names listed. A record identified only to genus
  (`Aspergillus sp.`) is not treated as a competing name.
* Below the species threshold the call is genus; below the genus threshold it is
  reported as "above genus".

Two caveats no threshold can remove. ITS does not separate species in
*Fusarium*, *Penicillium*, *Aspergillus*, *Colletotrichum* or *Trichoderma* -
add TEF1, RPB2 or benA for those. A 16S identity above 98.7 % is a screening
result, not a species claim; species-level bacterial assignment needs ANI from
genomes or a type-strain database.

## Tests

```bash
pip install -e ".[web,dev]"
pytest -q
```

Covers quality trimming, read merging, name assignment at each rank, and a full
offline run through the web API. GitHub Actions runs them on Linux and Windows
for Python 3.12 and 3.13.

## Deployment

The pipeline needs Python, so **GitHub Pages cannot host the running app** -
Pages serves static files only and would show the documentation, not the
interface. Options that do run it:

* **Locally** - `python -m webapp`, which is how most labs will use it.
* **Hugging Face Space** - free, no server administration, builds from the
  `Dockerfile` as it stands. Step-by-step instructions and the Space front-matter
  are in [`deploy/huggingface/`](deploy/huggingface/DEPLOY.md); a subdomain such
  as `barcode.example.org` can be pointed at the Space with a CNAME.
* **Container** - `docker build -t barcode-pipeline . && docker run -p 8000:8000
  barcode-pipeline`; deploy that image to any container host.
* **Shared server** - `uvicorn webapp.app:app --host 0.0.0.0 --port 8000` behind
  a reverse proxy. There is no authentication, so put it behind institutional
  access control rather than on the open internet.

### The public profile

Any container built from the `Dockerfile` starts with `BARCODE_PUBLIC=1`, which
switches on a restricted profile meant for servers strangers can reach:
identification is **disabled** (NCBI's usage policy does not welcome a public
service relaying many automated BLAST searches from one address), uploads are
capped at 40 files and 25 MB, threads and bootstrap replicates are limited, and
job directories are deleted after six hours. A banner on the page explains this
to visitors.

Run `docker run -e BARCODE_PUBLIC=0 ...` for a private server that should behave
like a local install. Individual limits are set with `BARCODE_ALLOW_BLAST`,
`BARCODE_MAX_FILES`, `BARCODE_MAX_UPLOAD_MB`, `BARCODE_MAX_THREADS`,
`BARCODE_MAX_BOOTSTRAP`, `BARCODE_JOB_TTL_HOURS` and `BARCODE_NOTICE`; the table
in [`deploy/huggingface/DEPLOY.md`](deploy/huggingface/DEPLOY.md) lists the
defaults.

## Software used

FAMSA (via `pyfamsa`) for alignment; IQ-TREE 2 with ModelFinder and UFBoot (via
`piqtree`) for phylogeny; Biopython for sequence handling, NCBI BLAST and
Entrez; DendroPy for rooting; FastAPI for the interface.

MAFFT and trimAl have no Windows build, so FAMSA replaces MAFFT and the
occupancy filter replaces trimAl. Both stages read and write plain FASTA, so an
alignment produced by MAFFT or trimmed by trimAl elsewhere can be substituted
into `03_alignment/` and the tree stage will use it unchanged.

`barcode_pipeline/_iqtree_worker.py` exists because IQ-TREE's bundled OpenMP
runtime collides with the one NumPy loads on Windows; inference therefore runs
in its own process with `KMP_DUPLICATE_LIB_OK=TRUE` set before any import.

## Citing

See `CITATION.cff`. Please also cite FAMSA, IQ-TREE 2, ModelFinder, UFBoot and
Biopython, which do the underlying work.

## Licence

MIT - see `LICENSE`.
