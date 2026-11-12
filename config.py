"""Locus presets: identification thresholds and search settings.

Percent-identity thresholds are the conventional cut-offs used for barcode
identification. They are deliberately conservative: a hit below the species
threshold is reported at genus level rather than forced to a species name.

References for the defaults:
  16S  species 98.7 % (Kim et al. 2014), genus 94.5 % (Yarza et al. 2014)
  ITS  species 98.5 %, genus 95 % (Nilsson et al. 2008; UNITE SH practice)
  LSU / TEF1 / RPB2 / BTUB: locus-specific practice, see notes per entry.

Edit the values here to match your laboratory's own criteria - every
downstream call reads them from this table.
"""
from dataclasses import dataclass, field


@dataclass
class LocusConfig:
    name: str
    group: str                 # "fungi" | "bacteria"
    description: str
    species_identity: float    # % identity at/above which a species name is proposed
    genus_identity: float      # % identity at/above which a genus name is proposed
    min_query_coverage: float  # % of the query aligned to the subject
    blast_database: str = "nt"
    megablast: bool = True
    entrez_query: str | None = None   # e.g. '"Fungi"[Organism]' to restrict the search
    substitution_models: tuple[str, ...] = ("GTR", "HKY", "TIM", "TN", "K80", "JC")
    notes: str = ""


LOCI: dict[str, LocusConfig] = {
    "ITS": LocusConfig(
        name="ITS", group="fungi",
        description="Internal transcribed spacer (ITS1-5.8S-ITS2), fungal primary barcode",
        species_identity=98.5, genus_identity=95.0, min_query_coverage=85.0,
        entrez_query='"Fungi"[Organism]',
        notes="Primary fungal barcode (Schoch et al. 2012). Species cut-off varies by "
              "genus; in species-rich genera (Fusarium, Penicillium, Aspergillus, "
              "Colletotrichum) ITS alone is not sufficient - add TEF1/RPB2/BTUB.",
    ),
    "LSU": LocusConfig(
        name="LSU", group="fungi",
        description="Nuclear large subunit rDNA (28S, D1/D2)",
        species_identity=99.5, genus_identity=97.0, min_query_coverage=85.0,
        entrez_query='"Fungi"[Organism]',
        notes="Conserved; useful for family/order placement, rarely for species.",
    ),
    "SSU": LocusConfig(
        name="SSU", group="fungi",
        description="Nuclear small subunit rDNA (18S)",
        species_identity=99.8, genus_identity=98.0, min_query_coverage=85.0,
        entrez_query='"Fungi"[Organism]',
        notes="Very conserved; higher-rank placement only.",
    ),
    "TEF1": LocusConfig(
        name="TEF1", group="fungi",
        description="Translation elongation factor 1-alpha",
        species_identity=99.0, genus_identity=95.0, min_query_coverage=80.0,
        entrez_query='"Fungi"[Organism]',
        notes="High resolution; intron-rich, so alignment needs care. "
              "For Fusarium prefer the curated FUSARIUM-ID / Fusarioid-ID databases.",
    ),
    "RPB2": LocusConfig(
        name="RPB2", group="fungi",
        description="RNA polymerase II second largest subunit",
        species_identity=99.0, genus_identity=95.0, min_query_coverage=80.0,
        entrez_query='"Fungi"[Organism]',
        notes="Protein-coding, single copy; align as codons where possible.",
    ),
    "BTUB": LocusConfig(
        name="BTUB", group="fungi",
        description="Beta-tubulin (tub2 / benA)",
        species_identity=99.0, genus_identity=95.0, min_query_coverage=80.0,
        entrez_query='"Fungi"[Organism]',
        notes="Standard secondary barcode for Penicillium/Aspergillus (benA).",
    ),
    "16S": LocusConfig(
        name="16S", group="bacteria",
        description="Bacterial/archaeal small subunit rRNA gene",
        species_identity=98.7, genus_identity=94.5, min_query_coverage=85.0,
        entrez_query='"Bacteria"[Organism] OR "Archaea"[Organism]',
        notes="98.7 % is the screening threshold below which DNA-DNA hybridisation "
              "is unnecessary (Kim et al. 2014). Above it, species assignment still "
              "needs genome-level data (ANI) or a type-strain database (EzBioCloud).",
    ),
    "RPOB": LocusConfig(
        name="RPOB", group="bacteria",
        description="RNA polymerase beta subunit, bacterial secondary marker",
        species_identity=97.7, genus_identity=94.0, min_query_coverage=80.0,
        entrez_query='"Bacteria"[Organism]',
        notes="Higher resolution than 16S within many genera.",
    ),
}


# ---- Sanger QC defaults -------------------------------------------------
QC_DEFAULTS = dict(
    mott_cutoff=0.05,      # Phred/Mott trimming error-probability cut-off (Q13)
    min_trimmed_length=200,
    min_mean_quality=20.0,
    max_ambiguity_fraction=0.02,
)

# ---- Alignment / trimming defaults -------------------------------------
ALIGN_DEFAULTS = dict(
    min_occupancy=0.50,    # drop columns present in < 50 % of sequences
    max_gap_fraction=0.50,
    min_block_length=3,    # discard isolated retained columns shorter than this
)

# ---- Tree defaults ------------------------------------------------------
TREE_DEFAULTS = dict(
    bootstrap_replicates=1000,   # IQ-TREE ultrafast bootstrap
    model=None,                  # None -> ModelFinder chooses by BIC
)
