# MOMENTO

**MOMENTO** is a species-agnostic toolkit for comparative bacterial methylomics.

It is designed for studies in which methylation has already been called by a sequencing-platform-specific tool (for example PacBio SMRT Analysis or Oxford Nanopore modified-base calling) and the biological question is comparative:

- Which methylation motifs are conserved or variable across isolates?
- Are there reproducible methylation states ("methylotypes")?
- Which methyltransferases (MTases) or restriction-modification systems explain those methylotypes?
- Are phase-variable loci associated with switching methylation states?
- Does methylation associate with host, ecology, disease, phenotype, or population structure?
- Does methylation add predictive information beyond the genome itself?

MOMENTO deliberately does **not** aim to replace raw-signal methylation callers. Instead, it provides a common representation and reproducible downstream analysis layer for bacterial population methylomics.

The name is a subtle nod to *Memento*: bacterial methylation can encode a molecular state or "memory", while the software reconstructs that state across populations.

## Design principles

1. **Species agnostic** — no organism-specific assumptions in the core package.
2. **Platform aware** — PacBio and ONT inputs are imported separately, then harmonised.
3. **Long-table first** — all importers produce one tidy internal schema.
4. **Population-aware** — host/phenotype association should account for genomic lineage.
5. **Interpretable** — motifs, MTases, phase-variable loci, genes and genomic positions remain traceable.
6. **Reproducible** — analyses should be rerunnable from primary caller outputs, not hand-edited spreadsheets.
7. **Teaching friendly** — examples explain *why* each step is performed, not only which command to type.

## Internal data model

The central table has one row per isolate × motif summary, with fields such as `sample_id`, `platform`, `motif`, `modification`, `called_sites`, `genomic_sites`, `methylated_fraction`, `median_score`, `source_file`, and `qc_status`.

## Intended workflow

```text
PacBio motif GFF / modification CSV ─┐
                                    ├─> MOMENTO import
ONT modBAM / bedMethyl ─────────────┘
                                           |
                                           v
                                harmonised long table
                                           |
                    ┌──────────────────────┼──────────────────────┐
                    v                      v                      v
                   QC                 motif matrix          annotation
                    |                      |                 MTase / phase
                    └──────────────┬───────┴───────────────┬─────┘
                                   v                       v
                              methylotypes            associations
                          PCA / PCoA / clustering    host / phenotype
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
momento --help
```

## First working importer: PacBio

`momento import-pacbio` converts a PacBio methylation-site GFF plus the matching assembly FASTA into MOMENTO's canonical motif summary. The assembly is used to count genomic methylatable target opportunities, so raw called-site counts can be converted to methylated fractions.

```bash
momento import-pacbio \
  --gff RM1245motifs.gff \
  --fasta CjRM1245.fasta \
  --sample RM1245 \
  --output RM1245.momento.tsv \
  --sites-output RM1245.sites.tsv
```

The importer supports IUPAC motifs such as `RAATTY`, handles reverse complements and older comma-encoded MotifMaker fields, and preserves the original GFF attributes in the site-level output.

### Why motif position matters

Historical PacBio/MotifMaker GFFs can annotate more than one modified base inside the same recognition motif. A motif label alone therefore does **not** prove that the called base is the cognate methyltransferase target.

When a GFF contains a `context` attribute, MOMENTO locates the called base at the centre of that context, determines its 1-based position within the motif, and infers the cognate modified position for each motif/modification combination. If multiple positions are observed, MOMENTO only selects one when the `identificationQv` evidence clearly separates it from the alternatives; otherwise it leaves the position unresolved rather than guessing.

All site assignments remain in `--sites-output`, including excluded calls. The site table records `motif_position`, `cognate_position`, `identification_qv` and explicit pass/fail flags so the filtering decision is auditable.

### Optional identification-QV filtering

There is deliberately **no universal default** for `identificationQv`, because historical PacBio pipelines and reprocessing conventions differ. A threshold can be supplied explicitly:

```bash
momento import-pacbio \
  --gff RM1245motifs.gff \
  --fasta CjRM1245.fasta \
  --sample RM1245 \
  --min-identification-qv 80 \
  --output RM1245.momento.tsv \
  --sites-output RM1245.sites.tsv
```

`80` is a development/validation setting for the RM1245 example, **not** a recommended universal PacBio threshold. Use `--keep-all-motif-positions` to disable cognate-position filtering for diagnostics or legacy comparisons.

## CLI roadmap

```bash
momento validate
momento import-pacbio   # implemented
momento import-ont      # planned
momento qc
momento matrix
momento cluster
momento mtase           # planned
momento phase           # planned
momento associate       # planned
momento report          # planned
```

The MVP now implements validation, context-aware PacBio GFF import, genomic target-opportunity counting, matrix generation, basic core-motif QC, prevalence summaries, PCA and Jaccard utilities.

## v0.1 acceptance test

The first real-data validation target is the published methylome of the *Campylobacter jejuni* HS:19 strains **RM1245** and **RM1477**, using primary PacBio methylation GFFs plus assembly FASTAs.

For RM1245, the published methylated/total counts include:

```text
RAATTY             26,952 / 27,594
CATG                6,207 /  6,266
AGTNNNNNNRTTG         310 /    316
```

Historical GFFs may not be byte-for-byte identical to the processed files used for a publication. The acceptance criterion is therefore **semantic reproduction**: recover the cognate modified position, confidence behaviour and approximately the published motif occupancy, and investigate rather than hide version-dependent discrepancies.

## Campylobacter teaching example

See `examples/campylobacter/README.md`. The example uses synthetic/public-safe data and explains:

- why a near-core motif such as `RAATTY` should be handled separately;
- how to distinguish failed methylome runs from genuine biological absence;
- how to define methylotypes;
- how to overlay host and genomic lineage;
- how to connect motifs to candidate MTases and phase-variable loci;
- how to avoid phylogenetic leakage in phenotype prediction.

For the current 57-genome Campylobacter panel, the core-genome phylogeny should be built separately from assembly FASTAs using a bacterial annotation workflow (for example Bakta or Prokka) followed by Panaroo/PIRATE and IQ-TREE. The PacBio `*motifs.gff` files are methylation-call GFFs, **not** gene-annotation GFFs and should not be supplied directly to Panaroo/PIRATE.

## Data policy

Do not commit raw reads, large modification CSVs, unpublished private assemblies, or sensitive metadata. Keep primary data in Box/HPC and commit code, schemas, public-safe derived tables, synthetic examples, documentation, and release-cleared figures.
