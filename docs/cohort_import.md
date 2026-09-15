# Manifest-driven PacBio cohort import

`momento build-cohort` imports multiple PacBio methylation GFFs using an explicit TSV manifest. This avoids fragile filename globbing and makes assembly/GFF provenance auditable.

## Manifest

Required columns are `sample_id`, `fasta` and `gff`. An optional `note` column is copied into the cohort audit table for provenance decisions such as choosing an assembly version or using a header-corrected derived FASTA.

```text
sample_id	fasta	gff	note
RM1245	01_GENOMES_FASTA/CjRM1245.fasta	02_MOTIF_GFF/RM1245motifs.gff	
RM3427	01_GENOMES_FASTA/MOMENTO_CORRECTED/CjRM3427.momento.fasta	02_MOTIF_GFF/RM3427motifs.gff	corrected duplicate FASTA record IDs; sequence unchanged
```

Relative paths are resolved relative to the manifest file. `sample_id` values must be unique.

Run:

```bash
momento build-cohort \
  --manifest samples.tsv \
  --output-dir MOMENTO_COHORT
```

There is deliberately no default `identificationQv` cutoff. Historical PacBio processing regimes can have very different QV distributions. If a threshold is needed for a validated legacy comparison, it must be supplied explicitly with `--min-identification-qv`.

## FASTA/GFF context preflight

Before importing a sample, MOMENTO compares odd-length PacBio `context` strings with the nominated assembly FASTA.

Two sequence-identity metrics are recorded:

- `context_coverage = context_compared / context_usable`: how many context-bearing GFF records can be mapped to the nominated FASTA sequence IDs and coordinates;
- `context_concordance = context_exact / context_compared`: how many compared contexts match the assembly exactly.

Defaults are:

```text
minimum context coverage      0.95
minimum context concordance   0.99
```

These can be changed with `--min-context-coverage` and `--min-context-concordance`.

Low coverage is typical of a sequence-ID mismatch or the wrong assembly. Low concordance is typical of a different assembly version. Duplicate FASTA record IDs are a hard preflight failure rather than being silently overwritten.

If a GFF contains no usable context attributes, sequence-context preflight is recorded as `NO_CONTEXT` unless another preflight condition gives a more informative status; import can still continue for legacy files.

## Partial/truncated GFF detection

A GFF can match its FASTA perfectly yet represent only part of the genome. If the whole FASTA is then used as the denominator, unassessed sequence can be mistaken for biologically unmethylated sequence. MOMENTO therefore records coordinate-coverage metrics using all GFF feature records:

- `gff_span_fraction`: genome-weighted min-to-max feature span across FASTA records;
- `gff_bin_coverage_fraction`: genome fraction in 20-bin windows containing at least one GFF feature.

By default, coordinate-span QC is enforced only for GFFs with at least 1,000 feature records. This avoids labelling genuinely sparse legacy inputs as truncated. The default thresholds are:

```text
minimum GFF coordinate span        0.90
minimum occupied-bin genome fraction 0.80
minimum records before span QC     1000
```

They can be changed with `--min-gff-span-fraction`, `--min-gff-bin-coverage` and `--min-gff-records-for-span`.

A sufficiently large GFF below either coordinate threshold is labelled `PARTIAL_GFF`. MOMENTO still writes that sample's motif and site tables for auditing, but propagates `qc_status=PARTIAL_GFF` so analyses restricted to `PASS` samples exclude it automatically.

This status is intended for cases such as a methylation GFF that stops abruptly part-way through an otherwise matching chromosome. It should not be interpreted as a biological low-methylation state.

## Outputs

The output directory contains:

```text
cohort.momento.tsv
cohort.qc.tsv
samples/
  SAMPLE_A/
    SAMPLE_A.momento.tsv
    SAMPLE_A.sites.tsv
  SAMPLE_B/
    SAMPLE_B.momento.tsv
    SAMPLE_B.sites.tsv
```

`cohort.momento.tsv` is the combined canonical long table for successfully imported samples, including imported warning states such as `LOW_DATA` and `PARTIAL_GFF` with their QC labels preserved.

`cohort.qc.tsv` contains one row per manifest sample with input paths, optional manifest provenance note, context metrics, GFF coordinate-coverage metrics, final QC status, motif-row count, total site assignments, included site assignments and any failure/warning message.

Sample-level errors do not abort the whole cohort. They are recorded as `FAIL` in the audit table and processing continues. Very small but technically parseable imports are labelled `LOW_DATA` using the same configurable threshold as `import-pacbio`.

## Provenance corrections

Do not silently edit archived primary inputs. If a FASTA has a recoverable metadata/header problem, create a corrected derived copy and point the manifest to that file, retaining a note in project provenance. The underlying biological sequence should remain unchanged.
