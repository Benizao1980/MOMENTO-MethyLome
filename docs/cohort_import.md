# Manifest-driven PacBio cohort import

`momento build-cohort` imports multiple PacBio methylation GFFs using an explicit TSV manifest. This avoids fragile filename globbing and makes assembly/GFF provenance auditable.

## Manifest

Required columns are:

```text
sample_id	fasta	gff
RM1245	01_GENOMES_FASTA/CjRM1245.fasta	02_MOTIF_GFF/RM1245motifs.gff
RM1477	01_GENOMES_FASTA/CjRM1477.fasta	02_MOTIF_GFF/RM1477motifs.gff
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

Two metrics are recorded:

- `context_coverage = context_compared / context_usable`: how many context-bearing GFF records can be mapped to the nominated FASTA sequence IDs and coordinates;
- `context_concordance = context_exact / context_compared`: how many compared contexts match the assembly exactly.

Defaults are:

```text
minimum context coverage      0.95
minimum context concordance   0.99
```

These can be changed with `--min-context-coverage` and `--min-context-concordance`.

Low coverage is typical of a sequence-ID mismatch or the wrong assembly. Low concordance is typical of a different assembly version. Duplicate FASTA record IDs are a hard preflight failure rather than being silently overwritten.

If a GFF contains no usable context attributes, preflight is recorded as `NO_CONTEXT` and import continues; this keeps older PacBio files usable while making the missing sequence-level validation explicit.

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

`cohort.momento.tsv` is the combined canonical long table for successfully imported samples.

`cohort.qc.tsv` contains one row per manifest sample with input paths, context metrics, final QC status, motif-row count, total site assignments, included site assignments and any failure/warning message.

Sample-level errors do not abort the whole cohort. They are recorded as `FAIL` in the audit table and processing continues. Very small but technically parseable imports are labelled `LOW_DATA` using the same configurable threshold as `import-pacbio`.

## Provenance corrections

Do not silently edit archived primary inputs. If a FASTA has a recoverable metadata/header problem, create a corrected derived copy and point the manifest to that file, retaining a note in project provenance. The underlying biological sequence should remain unchanged.
