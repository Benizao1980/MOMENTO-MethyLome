# Inputs

## Canonical methylation table

`sample_id, platform, motif, modification, called_sites, genomic_sites, methylated_fraction, median_score, source_file, qc_status`

For PacBio imports, `genomic_sites` is intended to represent methylatable target-site opportunities when a cognate modified position can be inferred. For self-reverse-complementary motifs this can differ from the number of physical motif strings in the assembly because the two DNA strands contain distinct methylatable targets.

## Sample metadata

Recommended: `sample_id, species, host, host_group, location, sampling_date, ST, clonal_complex, lineage, sequencing_platform`.

## PacBio

The first production route uses a motif/site GFF plus matching assembly FASTA. Large modification CSVs are optional richer inputs.

Historical PacBio/MotifMaker GFFs can contain:

- an explicit `motif` assignment;
- `context`, with the called base at the centre;
- `identificationQv` and other call-confidence fields;
- a GFF feature type such as `m6A`, `m4C` or `modified_base`;
- older motif encodings such as `RAATTY,2` or overlapping assignments such as `GCCAA,RAATTY`.

A motif label alone is not sufficient to establish that the called base is the cognate methyltransferase target. When `context` is present, MOMENTO records the called base's motif-relative position and conservatively infers the cognate position for each motif/modification combination.

The optional site-level output contains the original attributes plus:

`identification_qv, context, motif_position, cognate_position, cognate_position_source, passes_score, passes_identification_qv, passes_motif_position, included_in_summary`

Calls excluded from the motif summary remain in the site table for audit/QC.

`--min-identification-qv` has no default because historical PacBio pipelines can differ. Thresholds should be justified and validated for the dataset rather than assumed globally.

## ONT

First production route should target bedMethyl from a validated modBAM workflow; direct modBAM support can follow.

## Assemblies

Used for motif opportunities, genomic context, MTase annotation and candidate repeat detection. Important phase-state claims should be raw-read validated because homopolymer lengths can be assembly-sensitive.
