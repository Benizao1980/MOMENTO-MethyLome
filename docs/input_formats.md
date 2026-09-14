# Inputs

## Canonical methylation table
`sample_id, platform, motif, modification, called_sites, genomic_sites, methylated_fraction, median_score, source_file, qc_status`

## Sample metadata
Recommended: `sample_id, species, host, host_group, location, sampling_date, ST, clonal_complex, lineage, sequencing_platform`.

## PacBio
First production route should use motif/site GFF plus assembly FASTA. FASTA provides the denominator of genomic motif opportunities. Large modification CSVs are optional richer inputs.

## ONT
First production route should target bedMethyl from a validated modBAM workflow; direct modBAM support can follow.

## Assemblies
Used for motif opportunities, genomic context, MTase annotation and candidate repeat detection. Important phase-state claims should be raw-read validated because homopolymer lengths can be assembly-sensitive.
