# Campylobacter metadata used by the MOMENTO population analysis

This directory contains the small, derived metadata tables that are safe and useful to version with the code. Large primary files and private study data remain outside the repository.

## Current canonical table

`campylobacter_methylotypes_metadata_v0.2.tsv` is the canonical biological metadata table for the 52 QC-PASS isolates in the current Campylobacter methylome/phylogeny analysis.

The table was rebuilt from the investigator-curated workbook `Methylome Metadata.xlsx` and then merged with the useful provenance fields already present in the earlier MOMENTO v0.1 metadata table. The source workbook is not committed here; the exact workbook used for this derivation had SHA256:

```text
a8eab63ca59c8ef3cbaf64ebbe3c89d670c1072a3a5dbdc990baf23716d44b56
```

This checksum is also recorded in the `metadata_source` column so that downstream result folders preserve the provenance of the biological annotations.

## What changed from v0.1

The v0.1 table was sufficient for plotting a first phylogeny but contained many placeholder `unknown` values for host, ST and clonal complex. The curated workbook resolves nearly all of those fields. In v0.2:

- all 52 isolates have an ST;
- 39/52 have a named clonal complex and 13 are recorded as `ND`;
- host is resolved for 51/52 isolates;
- country is resolved for 51/52 isolates;
- year is present for 41/52 isolates;
- the detailed host field is kept separate from the broader ecological `host_group` used for analysis;
- study/provenance fields from v0.1 are retained where they were already known.

The only isolate with unresolved host/geography/year in the curated workbook is `RM12657`.

## Ecological host grouping

The population analysis uses a deliberately small set of ecological host categories. The transformation from the source workbook is explicit rather than inferred during plotting:

| Source host | `host_detail` | `host_group` |
| --- | --- | --- |
| Human | human | human |
| Bear | black bear | wild_mammal |
| Crow | crow | wild_bird |
| Blackbird | blackbird | wild_bird |
| Starling | starling | wild_bird |
| Chicken | chicken | poultry |
| Chicken Meat | chicken meat | poultry |
| Turkey | turkey | poultry |
| Water | water | environment |
| missing | unknown | unknown |

Current host-group counts are:

```text
human        31
wild_mammal   8
wild_bird     6
poultry       5
environment   1
unknown       1
```

`host_group` is the variable intended for broad ecological comparisons. `host_detail` should be used when the exact source matters.

## ST and clonal-complex conventions

The source workbook uses labels such as `ST89` and `ST-45`. The repository table stores the underlying identifiers without the redundant prefixes:

```text
ST89   -> st = 89
ST-45  -> clonal_complex = 45
ND     -> clonal_complex = ND
```

The plotting script adds the human-readable `ST`/`CC` prefixes back when drawing figure labels. This keeps the data column machine-friendly while preserving familiar nomenclature in figures.

## Provenance is not ecology

`provenance_class` and `study` are retained separately from host. For example, `reference_panel`, `salinas_collection`, `black_bear_2025` and `heikema_2021` describe where an isolate entered the dataset, not its ecological host.

Do **not** substitute provenance for ecology in association tests. Host and provenance are partially confounded in this small panel, so provenance should be included as a sensitivity covariate rather than interpreted as a biological host effect.

## Known limitations

The curated workbook did not contain a species/subspecies field. The v0.2 table therefore preserves species/subspecies only where they were already explicitly recorded in the v0.1 metadata. It does not silently fill missing taxonomy from general knowledge.

Similarly, missing years and `ND` clonal complexes are kept as missing/`ND`; they are not imputed.

## Rebuilding the table

The v0.2 table is a deterministic curation step, not a statistical output. If the workbook is revised later:

1. keep the isolate identifier as the join key;
2. regenerate `host_detail` and `host_group` using the mapping above;
3. normalise ST and CC labels as described above;
4. carry forward existing provenance fields only when they refer to the same isolate;
5. record the new workbook checksum in `metadata_source`;
6. increment the metadata version rather than overwriting this table;
7. rerun the phylogeny figure and all host/lineage analyses from the new version.

Keeping metadata versions immutable makes it possible to reproduce exactly which annotations were used for each result.
