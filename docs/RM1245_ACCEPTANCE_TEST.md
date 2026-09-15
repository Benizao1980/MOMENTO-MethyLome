# RM1245 real-data acceptance test

This note records the first real-data validation of `momento import-pacbio` against the legacy PacBio methylome for *Campylobacter jejuni* RM1245.

## Why this test matters

The initial importer correctly parsed synthetic GFFs but treated every motif-labelled modification call as equivalent. RM1245 showed that this is not sufficient for older PacBio/MotifMaker outputs.

## Observations from RM1245

The primary GFF contains motif assignments such as `RAATTY`, `CATG`, and `AGTNNNNNNRTTG`, together with a 41-bp `context` string and `identificationQv` values for many `m6A` calls.

For `RAATTY`, the called base occurs at two motif-relative positions:

| motif position | m6A calls | median identificationQv | calls with identificationQv >= 80 |
| ---: | ---: | ---: | ---: |
| 2 | 8,484 | 14 | 35 |
| 3 | 27,081 | 256 | 26,951 |

Position 3 is therefore the clear cognate methylated position. The published RM1245 methylome reported 26,952 methylated `RAATTY` sites, so the current GFF reproduces that result to within one call when the cognate position is used with an `identificationQv >= 80` diagnostic threshold.

For the other published motifs:

- `CATG`: all observed calls map to motif position 2; 6,149 `m6A` calls are present in the current GFF.
- `AGTNNNNNNRTTG`: all observed calls map to motif position 1; 312 `m6A` calls are present in the current GFF.

The current archived files are therefore not byte-for-byte identical to the processed files used for the published table, so validation should focus on reproducing the methylation semantics rather than requiring exact equality to every published count.

## Importer behaviour derived from this test

MOMENTO now:

1. preserves the raw GFF motif assignment;
2. parses the PacBio `context` field;
3. determines the called base's motif-relative position;
4. infers the cognate modified position conservatively from the data;
5. optionally filters on `identificationQv`;
6. retains excluded calls in the site-level table with explicit QC flags;
7. counts strand-specific target opportunities for self-reverse-complementary motifs once a cognate position is known.

There is deliberately **no universal default `identificationQv` threshold**. Historical PacBio pipelines differ, so thresholds must remain explicit and auditable until validated across more isolates.

## Next validation

Run the same workflow on RM1477, then test several biologically and technically distinct isolates before cohort-wide processing.
