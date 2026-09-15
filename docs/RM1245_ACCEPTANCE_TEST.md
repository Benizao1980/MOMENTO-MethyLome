# PacBio real-data acceptance tests: RM1245 and RM1477

This note records the first real-data validation of `momento import-pacbio` against legacy PacBio methylomes for *Campylobacter jejuni* RM1245 and RM1477.

## Why these tests matter

The initial importer correctly parsed synthetic GFFs but treated every motif-labelled modification call as equivalent. RM1245 showed that this is not sufficient for older PacBio/MotifMaker outputs. RM1477 provides an independent validation that the revised logic generalises beyond one isolate.

## RM1245: discovering cognate motif positions

The primary GFF contains motif assignments such as `RAATTY`, `CATG`, and `AGTNNNNNNRTTG`, together with a 41-bp `context` string and `identificationQv` values for many `m6A` calls.

For `RAATTY`, the called base occurs at two motif-relative positions:

| motif position | m6A calls | median identificationQv | calls with identificationQv >= 80 |
| ---: | ---: | ---: | ---: |
| 2 | 8,484 | 14 | 35 |
| 3 | 27,081 | 256 | 26,951 |

Position 3 is therefore the clear cognate methylated position. The published RM1245 methylome reported 26,952 methylated `RAATTY` sites, so the current GFF reproduces that result to within one call when the cognate position is used with an `identificationQv >= 80` diagnostic threshold.

With cognate-position filtering and `identificationQv >= 80`, MOMENTO reports:

| motif | MOMENTO methylated / genomic target sites | fraction | published methylated / total |
| --- | ---: | ---: | ---: |
| `RAATTY` | 26,951 / 27,106 | 0.9943 | 26,952 / 27,594 |
| `CATG` | 6,144 / 6,154 | 0.9984 | 6,207 / 6,266 |
| `AGTNNNNNNRTTG` | 312 / 314 | 0.9936 | 310 / 316 |

The current archived FASTA/GFF pair is not byte-for-byte identical to the processed files used for the published table, so validation should focus on reproducing the methylation semantics rather than requiring exact equality to every historical denominator.

## RM1477: independent validation

The same workflow was run without code changes on RM1477 using `CjRM1477.fasta`, `RM1477motifs.gff`, cognate-position filtering, and `identificationQv >= 80`.

MOMENTO reports:

| motif | MOMENTO methylated / genomic target sites | fraction | published methylated / total |
| --- | ---: | ---: | ---: |
| `RAATTY` | 27,524 / 27,646 | 0.9956 | 27,623 / 28,000 |
| `CATG` | 6,301 / 6,316 | 0.9976 | 6,377 / 6,404 |
| `AGTNNNNNNRTTG` | 322 / 322 | 1.0000 | 325 / 325 |

All three core published motifs produce biologically valid fractions at or below 1.0 and closely recapitulate the historical methylome. The modest count differences are consistent with archived inputs or processing versions differing from the exact files used to generate the published table.

## Importer behaviour derived from these tests

MOMENTO now:

1. preserves the raw GFF motif assignment;
2. parses the PacBio `context` field;
3. determines the called base's motif-relative position;
4. infers the cognate modified position conservatively from the data;
5. optionally filters on `identificationQv`;
6. retains excluded calls in the site-level table with explicit QC flags;
7. counts strand-specific target opportunities for self-reverse-complementary motifs once a cognate position is known.

There is deliberately **no universal default `identificationQv` threshold**. Historical PacBio pipelines differ, so thresholds must remain explicit and auditable until validated across more isolates.

## Acceptance status

- **RM1245: PASS** — cognate-position-aware import removes spurious `RAATTY` calls and reproduces the published methylated-site count to within one call under the diagnostic QV80 threshold.
- **RM1477: PASS** — the same importer and threshold produce coherent, near-complete methylation for the three published motifs without isolate-specific code changes.

These two isolates are sufficient to accept the core PacBio importer logic for v0.1. The next validation stage should deliberately include biologically and technically distinct isolates (for example a Salinas isolate, a black-bear isolate, *C. jejuni* subsp. *doylei* RM4099, and known tiny/partial GFF cases) before cohort-wide analysis.
