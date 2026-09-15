# PacBio real-data acceptance tests

This note records real-data validation of `momento import-pacbio` against legacy PacBio methylomes spanning published reference isolates, ecologically distinct *Campylobacter jejuni* isolates, *C. jejuni* subsp. *doylei*, and deliberately sparse/broken GFF inputs.

## Why these tests matter

The initial importer correctly parsed synthetic GFFs but treated every motif-labelled modification call as equivalent. RM1245 showed that this is not sufficient for older PacBio/MotifMaker outputs. Subsequent isolates test whether the revised logic generalises without isolate-specific code.

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

The current archived FASTA/GFF pair is not byte-for-byte identical to the processed files used for the published table, so validation should focus on reproducing methylation semantics rather than requiring exact equality to every historical denominator.

## RM1477: independent published-isolate validation

The same workflow was run without code changes on RM1477 using `CjRM1477.fasta`, `RM1477motifs.gff`, cognate-position filtering, and `identificationQv >= 80`.

| motif | MOMENTO methylated / genomic target sites | fraction | published methylated / total |
| --- | ---: | ---: | ---: |
| `RAATTY` | 27,524 / 27,646 | 0.9956 | 27,623 / 28,000 |
| `CATG` | 6,301 / 6,316 | 0.9976 | 6,377 / 6,404 |
| `AGTNNNNNNRTTG` | 322 / 322 | 1.0000 | 325 / 325 |

All three published motifs produce biologically valid fractions at or below 1.0 and closely recapitulate the historical methylome.

## RM10527: assembly provenance and generalisation

RM10527 had two candidate assemblies. GFF-context concordance resolved the correct input unambiguously:

- `CjRM10527.fasta`: 18,436 exact context matches, 46,852 mismatches, 2 out-of-range calls;
- `CjRM10527-PacBio-MiSeq.fasta`: 65,288 exact context matches, 0 mismatches, 2 out-of-range calls.

Using the PacBio-MiSeq assembly, the importer produced ten coherent motif/modification rows. `RAATTY` was 26,225 / 26,728 (0.9812); all motif fractions were <= 1.0. This establishes GFF-to-FASTA context concordance as an important future cohort preflight.

## SKBC41: ecologically distinct black-bear isolate

SKBC41 imported cleanly with seven motif/modification rows and 34,462 / 34,705 site-by-motif assignments included at the diagnostic QV80 setting. `RAATTY` was 27,131 / 27,220 (0.9967), while the accessory motif repertoire differed substantially from RM10527. This supports a conserved core methylation backbone plus variable accessory methylation systems.

## RM4099: confidence thresholds are dataset-specific

RM4099 (*C. jejuni* subsp. *doylei*) exposed an important confidence-threshold issue. Applying `identificationQv >= 80` retained only 6,921 / 45,292 site-by-motif assignments and reduced `RAATTY` to 2,635 / 29,330. This was not biologically credible for the primary GFF.

Without an extra identification-QV threshold, cognate-position-aware import was coherent:

| motif | methylated / genomic target sites | fraction | median identificationQv |
| --- | ---: | ---: | ---: |
| `CAACA` | 2,926 / 2,929 | 0.9990 | 63 |
| `CCAAC` | 1,317 / 1,320 | 0.9977 | 55 |
| `CCNNNNNNNGG` (m4C) | 1,921 / 2,230 | 0.8614 | 26 |
| `GATC` | 7,349 / 7,354 | 0.9993 | 71 |
| `GATGAA` | 3,017 / 3,020 | 0.9990 | 64 |
| `RAATTY` | 28,613 / 29,330 | 0.9756 | 48 |

This establishes a core design rule: **MOMENTO must not impose a universal `identificationQv` threshold.** The default cohort workflow should preserve the caller's accepted calls, use cognate-position semantics, retain identification-QV for QC/sensitivity analysis, and apply thresholds only when explicitly requested.

Modification type also matters: m4C calls can have substantially lower identification-QV distributions than m6A calls and should not be judged by the same confidence profile without validation.

## Sparse and broken-input tests

Two deliberately tiny GFFs were tested:

- **RM3433:** six site-by-motif assignments produced one motif row (`ANNGCNNANANNGACNG`, 6 / 6). Parsing is technically valid, but the sample is far too sparse to enter cohort analysis silently. MOMENTO therefore flags imports below a configurable site-count threshold as `LOW_DATA` rather than calling them an ordinary pass.
- **RM3439:** the 461-byte GFF contains no feature records. The importer correctly rejects it. The CLI converts this to a concise `FAIL` message rather than exposing a Python traceback.

These are QC outcomes rather than biological exclusions: sparse files remain auditable, while structurally empty GFFs fail explicitly.

## Importer behaviour derived from these tests

MOMENTO now:

1. preserves raw GFF motif assignments and attributes;
2. parses PacBio `context` and `identificationQv`;
3. determines the called base's motif-relative position;
4. infers the cognate modified position conservatively from the data;
5. makes identification-QV filtering optional rather than universal;
6. retains excluded calls in the site-level table with explicit QC flags;
7. counts strand-specific target opportunities for self-reverse-complementary motifs once a cognate position is known;
8. warns on extremely sparse imports using a configurable `LOW_DATA` threshold;
9. fails cleanly on structurally empty GFFs.

## Acceptance status

- **RM1245: PASS** — published-isolate validation; `RAATTY` published methylated-site count reproduced to within one call under the diagnostic QV80 threshold.
- **RM1477: PASS** — independent published-isolate validation.
- **RM10527: PASS** — non-publication generalisation and assembly-provenance test.
- **SKBC41: PASS** — ecologically distinct isolate with a different accessory methylome.
- **RM4099: PASS with no extra QV threshold** — demonstrates why confidence thresholds must remain dataset-specific.
- **RM3433: LOW_DATA** — technically parseable but insufficient evidence for ordinary cohort inclusion.
- **RM3439: FAIL** — no GFF feature records.

This validation set is sufficient to move from manual isolate testing to a manifest-driven cohort preflight and batch importer.
