# MOMENTO figure style: Life Aquatic

This is the default visual language for the Campylobacter methylome figures.
It is inspired by the cool, highly composed palette of *The Life Aquatic* while
remaining suitable for scientific publication.

## Core rules

- White figure and axes backgrounds (`#FFFFFF`).
- Cool blue / teal colours dominate.
- Coral is an accent only: use it for one highlighted class, feature or warning.
- Text is charcoal rather than pure black.
- Thin spines and dendrogram branches; avoid heavy boxes.
- Grid lines are omitted unless they materially aid interpretation.
- Prefer balanced, ordered layouts with generous whitespace.
- Export publication figures as PDF/SVG and a 300-dpi PNG preview.
- Never encode the same biological variable with two unrelated colour systems
  within a figure set.
- Missing / unknown metadata are pale grey and must remain visually distinct
  from true biological categories.

## Palette

| Name | Hex | Suggested use |
|---|---|---|
| deep_ocean | `#2F5D7E` | primary data / human |
| slate_blue | `#5E7E9B` | secondary category |
| muted_teal | `#5E8B8C` | environmental / general teal |
| pale_aqua | `#A8D5D1` | light category / annotation |
| seafoam | `#CFE8E5` | pale fills |
| steel_blue | `#7896AB` | additional cool category |
| powder_blue | `#B8CFDC` | additional cool category |
| sand | `#E8D9C5` | neutral warm contrast |
| coral | `#D97A6C` | sparing highlight / m4C accent |
| charcoal | `#2B2B2B` | text and fine outlines |
| light_grey | `#D9D9D9` | unknown / missing metadata |
| white | `#FFFFFF` | background |

## Preferred semantic colours

These mappings are defaults, not claims about biology. New categories should be
assigned from the cool palette before adding new warm colours.

| Category | Colour |
|---|---|
| human | deep_ocean |
| poultry / chicken | muted_teal |
| wild_bird | pale_aqua |
| black_bear | slate_blue |
| environment / water | seafoam |
| reference / laboratory | sand |
| unknown / missing | light_grey |
| m6A | deep_ocean |
| m4C | coral |

## Plot-specific guidance

### Accessory methylome heatmap

- Cluster isolates by Jaccard distance / average linkage.
- Keep the sample dendrogram in charcoal or blue-grey.
- Use white for absence and deep ocean blue for presence.
- Put categorical metadata in narrow annotation strips adjacent to the rows.
- Keep RAATTY out of the accessory matrix; show the core backbone separately.
- Use a small right-hand bar for the number of accessory features per isolate.

### PCoA

- White background, no grid.
- Colour points by the primary metadata variable (usually host group).
- Use marker shape for a second variable only when it remains legible.
- Use charcoal point edges and restrained alpha.
- Labels should be optional; use them for selected isolates rather than all
  points in the publication version.

### Prevalence and backbone plots

- Horizontal bars are preferred when motif labels are long.
- m6A is blue; m4C can use the coral accent.
- RAATTY is treated as the conserved backbone and should not be visually mixed
  into the accessory repertoire.

## Reproducibility

Figure scripts should take explicit input/output paths, avoid hard-coded HPC
locations, and write all outputs into a specified directory. Styling should be
implemented in code so re-running the analysis reproduces the same figures.
