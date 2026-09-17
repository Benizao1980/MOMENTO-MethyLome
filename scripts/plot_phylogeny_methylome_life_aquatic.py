#!/usr/bin/env python3
"""Plot a core-genome phylogeny with aligned MOMENTO methylome tracks.

Inputs are deliberately generic: an IQ-TREE-style Newick tree, a MOMENTO
methylotype directory, the cohort summary and isolate metadata. Extra tree tips
are pruned; every QC-clean methylotype sample must be present in the tree.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from momento.tree import TreeNode, ladderize, prune_tree, read_newick, tip_names, tree_layout


COLORS = {
    "deep_ocean": "#2F5D7E",
    "slate_blue": "#5E7E9B",
    "muted_teal": "#5E8B8C",
    "pale_aqua": "#A8D5D1",
    "seafoam": "#CFE8E5",
    "steel_blue": "#7896AB",
    "powder_blue": "#B8CFDC",
    "sand": "#E8D9C5",
    "coral": "#D97A6C",
    "charcoal": "#2B2B2B",
    "light_grey": "#D9D9D9",
    "white": "#FFFFFF",
}

COOL_CYCLE = [
    COLORS["deep_ocean"], COLORS["muted_teal"], COLORS["slate_blue"],
    COLORS["pale_aqua"], COLORS["steel_blue"], COLORS["seafoam"],
    COLORS["powder_blue"], COLORS["sand"], COLORS["coral"],
]

HOST_COLORS = {
    "human": COLORS["deep_ocean"],
    "poultry": COLORS["muted_teal"],
    "chicken": COLORS["muted_teal"],
    "wild_bird": COLORS["pale_aqua"],
    "black_bear": COLORS["slate_blue"],
    "environment": COLORS["seafoam"],
    "water": COLORS["seafoam"],
    "reference": COLORS["sand"],
    "unknown": COLORS["light_grey"],
    "": COLORS["light_grey"],
}

MOD_COLORS = {"m6A": COLORS["deep_ocean"], "m4C": COLORS["coral"]}


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": COLORS["white"],
            "axes.facecolor": COLORS["white"],
            "savefig.facecolor": COLORS["white"],
            "text.color": COLORS["charcoal"],
            "axes.labelcolor": COLORS["charcoal"],
            "axes.edgecolor": COLORS["charcoal"],
            "xtick.color": COLORS["charcoal"],
            "ytick.color": COLORS["charcoal"],
            "font.size": 8,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def normalise_category(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    text = str(value).strip()
    return text if text else "unknown"


def category_color_map(column: str, values: pd.Series) -> dict[str, str]:
    categories = sorted({normalise_category(x) for x in values})
    if column == "host_group":
        return {
            category: HOST_COLORS.get(category.lower(), COOL_CYCLE[i % len(COOL_CYCLE)])
            for i, category in enumerate(categories)
        }
    mapping: dict[str, str] = {}
    nonmissing = [x for x in categories if x.lower() not in {"unknown", "nan", "none"}]
    for i, category in enumerate(nonmissing):
        mapping[category] = COOL_CYCLE[i % len(COOL_CYCLE)]
    for category in categories:
        mapping.setdefault(category, COLORS["light_grey"])
    return mapping


def feature_modification(feature: str) -> str:
    return feature.split("|", 1)[0] if "|" in feature else "unknown"


def save_figure(fig: plt.Figure, stem: Path, formats: tuple[str, ...]) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight", "facecolor": COLORS["white"]}
        if fmt == "png":
            kwargs["dpi"] = 300
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)


def draw_tree(ax, node: TreeNode, coords: dict[int, tuple[float, float]]) -> None:
    x, _ = coords[id(node)]
    if node.children:
        ys = [coords[id(child)][1] for child in node.children]
        ax.plot([x, x], [min(ys), max(ys)], color=COLORS["slate_blue"], lw=0.8)
        for child in node.children:
            cx, cy = coords[id(child)]
            ax.plot([x, cx], [cy, cy], color=COLORS["slate_blue"], lw=0.8)
            draw_tree(ax, child, coords)


def load_raatty(summary_path: Path, samples: list[str]) -> pd.Series:
    df = pd.read_csv(summary_path, sep="\t")
    work = df.loc[
        (df["sample_id"].astype(str).isin(samples))
        & (df["motif"].astype(str) == "RAATTY")
        & (df["modification"].astype(str) == "m6A")
    ].copy()
    work["methylated_fraction"] = pd.to_numeric(work["methylated_fraction"], errors="coerce")
    dup = work["sample_id"].astype(str).duplicated()
    if dup.any():
        raise ValueError("Multiple RAATTY m6A summary rows for: " + ", ".join(work.loc[dup, "sample_id"].astype(str)))
    return work.set_index(work["sample_id"].astype(str))["methylated_fraction"].reindex(samples)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot phylogeny-aligned MOMENTO methylome tracks")
    parser.add_argument("--tree", required=True, type=Path, help="Newick tree, e.g. IQ-TREE .treefile")
    parser.add_argument("--methylotype-dir", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--cohort-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--annotations", default="host_group,provenance_class")
    parser.add_argument("--min-feature-prevalence", type=int, default=2)
    parser.add_argument("--max-features", type=int, default=30, help="0 means no maximum")
    parser.add_argument("--formats", default="png,pdf,svg")
    parser.add_argument("--no-ladderize", action="store_true")
    args = parser.parse_args()

    configure_matplotlib()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    formats = tuple(x.strip().lower() for x in args.formats.split(",") if x.strip())
    annotations = tuple(x.strip() for x in args.annotations.split(",") if x.strip())

    presence = pd.read_csv(args.methylotype_dir / "accessory_presence_absence.tsv", sep="\t", index_col=0)
    presence.index = presence.index.astype(str)
    metadata = pd.read_csv(args.metadata, sep="\t")
    if "sample_id" not in metadata.columns:
        raise SystemExit("metadata must contain sample_id")
    if metadata["sample_id"].astype(str).duplicated().any():
        raise SystemExit("metadata contains duplicate sample_id values")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id")

    missing_meta = sorted(set(presence.index) - set(metadata.index))
    if missing_meta:
        raise SystemExit("metadata missing samples: " + ", ".join(missing_meta))

    raw_tree = read_newick(args.tree)
    tree_tips = set(tip_names(raw_tree))
    samples = list(presence.index)
    missing_tree = sorted(set(samples) - tree_tips)
    extra_tree = sorted(tree_tips - set(samples))
    if missing_tree:
        raise SystemExit("tree missing MOMENTO samples: " + ", ".join(missing_tree))

    tree = prune_tree(raw_tree, samples)
    if not args.no_ladderize:
        tree = ladderize(tree)
    coords, tip_order = tree_layout(tree)

    presence = presence.loc[tip_order]
    metadata = metadata.loc[tip_order]
    raatty = load_raatty(args.cohort_summary, tip_order)

    prevalence = presence.sum(axis=0).astype(int)
    selected = [x for x in presence.columns if prevalence[x] >= args.min_feature_prevalence]
    selected = sorted(selected, key=lambda x: (-int(prevalence[x]), feature_modification(str(x)), str(x)))
    if args.max_features > 0:
        selected = selected[: args.max_features]
    if not selected:
        raise SystemExit("No accessory features remain after display filtering")

    n = len(tip_order)
    k = len(selected)
    valid_annotations = [x for x in annotations if x in metadata.columns]
    fig_w = max(15.0, 10.5 + 0.23 * k + 0.45 * len(valid_annotations))
    fig_h = max(9.0, 0.19 * n + 1.8)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=COLORS["white"])
    gs = fig.add_gridspec(
        1,
        6,
        width_ratios=[4.2, max(0.7, 0.38 * len(valid_annotations)), 1.35, max(4.5, 0.22 * k), 1.15, 2.2],
        wspace=0.05,
    )
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_ann = fig.add_subplot(gs[0, 1])
    ax_raatty = fig.add_subplot(gs[0, 2])
    ax_heat = fig.add_subplot(gs[0, 3])
    ax_count = fig.add_subplot(gs[0, 4])
    ax_leg = fig.add_subplot(gs[0, 5])

    # Tree
    draw_tree(ax_tree, tree, coords)
    xmax = max(x for x, _ in coords.values())
    label_pad = max(0.02 * xmax, 0.002) if xmax > 0 else 0.15
    for tip in [x for x in _walk(tree) if x.is_tip]:
        x, y = coords[id(tip)]
        ax_tree.text(x + label_pad, y, str(tip.name), va="center", ha="left", fontsize=6.5)
    ax_tree.set_ylim(-0.7, n - 0.3)
    ax_tree.invert_yaxis()
    ax_tree.set_yticks([])
    ax_tree.spines["left"].set_visible(False)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.set_xlabel("Core-genome branch length", fontsize=7)
    if xmax > 0:
        ax_tree.set_xlim(0, xmax + max(0.22 * xmax, 6 * label_pad))

    # Metadata strips
    ann_maps: dict[str, dict[str, str]] = {}
    if valid_annotations:
        rgb = np.zeros((n, len(valid_annotations), 3), dtype=float)
        for j, column in enumerate(valid_annotations):
            cmap = category_color_map(column, metadata[column])
            ann_maps[column] = cmap
            for i, value in enumerate(metadata[column]):
                rgb[i, j, :] = to_rgb(cmap[normalise_category(value)])
        ax_ann.imshow(rgb, aspect="auto", interpolation="nearest", origin="upper")
        ax_ann.set_xticks(np.arange(len(valid_annotations)))
        ax_ann.set_xticklabels(valid_annotations, rotation=90, fontsize=6.5)
        ax_ann.set_yticks([])
    else:
        ax_ann.axis("off")
    for spine in ax_ann.spines.values():
        spine.set_visible(False)

    # RAATTY occupancy
    y = np.arange(n)
    ax_raatty.scatter(raatty.to_numpy(), y, s=16, color=COLORS["deep_ocean"], edgecolors=COLORS["charcoal"], linewidths=0.3)
    ax_raatty.axvline(0.99, color=COLORS["powder_blue"], lw=0.7, ls="--")
    finite = raatty.dropna()
    xmin = max(0.0, min(0.94, float(finite.min()) - 0.005)) if len(finite) else 0.94
    ax_raatty.set_xlim(xmin, 1.002)
    ax_raatty.set_ylim(-0.7, n - 0.3)
    ax_raatty.invert_yaxis()
    ax_raatty.set_yticks([])
    ax_raatty.set_xlabel("RAATTY\noccupancy", fontsize=7)
    ax_raatty.tick_params(axis="x", labelsize=6)
    ax_raatty.spines["left"].set_visible(False)
    ax_raatty.spines["top"].set_visible(False)
    ax_raatty.spines["right"].set_visible(False)

    # Accessory methylome display tracks
    display = presence[selected]
    heat_rgb = np.ones((n, k, 3), dtype=float)
    for j, feature in enumerate(selected):
        color = to_rgb(MOD_COLORS.get(feature_modification(str(feature)), COLORS["deep_ocean"]))
        hits = display.iloc[:, j].to_numpy(dtype=bool)
        heat_rgb[hits, j, :] = color
    ax_heat.imshow(heat_rgb, aspect="auto", interpolation="nearest", origin="upper")
    ax_heat.set_xticks(np.arange(k))
    ax_heat.set_xticklabels(selected, rotation=90, fontsize=5.5)
    ax_heat.set_yticks([])
    ax_heat.set_xlabel("Accessory motif family", fontsize=7)
    ax_heat.set_title(
        f"Accessory methylome (top/common {k} features)", fontsize=10, pad=9
    )
    ax_heat.tick_params(length=0)
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.45)

    # Total accessory repertoire size always uses all features.
    counts = presence.sum(axis=1).to_numpy(dtype=float)
    ax_count.barh(y, counts, height=0.68, color=COLORS["muted_teal"], edgecolor="none")
    ax_count.set_ylim(-0.7, n - 0.3)
    ax_count.invert_yaxis()
    ax_count.set_yticks([])
    ax_count.set_xlabel("all accessory\nfeatures", fontsize=7)
    ax_count.tick_params(axis="x", labelsize=6)
    ax_count.spines["top"].set_visible(False)
    ax_count.spines["right"].set_visible(False)
    ax_count.spines["left"].set_visible(False)

    # Legends
    ax_leg.axis("off")
    legend_y = 0.98
    state_handles = [
        Patch(facecolor=MOD_COLORS["m6A"], label="m6A presence"),
        Patch(facecolor=MOD_COLORS["m4C"], label="m4C presence"),
        Patch(facecolor=COLORS["white"], edgecolor=COLORS["charcoal"], linewidth=0.5, label="absence"),
    ]
    leg = ax_leg.legend(handles=state_handles, title="motif state", frameon=False, loc="upper left", bbox_to_anchor=(0, legend_y), fontsize=7, title_fontsize=7)
    ax_leg.add_artist(leg)
    legend_y -= 0.22
    for column, cmap in ann_maps.items():
        if len(cmap) > 12:
            continue
        handles = [Patch(facecolor=color, edgecolor="none", label=cat) for cat, color in cmap.items()]
        leg = ax_leg.legend(handles=handles, title=column, frameon=False, loc="upper left", bbox_to_anchor=(0, legend_y), fontsize=6.5, title_fontsize=7)
        ax_leg.add_artist(leg)
        legend_y -= min(0.34, 0.07 + 0.038 * len(handles))

    fig.suptitle("Core-genome phylogeny and methylome repertoire", fontsize=13, y=0.995)
    save_figure(fig, args.output_dir / "figure_phylogeny_methylome", formats)
    plt.close(fig)

    # Audits for reproducibility.
    pd.DataFrame({"sample_id": tip_order}).to_csv(args.output_dir / "phylogeny_tip_order.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "feature": selected,
            "prevalence": [int(prevalence[x]) for x in selected],
            "modification": [feature_modification(str(x)) for x in selected],
        }
    ).to_csv(args.output_dir / "phylogeny_feature_order.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "category": ["methylotype_samples", "tree_tips_before_pruning", "extra_tree_tips_pruned", "tree_tips_plotted"],
            "n": [len(samples), len(tree_tips), len(extra_tree), len(tip_order)],
            "values": [";".join(sorted(samples)), ";".join(sorted(tree_tips)), ";".join(extra_tree), ";".join(tip_order)],
        }
    ).to_csv(args.output_dir / "phylogeny_alignment_audit.tsv", sep="\t", index=False)

    print(
        f"PHYLOGENY FIGURE: {len(tip_order)} samples; {k}/{presence.shape[1]} accessory features displayed; "
        f"{len(extra_tree)} extra tree tips pruned -> {args.output_dir}"
    )


def _walk(node: TreeNode):
    yield node
    for child in node.children:
        yield from _walk(child)


if __name__ == "__main__":
    main()
