#!/usr/bin/env python3
"""Publication-style exploratory figures for MOMENTO methylotype outputs.

The script consumes files written by ``momento methylotypes`` plus an isolate
metadata TSV. It deliberately keeps plotting separate from the core analysis so
that visual choices cannot alter the numerical outputs.
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
from scipy.cluster.hierarchy import dendrogram


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
    COLORS["deep_ocean"],
    COLORS["muted_teal"],
    COLORS["slate_blue"],
    COLORS["pale_aqua"],
    COLORS["steel_blue"],
    COLORS["seafoam"],
    COLORS["powder_blue"],
    COLORS["sand"],
    COLORS["coral"],
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

MOD_COLORS = {
    "m6A": COLORS["deep_ocean"],
    "m4C": COLORS["coral"],
}

MARKERS = ["o", "s", "^", "D", "P", "v", "X", "<", ">"]


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
            "font.size": 9,
            "axes.titleweight": "normal",
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tsv(path: Path, *, index_col=None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=index_col, dtype=None)


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
        if category not in mapping:
            mapping[category] = COLORS["light_grey"]
    return mapping


def save_figure(fig: plt.Figure, stem: Path, formats: tuple[str, ...]) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        kwargs = {"bbox_inches": "tight", "facecolor": COLORS["white"]}
        if fmt.lower() == "png":
            kwargs["dpi"] = 300
        fig.savefig(stem.with_suffix(f".{fmt}"), **kwargs)


def load_inputs(methylotype_dir: Path, metadata_path: Path):
    presence = read_tsv(methylotype_dir / "accessory_presence_absence.tsv", index_col=0)
    pcoa = read_tsv(methylotype_dir / "pcoa.tsv", index_col=0)
    variance = read_tsv(methylotype_dir / "pcoa_variance.tsv")
    linkage = read_tsv(methylotype_dir / "linkage.tsv")
    prevalence = read_tsv(methylotype_dir / "motif_family_prevalence.tsv")
    metadata = read_tsv(metadata_path)

    if "sample_id" not in metadata.columns:
        raise ValueError("metadata must contain a sample_id column")
    if metadata["sample_id"].duplicated().any():
        dup = sorted(metadata.loc[metadata["sample_id"].duplicated(), "sample_id"].astype(str).unique())
        raise ValueError("duplicate sample_id values in metadata: " + ", ".join(dup))

    metadata = metadata.set_index("sample_id")
    missing = sorted(set(presence.index.astype(str)) - set(metadata.index.astype(str)))
    if missing:
        raise ValueError(
            "metadata is missing methylotype samples: " + ", ".join(missing)
        )
    metadata = metadata.reindex(presence.index)
    return presence, pcoa, variance, linkage, prevalence, metadata


def feature_modification(feature: str) -> str:
    return feature.split("|", 1)[0] if "|" in feature else "unknown"


def plot_heatmap(
    presence: pd.DataFrame,
    linkage_df: pd.DataFrame,
    metadata: pd.DataFrame,
    annotations: tuple[str, ...],
    outdir: Path,
    formats: tuple[str, ...],
) -> list[str]:
    z = linkage_df[["left", "right", "distance", "n_members"]].to_numpy(dtype=float)

    n_samples, n_features = presence.shape
    fig_width = max(16.0, n_features * 0.16 + 6.0)
    fig_height = max(9.0, n_samples * 0.18 + 2.0)
    fig = plt.figure(figsize=(fig_width, fig_height), facecolor=COLORS["white"])

    valid_annotations = tuple(x for x in annotations if x in metadata.columns)
    ann_width = max(0.65, 0.35 * max(1, len(valid_annotations)))
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[2.3, ann_width, max(9.0, n_features * 0.16), 1.5],
        wspace=0.03,
    )
    ax_den = fig.add_subplot(gs[0, 0])
    ax_ann = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])
    ax_count = fig.add_subplot(gs[0, 3])

    dd = dendrogram(
        z,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color=COLORS["slate_blue"],
        link_color_func=lambda _k: COLORS["slate_blue"],
        ax=ax_den,
    )
    row_order = list(dd["leaves"])
    ordered_samples = presence.index[row_order]

    prevalence = presence.sum(axis=0)
    ordered_features = sorted(
        presence.columns,
        key=lambda x: (-int(prevalence[x]), feature_modification(str(x)), str(x)),
    )
    ordered = presence.loc[ordered_samples, ordered_features]
    y_max = n_samples * 10

    ax_den.set_ylim(0, y_max)
    ax_den.set_xticks([])
    ax_den.set_yticks([])
    for spine in ax_den.spines.values():
        spine.set_visible(False)

    if valid_annotations:
        ann_meta = metadata.reindex(ordered_samples)
        ann_rgb = np.zeros((n_samples, len(valid_annotations), 3), dtype=float)
        ann_maps: dict[str, dict[str, str]] = {}
        for j, column in enumerate(valid_annotations):
            cmap = category_color_map(column, ann_meta[column])
            ann_maps[column] = cmap
            for i, value in enumerate(ann_meta[column]):
                ann_rgb[i, j, :] = to_rgb(cmap[normalise_category(value)])
        ax_ann.imshow(
            ann_rgb,
            aspect="auto",
            interpolation="nearest",
            origin="lower",
            extent=[0, len(valid_annotations), 0, y_max],
        )
        ax_ann.set_xticks(np.arange(len(valid_annotations)) + 0.5)
        ax_ann.set_xticklabels(valid_annotations, rotation=90, fontsize=7)
        ax_ann.set_yticks([])
        for spine in ax_ann.spines.values():
            spine.set_visible(False)
    else:
        ann_maps = {}
        ax_ann.axis("off")

    rgb = np.ones((n_samples, n_features, 3), dtype=float)
    for j, feature in enumerate(ordered_features):
        mod = feature_modification(str(feature))
        present_color = to_rgb(MOD_COLORS.get(mod, COLORS["deep_ocean"]))
        hits = ordered.iloc[:, j].to_numpy(dtype=bool)
        rgb[hits, j, :] = present_color

    ax_heat.imshow(
        rgb,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
        extent=[0, n_features, 0, y_max],
    )
    ax_heat.set_xticks(np.arange(n_features) + 0.5)
    ax_heat.set_xticklabels(ordered_features, rotation=90, fontsize=5.5)
    ax_heat.set_yticks(np.arange(n_samples) * 10 + 5)
    ax_heat.set_yticklabels(ordered_samples, fontsize=6.5)
    ax_heat.set_ylabel("Isolate")
    ax_heat.set_xlabel("Accessory motif family")
    ax_heat.set_title("Accessory methylome repertoire", fontsize=13, pad=12)
    ax_heat.tick_params(length=0)
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.5)

    counts = ordered.sum(axis=1).to_numpy(dtype=float)
    y_centres = np.arange(n_samples) * 10 + 5
    ax_count.barh(y_centres, counts, height=7.0, color=COLORS["muted_teal"], edgecolor="none")
    ax_count.set_ylim(0, y_max)
    ax_count.set_yticks([])
    ax_count.set_xlabel("n accessory\nfeatures", fontsize=7)
    ax_count.spines["top"].set_visible(False)
    ax_count.spines["right"].set_visible(False)
    ax_count.spines["left"].set_visible(False)

    legend_items = [
        Patch(facecolor=MOD_COLORS["m6A"], edgecolor="none", label="m6A presence"),
        Patch(facecolor=MOD_COLORS["m4C"], edgecolor="none", label="m4C presence"),
        Patch(facecolor=COLORS["white"], edgecolor=COLORS["charcoal"], linewidth=0.5, label="absence"),
    ]
    fig.legend(
        handles=legend_items,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.57, -0.01),
    )

    # Small metadata legends are useful; large category sets are left as strips
    # only so the figure does not become dominated by legend furniture.
    legend_y = -0.055
    for column, cmap in ann_maps.items():
        if len(cmap) > 10:
            continue
        handles = [Patch(facecolor=color, edgecolor="none", label=category) for category, color in cmap.items()]
        fig.legend(
            handles=handles,
            title=column,
            loc="lower center",
            ncol=min(5, len(handles)),
            frameon=False,
            bbox_to_anchor=(0.57, legend_y),
            fontsize=7,
            title_fontsize=7,
        )
        legend_y -= 0.055

    fig.subplots_adjust(bottom=max(0.14, -legend_y + 0.04))
    save_figure(fig, outdir / "figure_accessory_heatmap", formats)
    plt.close(fig)
    return [str(x) for x in ordered_samples]


def pcoa_variance_fraction(variance: pd.DataFrame, axis: str) -> float:
    row = variance.loc[variance["axis"].astype(str) == axis]
    if row.empty:
        return float("nan")
    return float(row.iloc[0]["variance_fraction"])


def plot_pcoa(
    pcoa: pd.DataFrame,
    variance: pd.DataFrame,
    metadata: pd.DataFrame,
    color_by: str,
    shape_by: str | None,
    label_points: bool,
    outdir: Path,
    formats: tuple[str, ...],
) -> None:
    if color_by not in metadata.columns:
        raise ValueError(f"metadata column not found for --color-by: {color_by}")
    if shape_by and shape_by not in metadata.columns:
        raise ValueError(f"metadata column not found for --shape-by: {shape_by}")
    if "PCoA1" not in pcoa.columns or "PCoA2" not in pcoa.columns:
        raise ValueError("pcoa.tsv must contain PCoA1 and PCoA2")

    samples = pcoa.index.intersection(metadata.index)
    data = pcoa.loc[samples]
    meta = metadata.loc[samples]

    color_map = category_color_map(color_by, meta[color_by])
    if shape_by:
        shape_categories = sorted({normalise_category(x) for x in meta[shape_by]})
        shape_map = {category: MARKERS[i % len(MARKERS)] for i, category in enumerate(shape_categories)}
    else:
        shape_map = {"all": "o"}

    fig, ax = plt.subplots(figsize=(8.2, 7.2), facecolor=COLORS["white"])
    for sample in samples:
        color_cat = normalise_category(meta.loc[sample, color_by])
        shape_cat = normalise_category(meta.loc[sample, shape_by]) if shape_by else "all"
        ax.scatter(
            float(data.loc[sample, "PCoA1"]),
            float(data.loc[sample, "PCoA2"]),
            s=62,
            c=[color_map[color_cat]],
            marker=shape_map[shape_cat],
            edgecolors=COLORS["charcoal"],
            linewidths=0.55,
            alpha=0.95,
            zorder=3,
        )
        if label_points:
            ax.annotate(
                str(sample),
                (float(data.loc[sample, "PCoA1"]), float(data.loc[sample, "PCoA2"])),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
                color=COLORS["charcoal"],
            )

    v1 = pcoa_variance_fraction(variance, "PCoA1")
    v2 = pcoa_variance_fraction(variance, "PCoA2")
    label1 = f"PCoA1 ({v1:.1%})" if np.isfinite(v1) else "PCoA1"
    label2 = f"PCoA2 ({v2:.1%})" if np.isfinite(v2) else "PCoA2"
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    ax.set_title("Accessory methylome PCoA", fontsize=13)
    ax.axhline(0, color=COLORS["light_grey"], lw=0.6, zorder=0)
    ax.axvline(0, color=COLORS["light_grey"], lw=0.6, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    color_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=color, markeredgecolor=COLORS["charcoal"],
               markeredgewidth=0.5, label=category)
        for category, color in color_map.items()
    ]
    color_legend = ax.legend(
        handles=color_handles,
        title=color_by,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=8,
        title_fontsize=8,
    )
    ax.add_artist(color_legend)

    if shape_by:
        shape_handles = [
            Line2D([0], [0], marker=marker, linestyle="", markersize=7,
                   markerfacecolor=COLORS["white"], markeredgecolor=COLORS["charcoal"],
                   label=category)
            for category, marker in shape_map.items()
        ]
        if len(shape_handles) <= 12:
            ax.legend(
                handles=shape_handles,
                title=shape_by,
                frameon=False,
                loc="lower left",
                bbox_to_anchor=(1.02, 0.0),
                fontsize=8,
                title_fontsize=8,
            )

    save_figure(fig, outdir / "figure_pcoa", formats)
    plt.close(fig)


def plot_prevalence(
    prevalence: pd.DataFrame,
    top_n: int,
    outdir: Path,
    formats: tuple[str, ...],
) -> None:
    work = prevalence.loc[prevalence["motif_family"].astype(str) != "RAATTY"].copy()
    work = work.sort_values(["isolates", "modification", "motif_family"], ascending=[False, True, True]).head(top_n)
    work = work.iloc[::-1]

    fig_height = max(5.5, 0.30 * len(work) + 1.5)
    fig, ax = plt.subplots(figsize=(8.5, fig_height), facecolor=COLORS["white"])
    bar_colors = [MOD_COLORS.get(str(x), COLORS["deep_ocean"]) for x in work["modification"]]
    labels = [f"{motif}  [{mod}]" for motif, mod in zip(work["motif_family"], work["modification"])]
    ax.barh(np.arange(len(work)), work["isolates"].to_numpy(), color=bar_colors, edgecolor="none")
    ax.set_yticks(np.arange(len(work)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Number of PASS isolates")
    ax.set_title(f"Top {len(work)} accessory motif families", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        handles=[
            Patch(facecolor=MOD_COLORS["m6A"], label="m6A"),
            Patch(facecolor=MOD_COLORS["m4C"], label="m4C"),
        ],
        frameon=False,
        loc="lower right",
    )
    save_figure(fig, outdir / "figure_accessory_prevalence", formats)
    plt.close(fig)


def plot_raatty_backbone(
    cohort_summary_path: Path,
    outdir: Path,
    formats: tuple[str, ...],
) -> None:
    df = read_tsv(cohort_summary_path)
    work = df.loc[
        (df["qc_status"] == "PASS")
        & (df["motif"] == "RAATTY")
        & (df["modification"] == "m6A")
    ].copy()
    if work.empty:
        return
    work["methylated_fraction"] = pd.to_numeric(work["methylated_fraction"], errors="coerce")
    work = work.sort_values(["methylated_fraction", "sample_id"])

    fig_height = max(7.0, 0.16 * len(work) + 1.2)
    fig, ax = plt.subplots(figsize=(7.2, fig_height), facecolor=COLORS["white"])
    y = np.arange(len(work))
    ax.hlines(y, 0.9, work["methylated_fraction"], color=COLORS["powder_blue"], lw=1.0)
    ax.scatter(work["methylated_fraction"], y, s=30, color=COLORS["deep_ocean"], edgecolors=COLORS["charcoal"], linewidths=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(work["sample_id"], fontsize=6.5)
    ax.set_xlim(min(0.9, float(work["methylated_fraction"].min()) - 0.01), 1.005)
    ax.set_xlabel("RAATTY m6A called targets / genomic targets")
    ax.set_title("Conserved RAATTY methylation backbone", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_figure(fig, outdir / "figure_raatty_backbone", formats)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot MOMENTO methylotype outputs using the Life Aquatic figure style")
    parser.add_argument("--methylotype-dir", required=True, type=Path, help="directory produced by momento methylotypes")
    parser.add_argument("--metadata", required=True, type=Path, help="TSV with one row per sample and sample_id column")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cohort-summary", type=Path, help="optional cohort.momento.tsv for the RAATTY backbone panel")
    parser.add_argument(
        "--annotations",
        default="host_group,provenance_class,species,clonal_complex",
        help="comma-separated metadata columns shown beside the heatmap",
    )
    parser.add_argument("--color-by", default="host_group", help="metadata column used for PCoA point colour")
    parser.add_argument("--shape-by", default="provenance_class", help="metadata column used for PCoA marker shape; pass empty string to disable")
    parser.add_argument("--label-points", action="store_true", help="label every PCoA point (usually best only for diagnostics)")
    parser.add_argument("--top-prevalence", type=int, default=20)
    parser.add_argument("--formats", default="png,pdf,svg", help="comma-separated output formats")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_matplotlib()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    formats = tuple(x.strip().lower() for x in args.formats.split(",") if x.strip())
    annotations = tuple(x.strip() for x in args.annotations.split(",") if x.strip())
    shape_by = args.shape_by.strip() if args.shape_by else None

    presence, pcoa, variance, linkage, prevalence, metadata = load_inputs(
        args.methylotype_dir,
        args.metadata,
    )

    ordered_samples = plot_heatmap(
        presence,
        linkage,
        metadata,
        annotations,
        args.output_dir,
        formats,
    )
    plot_pcoa(
        pcoa,
        variance,
        metadata,
        args.color_by,
        shape_by,
        args.label_points,
        args.output_dir,
        formats,
    )
    plot_prevalence(prevalence, args.top_prevalence, args.output_dir, formats)
    if args.cohort_summary:
        plot_raatty_backbone(args.cohort_summary, args.output_dir, formats)

    pd.DataFrame({"sample_id": ordered_samples}).to_csv(
        args.output_dir / "heatmap_sample_order.tsv",
        sep="\t",
        index=False,
    )
    print(
        f"FIGURES: {presence.shape[0]} samples; {presence.shape[1]} accessory features -> {args.output_dir}"
    )


if __name__ == "__main__":
    main()
