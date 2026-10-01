#!/usr/bin/env python3
"""Quantify lineage and host contributions to accessory methylome variation.

This script is intentionally an exploratory population-genetic analysis rather
than a claim of causal host effects. It answers two concrete questions:

1. Do more genetically distant isolates tend to have more different accessory
   methylome repertoires?
2. After representing genomic lineage with phylogenetic principal coordinates,
   how much additional accessory-methylome variation is associated with host?

The primary response is Jaccard distance across MOMENTO accessory motif
families. The genomic-lineage representation is derived from patristic
distances on the supplied core-genome tree. Small host groups are excluded from
inferential variance partitioning by default but remain in descriptive outputs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from momento.population import (
    jaccard_distance_matrix,
    model_r2,
    patristic_distance_matrix,
    pcoa,
    upper_triangle,
)
from momento.tree import prune_tree, read_newick, tip_names


def categorical_design(series: pd.Series) -> tuple[np.ndarray, list[str]]:
    values = series.fillna("unknown").astype(str)
    dummies = pd.get_dummies(values, prefix=series.name, dtype=float)
    return dummies.to_numpy(dtype=float), list(dummies.columns)


def with_intercept(*blocks: np.ndarray) -> np.ndarray:
    if not blocks:
        raise ValueError("at least one block is required")
    n = blocks[0].shape[0]
    cols = [np.ones((n, 1), dtype=float)]
    for block in blocks:
        if block.shape[0] != n:
            raise ValueError("design blocks have different row counts")
        if block.shape[1] > 0:
            cols.append(np.asarray(block, dtype=float))
    return np.column_stack(cols)


def mantel_permutation(
    a: np.ndarray,
    b: np.ndarray,
    permutations: int,
    seed: int,
) -> dict[str, float]:
    va = upper_triangle(a)
    vb = upper_triangle(b)
    pearson = float(pearsonr(va, vb).statistic)
    spearman = float(spearmanr(va, vb).statistic)
    rng = np.random.default_rng(seed)
    ge_pearson = 0
    ge_spearman = 0
    for _ in range(permutations):
        order = rng.permutation(b.shape[0])
        perm = b[np.ix_(order, order)]
        vp = upper_triangle(perm)
        rp = float(pearsonr(va, vp).statistic)
        rs = float(spearmanr(va, vp).statistic)
        if abs(rp) >= abs(pearson):
            ge_pearson += 1
        if abs(rs) >= abs(spearman):
            ge_spearman += 1
    return {
        "pearson_r": pearson,
        "pearson_p_perm": (ge_pearson + 1) / (permutations + 1),
        "spearman_rho": spearman,
        "spearman_p_perm": (ge_spearman + 1) / (permutations + 1),
    }


def choose_lineage_axes(
    eigenvalues: np.ndarray,
    target_variance: float,
    max_axes: int,
) -> int:
    if len(eigenvalues) == 0:
        raise ValueError("phylogenetic PCoA produced no positive axes")
    cumulative = np.cumsum(eigenvalues) / eigenvalues.sum()
    k = int(np.searchsorted(cumulative, target_variance) + 1)
    return max(1, min(k, max_axes, len(eigenvalues)))


def add_model_row(rows: list[dict], name: str, distance: np.ndarray, design: np.ndarray) -> None:
    r2, adjusted, p = model_r2(distance, design)
    rows.append({"model": name, "r2": r2, "adjusted_r2": adjusted, "predictor_df": p})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--methylotype-dir", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--host-column", default="host_group")
    parser.add_argument("--provenance-column", default="provenance_class")
    parser.add_argument("--min-host-group-size", type=int, default=3)
    parser.add_argument("--lineage-variance", type=float, default=0.80)
    parser.add_argument("--max-lineage-axes", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260917)
    args = parser.parse_args()

    if not 0 < args.lineage_variance <= 1:
        raise SystemExit("--lineage-variance must be in (0, 1]")
    if args.min_host_group_size < 2:
        raise SystemExit("--min-host-group-size must be >=2")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    presence = pd.read_csv(
        args.methylotype_dir / "accessory_presence_absence.tsv",
        sep="\t",
        index_col=0,
    )
    presence.index = presence.index.astype(str)
    presence = presence.apply(pd.to_numeric, errors="raise")

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    if "sample_id" not in metadata.columns:
        raise SystemExit("metadata must contain sample_id")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    if metadata["sample_id"].duplicated().any():
        raise SystemExit("metadata contains duplicate sample_id values")
    metadata = metadata.set_index("sample_id")

    missing_meta = sorted(set(presence.index) - set(metadata.index))
    if missing_meta:
        raise SystemExit("metadata missing samples: " + ", ".join(missing_meta))
    if args.host_column not in metadata.columns:
        raise SystemExit(f"metadata missing host column: {args.host_column}")

    raw_tree = read_newick(args.tree)
    missing_tree = sorted(set(presence.index) - set(tip_names(raw_tree)))
    if missing_tree:
        raise SystemExit("tree missing MOMENTO samples: " + ", ".join(missing_tree))
    tree = prune_tree(raw_tree, presence.index)

    samples = list(presence.index)
    _, phylo = patristic_distance_matrix(tree, samples)
    methylome = jaccard_distance_matrix(presence.to_numpy(dtype=bool))

    pd.DataFrame(phylo, index=samples, columns=samples).to_csv(
        args.output_dir / "core_phylogeny_patristic_distance.tsv", sep="\t"
    )
    pd.DataFrame(methylome, index=samples, columns=samples).to_csv(
        args.output_dir / "accessory_methylome_jaccard_distance.tsv", sep="\t"
    )

    mantel = mantel_permutation(
        phylo,
        methylome,
        permutations=args.permutations,
        seed=args.seed,
    )
    pd.DataFrame([mantel]).to_csv(
        args.output_dir / "lineage_methylome_distance_correlation.tsv",
        sep="\t",
        index=False,
    )

    host = metadata.loc[samples, args.host_column].fillna("unknown").astype(str)
    host_counts = host.value_counts().rename_axis("host_group").reset_index(name="n")
    host_counts["included_primary"] = host_counts["n"] >= args.min_host_group_size
    host_counts.loc[host_counts["host_group"].str.lower() == "unknown", "included_primary"] = False
    host_counts.to_csv(args.output_dir / "host_group_counts.tsv", sep="\t", index=False)

    included_groups = set(
        host_counts.loc[host_counts["included_primary"], "host_group"].astype(str)
    )
    include = host.isin(included_groups)
    primary_samples = [s for s in samples if bool(include.loc[s])]
    if len(primary_samples) < 8 or len(included_groups) < 2:
        raise SystemExit("too few samples/host groups remain for host-lineage analysis")

    idx = np.array([samples.index(s) for s in primary_samples], dtype=int)
    d_m = methylome[np.ix_(idx, idx)]
    d_p = phylo[np.ix_(idx, idx)]
    meta = metadata.loc[primary_samples].copy()

    lineage_coords, lineage_eig = pcoa(d_p)
    k = choose_lineage_axes(
        lineage_eig,
        target_variance=args.lineage_variance,
        max_axes=args.max_lineage_axes,
    )
    lineage = lineage_coords[:, :k]
    lineage_var = float(lineage_eig[:k].sum() / lineage_eig.sum())

    host_design, host_names = categorical_design(meta[args.host_column])

    rows: list[dict] = []
    add_model_row(rows, "host", d_m, with_intercept(host_design))
    add_model_row(rows, "lineage", d_m, with_intercept(lineage))
    add_model_row(rows, "host+lineage", d_m, with_intercept(lineage, host_design))

    provenance_used = False
    if args.provenance_column in meta.columns:
        prov = meta[args.provenance_column].fillna("unknown").astype(str)
        if prov.nunique() > 1:
            prov_design, prov_names = categorical_design(prov)
            provenance_used = True
            add_model_row(rows, "provenance", d_m, with_intercept(prov_design))
            add_model_row(rows, "lineage+provenance", d_m, with_intercept(lineage, prov_design))
            add_model_row(
                rows,
                "host+lineage+provenance",
                d_m,
                with_intercept(lineage, host_design, prov_design),
            )
        else:
            prov_names = []
    else:
        prov_names = []

    models = pd.DataFrame(rows)
    models.to_csv(args.output_dir / "distance_variance_models.tsv", sep="\t", index=False)

    lookup = models.set_index("model")
    components = [
        {
            "component": "unique_host_after_lineage",
            "raw_r2": lookup.loc["host+lineage", "r2"] - lookup.loc["lineage", "r2"],
            "adjusted_r2": lookup.loc["host+lineage", "adjusted_r2"] - lookup.loc["lineage", "adjusted_r2"],
        },
        {
            "component": "unique_lineage_after_host",
            "raw_r2": lookup.loc["host+lineage", "r2"] - lookup.loc["host", "r2"],
            "adjusted_r2": lookup.loc["host+lineage", "adjusted_r2"] - lookup.loc["host", "adjusted_r2"],
        },
        {
            "component": "shared_host_lineage",
            "raw_r2": lookup.loc["host", "r2"] + lookup.loc["lineage", "r2"] - lookup.loc["host+lineage", "r2"],
            "adjusted_r2": lookup.loc["host", "adjusted_r2"] + lookup.loc["lineage", "adjusted_r2"] - lookup.loc["host+lineage", "adjusted_r2"],
        },
    ]
    if provenance_used:
        components.append(
            {
                "component": "unique_host_after_lineage_and_provenance",
                "raw_r2": lookup.loc["host+lineage+provenance", "r2"] - lookup.loc["lineage+provenance", "r2"],
                "adjusted_r2": lookup.loc["host+lineage+provenance", "adjusted_r2"] - lookup.loc["lineage+provenance", "adjusted_r2"],
            }
        )
    pd.DataFrame(components).to_csv(
        args.output_dir / "variance_partition_components.tsv",
        sep="\t",
        index=False,
    )

    reset_meta = meta.reset_index()
    sample_audit = reset_meta[
        [
            c
            for c in [
                "sample_id",
                args.host_column,
                args.provenance_column,
                "host_detail",
                "st",
                "clonal_complex",
                "country",
                "year",
            ]
            if c in reset_meta.columns
        ]
    ]
    sample_audit.to_csv(args.output_dir / "analysis_samples.tsv", sep="\t", index=False)

    settings = pd.DataFrame(
        [
            ["n_all_samples", len(samples)],
            ["n_primary_samples", len(primary_samples)],
            ["primary_host_groups", ";".join(sorted(included_groups))],
            ["min_host_group_size", args.min_host_group_size],
            ["lineage_axes", k],
            ["lineage_positive_eigenvalue_variance_represented", lineage_var],
            ["lineage_variance_target", args.lineage_variance],
            ["max_lineage_axes", args.max_lineage_axes],
            ["mantel_permutations", args.permutations],
            ["random_seed", args.seed],
            ["host_design_columns", ";".join(host_names)],
            ["provenance_design_columns", ";".join(prov_names)],
        ],
        columns=["setting", "value"],
    )
    settings.to_csv(args.output_dir / "analysis_settings.tsv", sep="\t", index=False)

    print(f"LINEAGE/HOST ANALYSIS: {len(primary_samples)}/{len(samples)} samples in primary host comparison")
    print(f"  host groups: {', '.join(sorted(included_groups))}")
    print(f"  lineage axes: {k} ({lineage_var:.1%} of positive phylogenetic PCoA variance)")
    print(f"  distance correlation: Pearson r={mantel['pearson_r']:.3f}, permutation p={mantel['pearson_p_perm']:.4g}")
    print(f"  outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
