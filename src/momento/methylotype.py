"""Motif-family cleanup and exploratory accessory methylotype analysis.

The routines in this module are deliberately conservative. Reverse-complement
representations are collapsed because they describe the same recognition
specificity. IUPAC motif nesting is *flagged* rather than automatically removed:
a narrower motif can represent a caller artefact, but it can also reflect a real
additional methyltransferase specificity.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from .cluster import jaccard_distance
from .io.pacbio import canonicalise_motif, reverse_complement_iupac

IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def reverse_complement_motif(motif: str) -> str:
    """Reverse-complement an IUPAC motif, including slash-separated partners."""
    motif = canonicalise_motif(motif)
    parts = motif.split("/")
    return "/".join(reverse_complement_iupac(part) for part in reversed(parts))


def canonical_motif_family(motif: str) -> str:
    """Return a stable family label shared by a motif and its reverse complement."""
    motif = canonicalise_motif(motif)
    rc = reverse_complement_motif(motif)
    return min(motif, rc)


def motif_orientations(motif: str) -> tuple[str, ...]:
    """Return unique single-component motif orientations."""
    motif = canonicalise_motif(motif)
    if "/" in motif:
        return ()
    rc = reverse_complement_iupac(motif)
    return tuple(dict.fromkeys((motif, rc)))


def nested_alignment(parent: str, child: str) -> dict | None:
    """Return an IUPAC containment alignment when ``child`` is nested in ``parent``.

    ``parent`` is the broader recognition language. Every concrete sequence
    matching ``child`` must also match ``parent`` at the returned window. Both
    strand orientations are tested. Compound slash motifs are not compared.
    ``child_start_1based`` reports where the parent window starts within the
    oriented child string.
    """
    parent_orientations = motif_orientations(parent)
    child_orientations = motif_orientations(child)
    if not parent_orientations or not child_orientations:
        return None

    for p in parent_orientations:
        for c in child_orientations:
            if len(p) > len(c):
                continue
            for offset in range(len(c) - len(p) + 1):
                if all(
                    IUPAC_BASES[c[offset + i]].issubset(IUPAC_BASES[p[i]])
                    for i in range(len(p))
                ):
                    return {
                        "parent_orientation": p,
                        "child_orientation": c,
                        "child_start_1based": offset + 1,
                    }
    return None


def collapse_motif_families(
    summary: pd.DataFrame,
    *,
    qc_statuses: Sequence[str] = ("PASS",),
    modifications: Sequence[str] = ("m6A", "m4C"),
    min_called_sites: int = 1,
) -> pd.DataFrame:
    """Collapse reverse-complement motif labels within each sample.

    Counts are retained for diagnostics, but family-level occupancy is not
    calculated here because reverse-complement summary rows can share genomic
    opportunity denominators. Presence/absence is the intended first-pass use.
    """
    required = {"sample_id", "motif", "modification", "called_sites", "qc_status"}
    missing = sorted(required - set(summary.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))

    work = summary.loc[
        summary["qc_status"].isin(qc_statuses)
        & summary["modification"].isin(modifications)
    ].copy()
    work["called_sites"] = pd.to_numeric(work["called_sites"], errors="coerce").fillna(0)
    work = work.loc[work["called_sites"] >= min_called_sites].copy()
    work["motif_family"] = work["motif"].map(canonical_motif_family)

    rows = []
    for (sample, modification, family), group in work.groupby(
        ["sample_id", "modification", "motif_family"], sort=True
    ):
        motifs = sorted(set(group["motif"].astype(str)))
        rows.append({
            "sample_id": sample,
            "modification": modification,
            "motif_family": family,
            "raw_motifs": ";".join(motifs),
            "n_raw_motifs": len(motifs),
            "max_called_sites": int(group["called_sites"].max()),
            "sum_called_sites": int(group["called_sites"].sum()),
        })
    return pd.DataFrame(rows)


def build_family_presence(
    family_long: pd.DataFrame,
    *,
    exclude_families: Iterable[str] = ("RAATTY",),
) -> pd.DataFrame:
    """Build a sample x modification|motif-family binary matrix.

    Samples carrying only excluded/core families are retained as all-zero rows,
    which is important for unbiased cohort-level distance calculations.
    """
    all_samples = sorted(set(family_long["sample_id"].astype(str)))
    exclude = {canonical_motif_family(x) for x in exclude_families}
    work = family_long.loc[~family_long["motif_family"].isin(exclude)].copy()
    work["feature"] = work["modification"].astype(str) + "|" + work["motif_family"].astype(str)
    work["present"] = 1
    matrix = work.pivot_table(
        index="sample_id",
        columns="feature",
        values="present",
        aggfunc="max",
        fill_value=0,
    )
    matrix = matrix.reindex(all_samples, fill_value=0)
    return matrix.sort_index().sort_index(axis=1).astype(int)


def family_prevalence(family_long: pd.DataFrame) -> pd.DataFrame:
    """Count isolates carrying each modification/motif family."""
    return (
        family_long.groupby(["modification", "motif_family"], as_index=False)
        .agg(
            isolates=("sample_id", "nunique"),
            raw_motif_labels=("raw_motifs", lambda x: ";".join(sorted(set(x)))),
        )
        .sort_values(["isolates", "modification", "motif_family"], ascending=[False, True, True])
        .reset_index(drop=True)
    )


def find_nested_motif_pairs(
    family_long: pd.DataFrame,
    *,
    excluded_families: Iterable[str] = ("RAATTY",),
) -> pd.DataFrame:
    """Flag IUPAC motif-family containment and cohort co-occurrence.

    Nesting is reported, not removed. Co-occurrence fractions help distinguish
    likely redundant caller representations from potentially independent
    specificities, but they are not treated as proof either way.
    """
    excluded = {canonical_motif_family(x) for x in excluded_families}
    rows: list[dict] = []

    for modification, group in family_long.groupby("modification"):
        sample_sets = {
            family: set(sub["sample_id"].astype(str))
            for family, sub in group.groupby("motif_family")
        }
        families = sorted(sample_sets)

        for a, b in combinations(families, 2):
            ab = nested_alignment(a, b)
            ba = nested_alignment(b, a)
            if ab is not None and ba is not None:
                # Equivalent recognition languages represented by different labels.
                continue
            if ab is not None:
                parent, child, alignment = a, b, ab
            elif ba is not None:
                parent, child, alignment = b, a, ba
            else:
                continue

            parent_samples = sample_sets[parent]
            child_samples = sample_sets[child]
            overlap = parent_samples & child_samples
            rows.append({
                "modification": modification,
                "parent_family": parent,
                "child_family": child,
                "parent_orientation": alignment["parent_orientation"],
                "child_orientation": alignment["child_orientation"],
                "child_start_1based": alignment["child_start_1based"],
                "parent_isolates": len(parent_samples),
                "child_isolates": len(child_samples),
                "cooccurring_isolates": len(overlap),
                "child_with_parent_fraction": (
                    len(overlap) / len(child_samples) if child_samples else np.nan
                ),
                "parent_with_child_fraction": (
                    len(overlap) / len(parent_samples) if parent_samples else np.nan
                ),
                "parent_excluded_from_accessory": parent in excluded,
                "child_excluded_from_accessory": child in excluded,
            })

    columns = [
        "modification", "parent_family", "child_family", "parent_orientation",
        "child_orientation", "child_start_1based", "parent_isolates",
        "child_isolates", "cooccurring_isolates", "child_with_parent_fraction",
        "parent_with_child_fraction", "parent_excluded_from_accessory",
        "child_excluded_from_accessory",
    ]
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["child_with_parent_fraction", "child_isolates", "modification", "parent_family", "child_family"],
        ascending=[False, False, True, True, True],
    ).reset_index(drop=True)


def classical_pcoa(distance: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classical PCoA of a symmetric distance matrix."""
    if distance.shape[0] != distance.shape[1]:
        raise ValueError("distance matrix must be square")
    if list(distance.index) != list(distance.columns):
        raise ValueError("distance matrix rows and columns must have the same sample order")

    n = len(distance)
    if n == 0:
        return pd.DataFrame(), pd.DataFrame(columns=["axis", "eigenvalue", "variance_fraction"])

    d2 = np.square(distance.to_numpy(dtype=float))
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ d2 @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = np.clip(eigenvalues, 0, None)
    total_positive = positive.sum()
    k = min(n_components, n)
    coords = eigenvectors[:, :k] * np.sqrt(positive[:k])
    coord_df = pd.DataFrame(
        coords,
        index=distance.index,
        columns=[f"PCoA{i + 1}" for i in range(k)],
    )
    coord_df.index.name = "sample_id"

    variance = pd.DataFrame({
        "axis": [f"PCoA{i + 1}" for i in range(k)],
        "eigenvalue": eigenvalues[:k],
        "variance_fraction": (
            positive[:k] / total_positive if total_positive > 0 else np.zeros(k)
        ),
    })
    return coord_df, variance


def average_linkage_order(distance: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Average-linkage hierarchical clustering from a precomputed distance matrix."""
    n = len(distance)
    if n == 0:
        return pd.DataFrame(columns=["sample_id", "cluster_order"]), pd.DataFrame(
            columns=["left", "right", "distance", "n_members"]
        )
    if n == 1:
        return pd.DataFrame({"sample_id": distance.index, "cluster_order": [1]}), pd.DataFrame(
            columns=["left", "right", "distance", "n_members"]
        )

    condensed = squareform(distance.to_numpy(dtype=float), checks=False)
    z = linkage(condensed, method="average")
    leaves = leaves_list(z)
    order_df = pd.DataFrame({
        "sample_id": [distance.index[i] for i in leaves],
        "cluster_order": np.arange(1, n + 1),
    })
    linkage_df = pd.DataFrame(z, columns=["left", "right", "distance", "n_members"])
    return order_df, linkage_df


def run_methylotype_analysis(
    summary: pd.DataFrame,
    output_dir: str | Path,
    *,
    qc_statuses: Sequence[str] = ("PASS",),
    modifications: Sequence[str] = ("m6A", "m4C"),
    exclude_families: Sequence[str] = ("RAATTY",),
    min_called_sites: int = 1,
) -> dict[str, object]:
    """Run the first-pass binary accessory methylotype workflow and write TSV outputs."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    family_long = collapse_motif_families(
        summary,
        qc_statuses=qc_statuses,
        modifications=modifications,
        min_called_sites=min_called_sites,
    )
    prevalence = family_prevalence(family_long)
    nested = find_nested_motif_pairs(family_long, excluded_families=exclude_families)
    presence = build_family_presence(family_long, exclude_families=exclude_families)
    if presence.empty or presence.shape[1] == 0:
        raise ValueError("No accessory motif-family features remain after filtering")

    distance = jaccard_distance(presence)
    pcoa, variance = classical_pcoa(distance, n_components=2)
    cluster_order, linkage_df = average_linkage_order(distance)

    family_long.to_csv(outdir / "motif_families_long.tsv", sep="\t", index=False)
    prevalence.to_csv(outdir / "motif_family_prevalence.tsv", sep="\t", index=False)
    nested.to_csv(outdir / "nested_motif_pairs.tsv", sep="\t", index=False)
    presence.to_csv(outdir / "accessory_presence_absence.tsv", sep="\t")
    distance.to_csv(outdir / "jaccard_distance.tsv", sep="\t")
    pcoa.to_csv(outdir / "pcoa.tsv", sep="\t")
    variance.to_csv(outdir / "pcoa_variance.tsv", sep="\t", index=False)
    cluster_order.to_csv(outdir / "cluster_order.tsv", sep="\t", index=False)
    linkage_df.to_csv(outdir / "linkage.tsv", sep="\t", index=False)

    return {
        "family_long": family_long,
        "prevalence": prevalence,
        "nested": nested,
        "presence": presence,
        "distance": distance,
        "pcoa": pcoa,
        "variance": variance,
        "cluster_order": cluster_order,
        "linkage": linkage_df,
    }
