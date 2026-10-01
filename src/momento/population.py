"""Population-level helpers for lineage-aware methylome analysis.

The functions in this module are intentionally small and auditable. They avoid
hiding lineage correction inside a black-box model so that the Campylobacter
worked example can expose each intermediate quantity (distance matrices,
principal coordinates and variance components).
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from momento.tree import TreeNode, iter_tips


def jaccard_distance_matrix(binary: np.ndarray) -> np.ndarray:
    """Return pairwise Jaccard distances for a binary presence/absence matrix.

    Rows are isolates and columns are features. Two all-zero rows are assigned
    distance 0 because their observed accessory repertoires are identical.
    """
    x = np.asarray(binary, dtype=bool)
    if x.ndim != 2:
        raise ValueError("binary matrix must be two-dimensional")
    n = x.shape[0]
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            union = np.logical_or(x[i], x[j]).sum()
            if union == 0:
                d = 0.0
            else:
                intersection = np.logical_and(x[i], x[j]).sum()
                d = 1.0 - (intersection / union)
            out[i, j] = out[j, i] = float(d)
    return out


def patristic_distance_matrix(
    root: TreeNode,
    tip_order: Sequence[str] | None = None,
) -> tuple[list[str], np.ndarray]:
    """Return tip names and pairwise patristic distances from a Newick tree.

    Patristic distance is the sum of branch lengths along the path joining two
    tips. It is independent of the arbitrary display root used for an unrooted
    IQ-TREE Newick file.
    """
    tips = list(iter_tips(root))
    by_name = {str(t.name): t for t in tips}
    if len(by_name) != len(tips):
        raise ValueError("tree contains duplicate tip names")

    if tip_order is None:
        names = [str(t.name) for t in tips]
    else:
        names = [str(x) for x in tip_order]
        missing = [x for x in names if x not in by_name]
        if missing:
            raise ValueError("tree missing requested tips: " + ", ".join(missing))

    root_distance: dict[int, float] = {}
    tip_paths: dict[str, list[int]] = {}

    def walk(node: TreeNode, distance: float, ancestors: list[int]) -> None:
        node_id = id(node)
        root_distance[node_id] = float(distance)
        here = ancestors + [node_id]
        if node.is_tip:
            tip_paths[str(node.name)] = here
            return
        for child in node.children:
            walk(child, distance + max(0.0, float(child.length)), here)

    walk(root, 0.0, [])

    n = len(names)
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        p_i = tip_paths[names[i]]
        set_i = set(p_i)
        d_i = root_distance[p_i[-1]]
        for j in range(i + 1, n):
            p_j = tip_paths[names[j]]
            d_j = root_distance[p_j[-1]]
            shared = [node_id for node_id in p_j if node_id in set_i]
            if not shared:
                raise RuntimeError("tree tips have no shared root")
            mrca = max(shared, key=lambda node_id: root_distance[node_id])
            d = d_i + d_j - 2.0 * root_distance[mrca]
            out[i, j] = out[j, i] = max(0.0, float(d))
    return names, out


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Classical principal coordinates analysis of a distance matrix.

    Returns coordinates and positive eigenvalues only. Small negative
    eigenvalues can occur for non-Euclidean distances such as Jaccard; they are
    not silently converted into positive axes.
    """
    d = np.asarray(distance, dtype=float)
    if d.ndim != 2 or d.shape[0] != d.shape[1]:
        raise ValueError("distance matrix must be square")
    if not np.allclose(d, d.T, atol=1e-10):
        raise ValueError("distance matrix must be symmetric")
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    tol = max(1e-12, abs(eigenvalues[0]) * 1e-10) if len(eigenvalues) else 1e-12
    keep = eigenvalues > tol
    positive = eigenvalues[keep]
    coordinates = eigenvectors[:, keep] * np.sqrt(positive)
    return coordinates, positive


def gower_center(distance: np.ndarray) -> np.ndarray:
    """Return the Gower-centred matrix used by distance-based R2."""
    d = np.asarray(distance, dtype=float)
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ (d ** 2) @ j


def model_r2(distance: np.ndarray, design: np.ndarray) -> tuple[float, float, int]:
    """Return raw R2, adjusted R2 and predictor degrees of freedom.

    `design` may contain an intercept or redundant dummy columns. Matrix rank,
    rather than nominal column count, is used for degrees of freedom.
    """
    g = gower_center(distance)
    x = np.asarray(design, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    if x.shape[0] != g.shape[0]:
        raise ValueError("design rows must match distance matrix")
    if not np.all(np.isfinite(x)):
        raise ValueError("design contains non-finite values")

    h = x @ np.linalg.pinv(x)
    total = float(np.trace(g))
    fitted = float(np.trace(h @ g))
    r2 = fitted / total if total > 0 else 0.0

    rank = int(np.linalg.matrix_rank(x))
    p = max(0, rank - 1)
    n = x.shape[0]
    if n - p - 1 <= 0:
        adjusted = float("nan")
    else:
        adjusted = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)
    return float(r2), float(adjusted), p


def upper_triangle(distance: np.ndarray) -> np.ndarray:
    """Return the upper-triangle distances as a 1D vector."""
    d = np.asarray(distance, dtype=float)
    idx = np.triu_indices(d.shape[0], k=1)
    return d[idx]
