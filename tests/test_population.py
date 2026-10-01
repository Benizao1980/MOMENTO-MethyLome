import numpy as np

from momento.population import jaccard_distance_matrix, model_r2, patristic_distance_matrix, pcoa
from momento.tree import parse_newick


def test_jaccard_distance_matrix():
    x = np.array([[1, 0, 1], [1, 1, 0], [1, 0, 1]], dtype=int)
    d = jaccard_distance_matrix(x)
    assert np.allclose(np.diag(d), 0)
    assert np.isclose(d[0, 1], 2 / 3)
    assert np.isclose(d[0, 2], 0)
    assert np.allclose(d, d.T)


def test_patristic_distance_matrix():
    tree = parse_newick("((A:1,B:2):3,C:4);")
    names, d = patristic_distance_matrix(tree, ["A", "B", "C"])
    assert names == ["A", "B", "C"]
    assert np.isclose(d[0, 1], 3.0)
    assert np.isclose(d[0, 2], 8.0)
    assert np.isclose(d[1, 2], 9.0)


def test_pcoa_and_model_r2_are_finite():
    d = np.array(
        [
            [0.0, 1.0, 2.0, 2.0],
            [1.0, 0.0, 1.0, 2.0],
            [2.0, 1.0, 0.0, 1.0],
            [2.0, 2.0, 1.0, 0.0],
        ]
    )
    coords, eig = pcoa(d)
    assert coords.shape[0] == 4
    assert len(eig) >= 1
    design = np.column_stack([np.ones(4), [0, 0, 1, 1]])
    r2, adj, p = model_r2(d, design)
    assert np.isfinite(r2)
    assert np.isfinite(adj)
    assert p == 1
