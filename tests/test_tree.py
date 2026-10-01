import pytest

from momento.tree import ladderize, parse_newick, prune_tree, tip_names, tree_layout


def test_parse_iqtree_style_newick():
    tree = parse_newick("((A:0.1,B:0.2)99:0.3,C:0.4);")
    assert sorted(tip_names(tree)) == ["A", "B", "C"]
    assert tree.children[0].name == "99"
    assert tree.children[0].length == pytest.approx(0.3)


def test_quoted_tip_labels_are_supported():
    tree = parse_newick("('sample one':0.1,B:0.2);")
    assert set(tip_names(tree)) == {"sample one", "B"}


def test_prune_collapses_unary_internal_nodes_and_preserves_distances():
    tree = parse_newick("(((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);")
    pruned = prune_tree(tree, {"A", "C", "D"})
    assert set(tip_names(pruned)) == {"A", "C", "D"}

    coords, order = tree_layout(pruned)
    tips = {node.name: node for node in _walk(pruned) if node.is_tip}
    # A keeps its terminal 0.1 plus the collapsed 0.3 branch.
    assert coords[id(tips["A"])][0] == pytest.approx(0.9)
    assert set(order) == {"A", "C", "D"}


def test_prune_requires_all_requested_samples():
    tree = parse_newick("(A:0.1,B:0.2);")
    with pytest.raises(ValueError, match="missing requested samples"):
        prune_tree(tree, {"A", "C"})


def test_ladderize_and_layout_return_all_tips_once():
    tree = parse_newick("((A:0.1,B:0.1):0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1);")
    ordered = ladderize(tree)
    _, tips = tree_layout(ordered)
    assert len(tips) == 5
    assert set(tips) == {"A", "B", "C", "D", "E"}


def test_duplicate_tip_labels_rejected():
    with pytest.raises(ValueError, match="Duplicate tree tip labels"):
        parse_newick("(A:0.1,A:0.2);")


def _walk(node):
    yield node
    for child in node.children:
        yield from _walk(child)
