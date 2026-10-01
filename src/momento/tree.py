"""Lightweight Newick parsing and tree-layout helpers for MOMENTO figures.

The parser is intentionally small and dependency-free. It supports the simple
Newick written by common bacterial phylogeny tools such as IQ-TREE: tip names,
internal labels/supports and branch lengths. The plotting workflow prunes trees
to the QC-clean MOMENTO sample set before rendering.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class TreeNode:
    name: str | None = None
    length: float = 0.0
    children: list["TreeNode"] = field(default_factory=list)

    @property
    def is_tip(self) -> bool:
        return not self.children


def _skip_ws(text: str, pos: int) -> int:
    while pos < len(text) and text[pos].isspace():
        pos += 1
    return pos


def _parse_label(text: str, pos: int) -> tuple[str | None, int]:
    pos = _skip_ws(text, pos)
    if pos >= len(text) or text[pos] in ",():;":
        return None, pos
    if text[pos] == "'":
        pos += 1
        chars: list[str] = []
        while pos < len(text):
            if text[pos] == "'":
                if pos + 1 < len(text) and text[pos + 1] == "'":
                    chars.append("'")
                    pos += 2
                    continue
                return "".join(chars), pos + 1
            chars.append(text[pos])
            pos += 1
        raise ValueError("Unterminated quoted Newick label")

    start = pos
    while pos < len(text) and text[pos] not in ",():;":
        pos += 1
    label = text[start:pos].strip()
    return (label or None), pos


def _parse_length(text: str, pos: int) -> tuple[float, int]:
    pos = _skip_ws(text, pos)
    if pos >= len(text) or text[pos] != ":":
        return 0.0, pos
    pos += 1
    start = pos
    while pos < len(text) and text[pos] not in ",();":
        pos += 1
    token = text[start:pos].strip()
    if not token:
        raise ValueError("Empty Newick branch length")
    try:
        return float(token), pos
    except ValueError as exc:
        raise ValueError(f"Invalid Newick branch length: {token!r}") from exc


def _parse_subtree(text: str, pos: int) -> tuple[TreeNode, int]:
    pos = _skip_ws(text, pos)
    if pos >= len(text):
        raise ValueError("Unexpected end of Newick")

    if text[pos] == "(":
        pos += 1
        children: list[TreeNode] = []
        while True:
            child, pos = _parse_subtree(text, pos)
            children.append(child)
            pos = _skip_ws(text, pos)
            if pos >= len(text):
                raise ValueError("Unterminated Newick internal node")
            if text[pos] == ",":
                pos += 1
                continue
            if text[pos] == ")":
                pos += 1
                break
            raise ValueError(f"Unexpected Newick character {text[pos]!r} at {pos}")
        name, pos = _parse_label(text, pos)
        length, pos = _parse_length(text, pos)
        return TreeNode(name=name, length=length, children=children), pos

    name, pos = _parse_label(text, pos)
    if name is None:
        raise ValueError(f"Missing Newick tip label at {pos}")
    length, pos = _parse_length(text, pos)
    return TreeNode(name=name, length=length), pos


def parse_newick(text: str) -> TreeNode:
    """Parse one Newick tree into :class:`TreeNode`."""
    text = text.strip()
    if not text:
        raise ValueError("Empty Newick input")
    root, pos = _parse_subtree(text, 0)
    pos = _skip_ws(text, pos)
    if pos < len(text) and text[pos] == ";":
        pos += 1
    pos = _skip_ws(text, pos)
    if pos != len(text):
        raise ValueError(f"Unexpected trailing Newick content: {text[pos:pos+40]!r}")
    validate_tree(root)
    return root


def read_newick(path: str | Path) -> TreeNode:
    return parse_newick(Path(path).read_text())


def iter_tips(node: TreeNode):
    if node.is_tip:
        yield node
        return
    for child in node.children:
        yield from iter_tips(child)


def tip_names(node: TreeNode) -> list[str]:
    return [str(tip.name) for tip in iter_tips(node)]


def validate_tree(node: TreeNode) -> None:
    tips = tip_names(node)
    if not tips:
        raise ValueError("Tree contains no tips")
    missing = [x for x in tips if not x or x == "None"]
    if missing:
        raise ValueError("Tree contains unnamed tips")
    dup = sorted({x for x in tips if tips.count(x) > 1})
    if dup:
        raise ValueError("Duplicate tree tip labels: " + ", ".join(dup))


def clone_tree(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        length=float(node.length),
        children=[clone_tree(child) for child in node.children],
    )


def prune_tree(node: TreeNode, keep: Iterable[str]) -> TreeNode:
    """Return a copy pruned to ``keep`` tip names.

    Unary internal nodes created by pruning are collapsed and their branch
    lengths are added to the retained child. The returned root has length 0.
    """
    keep_set = {str(x) for x in keep}

    def rec(current: TreeNode) -> TreeNode | None:
        if current.is_tip:
            if str(current.name) in keep_set:
                return TreeNode(name=str(current.name), length=float(current.length))
            return None
        children = [x for child in current.children if (x := rec(child)) is not None]
        if not children:
            return None
        if len(children) == 1:
            child = children[0]
            child.length += float(current.length)
            return child
        return TreeNode(name=current.name, length=float(current.length), children=children)

    result = rec(node)
    if result is None:
        raise ValueError("No requested samples remain after pruning tree")
    result.length = 0.0
    validate_tree(result)
    observed = set(tip_names(result))
    missing = sorted(keep_set - observed)
    if missing:
        raise ValueError("Tree is missing requested samples: " + ", ".join(missing))
    return result


def descendant_tip_count(node: TreeNode) -> int:
    if node.is_tip:
        return 1
    return sum(descendant_tip_count(child) for child in node.children)


def ladderize(node: TreeNode, *, reverse: bool = True) -> TreeNode:
    """Return a cloned tree with children sorted by descendant-tip count."""
    result = clone_tree(node)

    def rec(current: TreeNode) -> None:
        for child in current.children:
            rec(child)
        current.children.sort(
            key=lambda child: (descendant_tip_count(child), min(tip_names(child))),
            reverse=reverse,
        )

    rec(result)
    return result


def tree_layout(node: TreeNode) -> tuple[dict[int, tuple[float, float]], list[str]]:
    """Return node coordinates keyed by object id plus tip order.

    X is cumulative branch length when any positive branch lengths exist;
    otherwise it is cladogram depth. Y positions are consecutive tip rows.
    """
    tips = tip_names(node)
    use_lengths = any(float(n.length) > 0 for n in _walk(node) if n is not node)
    coords: dict[int, tuple[float, float]] = {}
    ordered_tips: list[str] = []

    def rec(current: TreeNode, x: float, depth: int) -> float:
        here_x = x if use_lengths else float(depth)
        if current.is_tip:
            y = float(len(ordered_tips))
            ordered_tips.append(str(current.name))
            coords[id(current)] = (here_x, y)
            return y
        child_y: list[float] = []
        for child in current.children:
            child_x = x + max(0.0, float(child.length)) if use_lengths else x
            child_y.append(rec(child, child_x, depth + 1))
        y = sum(child_y) / len(child_y)
        coords[id(current)] = (here_x, y)
        return y

    rec(node, 0.0, 0)
    if len(ordered_tips) != len(tips):
        raise RuntimeError("Tree layout tip count mismatch")
    return coords, ordered_tips


def _walk(node: TreeNode):
    yield node
    for child in node.children:
        yield from _walk(child)
