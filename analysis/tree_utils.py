"""Newick utilities (standard library only)."""

import re


def read_tree_file(path):
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"empty tree file: {path}")
    if re.fullmatch(r"\d+\s+\d+", lines[0].strip()):
        lines = lines[1:]
    m = re.search(r"[^;]+;", re.sub(r"\s+", "", "".join(lines)))
    if not m:
        raise ValueError(f"no Newick tree in {path}")
    return m.group(0)


class Tree:
    def __init__(self):
        self.nodes = {}
        self.root = None


def parse_newick(txt):
    txt = re.sub(r"\s+", "", txt)
    tree = Tree()
    pos = [0]
    n = len(txt)

    def parse_clade():
        node_id = len(tree.nodes) + 1
        node = {"children": [], "name": None, "length": None, "calib": None}
        tree.nodes[node_id] = node
        if pos[0] < n and txt[pos[0]] == "(":
            pos[0] += 1
            while True:
                child = parse_clade()
                node["children"].append(child)
                ch = txt[pos[0]] if pos[0] < n else ""
                if ch == ",":
                    pos[0] += 1
                elif ch == ")":
                    pos[0] += 1
                    break
                else:
                    raise ValueError(f"malformed newick at {pos[0]}")
        else:
            m = re.match(r"'[^']*'|[^(),:;]+", txt[pos[0]:])
            if not m:
                raise ValueError(f"malformed newick at {pos[0]}")
            node["name"] = m.group(0).replace("'", "")
            pos[0] += m.end()
        if pos[0] < n and txt[pos[0]] == "'":
            close = txt.find("'", pos[0] + 1)
            if close > 0:
                node["calib"] = txt[pos[0] + 1:close]
                pos[0] = close + 1
        if pos[0] < n and txt[pos[0]] == ":":
            pos[0] += 1
            m = re.match(r"[0-9.eE+-]+", txt[pos[0]:])
            if m:
                node["length"] = float(m.group(0))
                pos[0] += m.end()
        return node_id

    tree.root = parse_clade()
    return tree


def tip_names(tree):
    return [n["name"] for n in tree.nodes.values() if not n["children"]]


def descendants(tree, node_id):
    node = tree.nodes[node_id]
    if not node["children"]:
        return [node["name"]]
    out = []
    for child in node["children"]:
        out.extend(descendants(tree, child))
    return sorted(out)


def paml_node_ids(tree, tip_order):
    ids = {}
    internal = [0]

    def assign(node_id):
        node = tree.nodes[node_id]
        if not node["children"]:
            ids[node_id] = tip_order.index(node["name"]) + 1
        else:
            internal[0] += 1
            ids[node_id] = len(tip_order) + internal[0]
        for child in node["children"]:
            assign(child)

    assign(tree.root)
    return ids


def preorder_edges(tree):
    edges = []

    def visit(node_id):
        for child in tree.nodes[node_id]["children"]:
            edges.append((node_id, child))
            visit(child)

    visit(tree.root)
    return edges


def serialize_newick(tree, lengths=None, keep_calib=False):
    def rec(node_id):
        node = tree.nodes[node_id]
        if node["children"]:
            text = "(" + ",".join(rec(c) for c in node["children"]) + ")"
            if keep_calib and node["calib"]:
                text += "'" + node["calib"] + "'"
        else:
            text = node["name"]
        if lengths is not None and node_id != tree.root and lengths.get(node_id) is not None:
            text += ":%s" % format(lengths[node_id], ".12g")
        return text

    return rec(tree.root) + ";"
