"""Convert a GFmixF Hessian into a non-reversible MCMCtree Hessian and control file."""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
sys.path.insert(0, HERE)
from tree_utils import (descendants, preorder_edges, read_tree_file, parse_newick,
                        serialize_newick)
from gfmix_common import read_node_clades


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--ctl-template", required=True)
    ap.add_argument("--ctl-out", required=True)
    ap.add_argument("--hessian-out", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--base-prefix", required=True)
    a = ap.parse_args()
    theta = [float(x) for x in open(os.path.join(a.fit_dir, "theta.txt")).read().split()]
    grad = [float(x) for x in open(os.path.join(a.fit_dir, "gradient.txt")).read().split()]
    p = len(theta)
    hess = [[float(x) for x in ln.split()]
            for ln in open(os.path.join(a.fit_dir, "hessian.txt")) if ln.strip()]
    if len(hess) != p or len(grad) != p:
        raise SystemExit("inconsistent Hessian dimensions")
    nc = read_node_clades(os.path.join(a.fit_dir, "node_clades.txt"))
    clade_key = {}
    for node in range(p):
        clade_key[node] = ",".join(sorted(nc[node]["clade"]))
    tree = parse_newick(read_tree_file(a.tree))
    tips = [re.sub(r"^>\s*", "", ln) for ln in open(a.aln, encoding="utf-8") if ln.startswith(">")]
    tips = [re.sub(r"\s.*$", "", t) for t in tips]
    edges = preorder_edges(tree)
    edge_clade = [",".join(sorted(descendants(tree, ch))) for _, ch in edges]
    perm = []
    for key in edge_clade:
        hits = [node for node, val in clade_key.items() if val == key]
        if len(hits) != 1:
            raise SystemExit("cannot map dating-tree edge to GFmix branch: " + key)
        perm.append(hits[0])
    if len(set(perm)) != p or p != len(edges):
        raise SystemExit("edge mapping is not one-to-one")
    bl = [theta[i] for i in perm]
    g = [grad[i] for i in perm]
    ho = [[hess[i][j] for j in perm] for i in perm]
    lengths = {ch: bl[k] for k, (_, ch) in enumerate(edges)}
    tree_line = serialize_newick(tree, lengths=lengths)
    with open(a.hessian_out, "w", encoding="utf-8") as fh:
        fh.write("\n%s\n\n%s\n\n%s\n\n%s\n\nHessian\n\n" % (
            len(tips), tree_line,
            " ".join(format(x, ".17e") for x in bl),
            " ".join(format(x, ".17e") for x in g)))
        for row in ho:
            fh.write(" ".join(format(x, ".17e") for x in row) + "\n")
    ctl = open(a.ctl_template, encoding="utf-8").read()
    ctl = ctl.replace("@PREFIX@", a.out_prefix)
    ctl = ctl.replace("@SEQFILE@", a.base_prefix + ".dummy.phy")
    ctl = ctl.replace("@TREEFILE@", a.base_prefix + ".rooted.nwk")
    ctl = ctl.replace("@HESSIANFILE@", os.path.basename(a.hessian_out))
    open(a.ctl_out, "w", encoding="utf-8").write(ctl)
    print("GFmixF Hessian: %s (%d branches)" % (a.hessian_out, p))


if __name__ == "__main__":
    main()
