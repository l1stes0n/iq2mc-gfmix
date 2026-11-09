"""Build GFmixF inputs from the workflow FASTA alignment, tree and IQ-TREE report."""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
sys.path.insert(0, HERE)
from tree_utils import descendants, parse_newick, read_tree_file, serialize_newick, tip_names
from gfmix_common import AA, read_fasta


def seed_lengths(tree, tree_short, seed_file):
    ut = parse_newick(read_tree_file(seed_file))
    all_tips = sorted(tip_names(ut))
    ut_key = {}
    for nid in ut.nodes:
        for ch in ut.nodes[nid]["children"]:
            ut_key[",".join(descendants(ut, ch))] = ut.nodes[ch]["length"]
    for nid in tree.nodes:
        if nid == tree.root:
            continue
        clade = ",".join(descendants(tree, nid))
        comp = ",".join(sorted(set(all_tips) - set(descendants(tree, nid))))
        length = ut_key.get(clade, ut_key.get(comp))
        if length is None:
            continue
        root_children = tree.nodes[tree.root]["children"]
        tree_short.nodes[nid]["length"] = length / (2 if nid in root_children else 1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--iqtree", required=True)
    ap.add_argument("--gfmixf-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--classes", type=int, default=20)
    ap.add_argument("--groups", default="GARP,FYMINK")
    ap.add_argument("--seed-tree", default="")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    names, seqs = read_fasta(a.aln)
    alphabet = set("".join(seqs).upper())
    if alphabet <= set("ACGTUN-?"):
        raise SystemExit("GFmixF requires an amino-acid alignment")
    with open(os.path.join(a.outdir, "aln.fa"), "w", encoding="utf-8") as fh:
        for nm, s in zip(names, seqs):
            fh.write(">" + nm + "\n" + s + "\n")
    tree = parse_newick(read_tree_file(a.tree))
    if a.seed_tree:
        seed_lengths(tree, tree, a.seed_tree)
    lengths = {nid: node["length"] for nid, node in tree.nodes.items() if node["length"] is not None}
    with open(os.path.join(a.outdir, "tree.nwk"), "w", encoding="utf-8") as fh:
        fh.write(serialize_newick(tree, lengths=lengths or None) + "\n")
    root_children = tree.nodes[tree.root]["children"]
    with open(os.path.join(a.outdir, "rootfile.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(descendants(tree, root_children[0])) + "\n")
    report = open(a.iqtree, encoding="utf-8").read().splitlines()
    alpha = float(re.search(r"Gamma shape alpha:\s*([0-9.eE+-]+)", "\n".join(report)).group(1))
    comps = [ln.split() for ln in report if re.match(r"^\s*\d+\s+LG\S+\s", ln)]
    if len(comps) != a.classes + 1:
        raise SystemExit(f"IQ-TREE report has {len(comps)} components, expected {a.classes + 1}; "
                         f"make the C-class count in --model match the IQ-TREE report")
    fcomp = next(i for i, c in enumerate(comps) if re.fullmatch(r"[A-Za-z0-9]+\+F", c[1]))
    weights = [float(c[3]) for c in comps]
    w_start = [weights[i] for i in range(len(comps)) if i != fcomp] + [weights[fcomp]]
    total = sum(w_start)
    with open(os.path.join(a.outdir, "weights.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(f"{w / total:.12g}" for w in w_start) + "\n")
    with open(os.path.join(a.outdir, "alpha.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"{alpha:.12g}\n")
    prof = [ln.split() for ln in open(os.path.join(a.gfmixf_dir, f"C{a.classes}.aafreq.dat")) if ln.strip()]
    counts = {aa: 0 for aa in AA}
    for s in seqs:
        for ch in s.upper():
            if ch in counts:
                counts[ch] += 1
    total_count = sum(counts.values()) + 10
    emp = [(counts[aa] + 0.5) / total_count for aa in AA]
    with open(os.path.join(a.outdir, "freq.dat"), "w", encoding="utf-8") as fh:
        for row in prof:
            fh.write(" ".join(row) + "\n")
        fh.write(" ".join(f"{x:.12g}" for x in emp) + "\n")
    g1, g2 = a.groups.split(",")
    with open(os.path.join(a.outdir, "groups.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"{g1},{g2}\n")


if __name__ == "__main__":
    main()
