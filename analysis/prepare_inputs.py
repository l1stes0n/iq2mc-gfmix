"""Convert a FASTA alignment and a Newick calibration tree for the workflow."""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tree_utils import read_tree_file


def read_fasta(path):
    names, seqs, cur = [], [], -1
    for ln in open(path, encoding="utf-8"):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            names.append(re.sub(r"^>\s*", "", ln))
            seqs.append("")
            cur += 1
        elif cur >= 0:
            seqs[cur] += ln.replace(" ", "")
    if not names:
        raise SystemExit("empty FASTA alignment: " + path)
    return names, seqs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("aln")
    ap.add_argument("tree")
    ap.add_argument("outdir")
    ap.add_argument("prefix")
    a = ap.parse_args(argv)
    names, seqs = read_fasta(a.aln)
    if len({len(s) for s in seqs}) != 1:
        raise SystemExit("sequences have unequal lengths")
    if len(set(names)) != len(names):
        raise SystemExit("duplicate taxon names in alignment")
    tree_txt = read_tree_file(a.tree)
    plain = re.sub(r"'[^']*'", "", tree_txt)
    depth = 0
    for t in re.findall(r"[(),;:]", plain):
        depth += 1 if t == "(" else -1 if t == ")" else 0
        if depth < 0:
            raise SystemExit("unbalanced parentheses in tree")
    if depth != 0:
        raise SystemExit("unterminated tree")
    tips = [x.strip() for x in re.findall(r"[(,]([^(),:]+)", re.sub(r":[0-9.eE+-]+", "", plain))]
    if sorted(tips) != sorted(names):
        raise SystemExit("taxa in alignment and tree differ")
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, a.prefix + ".aln.fa"), "w", encoding="utf-8") as fh:
        for nm, s in zip(names, seqs):
            fh.write(">" + nm + "\n" + s + "\n")
    with open(os.path.join(a.outdir, a.prefix + ".tree.nwk"), "w", encoding="utf-8") as fh:
        fh.write("%d 1\n%s\n" % (len(names), tree_txt))
    with open(os.path.join(a.outdir, a.prefix + ".tree.plain.nwk"), "w", encoding="utf-8") as fh:
        fh.write(plain + "\n")
    with open(os.path.join(a.outdir, a.prefix + ".dummy.phy"), "w", encoding="utf-8") as fh:
        fh.write("%d 1\n" % len(names))
        for nm in names:
            pad = max(2, 17 - len(nm))
            fh.write(nm + " " * pad + "A\n")
    ncal = len(re.findall(r"'[^']*'", tree_txt))
    print("taxa=%d sites=%d calib=%d" % (len(names), len(seqs[0]), ncal))


if __name__ == "__main__":
    main()
