"""Summarise an MCMCtree posterior table into node ages with clades."""

import argparse
import csv
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from tree_utils import (descendants, paml_node_ids, parse_newick, read_tree_file)


def read_tips(fasta):
    names = [re.sub(r"^>\s*", "", ln) for ln in open(fasta, encoding="utf-8") if ln.startswith(">")]
    return [re.sub(r"\s.*$", "", n) for n in names]


def ess(values):
    v = [x for x in values if math.isfinite(x)]
    n = len(v)
    if n < 10:
        return float("nan")
    mean = sum(v) / n
    var = sum((x - mean) ** 2 for x in v) / n
    if var == 0:
        return float("nan")
    ac = []
    maxlag = min(200, n - 1)
    for lag in range(1, maxlag + 1):
        c = sum((v[i] - mean) * (v[i + lag] - mean) for i in range(n - lag)) / n / var
        ac.append(c)
        if c < 0.05:
            break
    return n / (1 + 2 * sum(ac))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("prefix")
    a = ap.parse_args()
    base = re.sub(r"\.gfmix$", "", a.prefix)
    out_txt = os.path.join(a.run_dir, a.prefix + ".mcmctree.out")
    tree_file = os.path.join(a.run_dir, a.prefix + ".rooted.nwk")
    if not os.path.exists(tree_file):
        tree_file = os.path.join(a.run_dir, base + ".rooted.nwk")
    aln = os.path.join(a.run_dir, "inputs", a.prefix + ".aln.fa")
    if not os.path.exists(aln):
        aln = os.path.join(a.run_dir, "inputs", base + ".aln.fa")
    tips = read_tips(aln)
    tree = parse_newick(read_tree_file(tree_file))
    ids = paml_node_ids(tree, tips)
    id_to_node = {v: k for k, v in ids.items()}
    clade = {ids[nid]: ",".join(descendants(tree, nid)) for nid in tree.nodes}
    n_taxa = {ids[nid]: len(descendants(tree, nid)) for nid in tree.nodes}
    calib = {ids[nid]: (tree.nodes[nid]["calib"] or "") for nid in tree.nodes}
    rows = []
    with open(out_txt, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    start = next((i for i, ln in enumerate(lines) if "Posterior mean (95% Equal-tail CI)" in ln), None)
    if start is None:
        raise SystemExit("no posterior table in " + out_txt)
    for ln in lines[start + 1:]:
        m = re.match(r"^t_n(\d+)\s+(.*)$", ln)
        if not m:
            continue
        vals = [float(x) for x in re.findall(r"-?\d+\.\d+", m.group(2))]
        node = int(m.group(1))
        if node not in id_to_node or len(vals) < 6:
            continue
        lo = hi = None
        lab = calib.get(node, "")
        mb = re.search(r"B\(([-0-9.eE]+),([-0-9.eE]+)", lab)
        if mb:
            lo, hi = float(mb.group(1)), float(mb.group(2))
        else:
            mg = re.search(r">([-0-9.eE]+)", lab)
            mh = re.search(r"<([-0-9.eE]+)", lab)
            if mg:
                lo = float(mg.group(1))
            if mh:
                hi = float(mh.group(1))
        rows.append({"node": node, "mean": vals[0], "et_low": vals[1], "et_high": vals[2],
                     "hpd_low": vals[3], "hpd_high": vals[4], "hpd_width": vals[5],
                     "clade": clade[node], "n_taxa": n_taxa[node], "calibration": lab,
                     "cal_low": lo, "cal_high": hi})
    rows.sort(key=lambda r: -r["mean"])
    out_csv = os.path.join(a.run_dir, a.prefix + ".dates.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    log = os.path.join(a.run_dir, a.prefix + ".mcmctree.log")
    if os.path.exists(log) and os.path.getsize(log) > 0:
        lines = open(log, encoding="utf-8").read().splitlines()
        header = lines[0].split()
        data = []
        for ln in lines[1:]:
            vals = [float(x) for x in ln.split()]
            if len(vals) == len(header):
                data.append(vals)
        if data:
            cols = {name: [row[i] for row in data] for i, name in enumerate(header)}
            root = max((c for c in cols if c.startswith("t_n")),
                       key=lambda c: sum(cols[c]) / len(cols[c]))
            diag = os.path.join(a.run_dir, a.prefix + ".mcmc_diagnostics.csv")
            with open(diag, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["parameter", "mean", "ess"])
                for key in (root, "mu", "sigma2", "lnL"):
                    if key in cols:
                        w.writerow([key, sum(cols[key]) / len(cols[key]), ess(cols[key])])
    print(out_csv)


if __name__ == "__main__":
    main()
