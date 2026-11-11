"""Full maximum-likelihood GFmixF fit through the C engine."""

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
sys.path.insert(0, HERE)
from gfmix_common import (node_clades, read_fasta, read_param, write_gpfile,
                          write_integer_seqfile, write_mfile, write_node_clades,
                          write_utreec, write_weights)
from tree_utils import parse_newick, read_tree_file, serialize_newick


def split_leading1(splits):
    for row in splits:
        if row[0] == 0:
            row[:] = [1 - x for x in row]
    return splits


def utreec_splits(ut, ntaxa):
    splits = [[0] * ntaxa for _ in range(ntaxa - 3)]
    for s in range(ntaxa - 3):
        for k in (0, 1):
            b = ut[s][k]
            if b < ntaxa:
                splits[s][b] = 1
            else:
                for i in range(ntaxa):
                    splits[s][i] += splits[b - ntaxa][i]
    return split_leading1(splits)


def locate_root(ut, taxa):
    ntaxa = len(ut) + 1
    if len(taxa) == 1:
        return taxa[0]
    if 0 in taxa:
        want = [1 if i in taxa else 0 for i in range(ntaxa)]
    else:
        want = [0 if i in taxa else 1 for i in range(ntaxa)]
    for i, row in enumerate(utreec_splits(ut, ntaxa)):
        if row == want:
            return ntaxa + i
    raise SystemExit("root split not found in tree")


def tree_to_int(tree_path, names, path):
    tree = parse_newick(read_tree_file(tree_path))
    ids = {nm: i for i, nm in enumerate(names)}
    lengths = {}
    for nid, node in tree.nodes.items():
        if node["length"] is not None:
            lengths[nid] = node["length"]
        node["calib"] = None
        if node["children"]:
            node["name"] = None
        else:
            node["name"] = str(ids[node["name"]])
    open(path, "w", encoding="utf-8").write(serialize_newick(tree, lengths=lengths or None) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--gfmixf-dir", required=True)
    ap.add_argument("--classes", type=int, required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--threads", default="2")
    ap.add_argument("--optimize", type=int, default=6, choices=range(1, 7))
    ap.add_argument("--start-from", default="")
    a = ap.parse_args(argv)
    d = a.fit_dir
    names, seqs = read_fasta(os.path.join(d, "aln.fa"))
    m = len(names)
    ids = {nm: i for i, nm in enumerate(names)}
    seqfile = os.path.join(d, "eval.seqfile")
    write_integer_seqfile(seqfile, names, seqs, ids)
    tree_int = os.path.join(d, "tree.int")
    tree_to_int(os.path.join(d, "tree.nwk"), names, tree_int)
    utcu = os.path.join(d, "eval.utreecu")
    treecns = os.path.join(a.gfmixf_dir, "treecns")
    rert = os.path.join(a.gfmixf_dir, "rert")
    gfmix = os.path.join(a.gfmixf_dir, "gfmix_lnl_align")
    names_out = os.path.join(d, "eval.names.treecns")
    with open(seqfile, "rb") as fin, open(names_out, "wb") as fout:
        rc = subprocess.run([treecns, tree_int, utcu, "-1"], stdin=fin, stdout=fout).returncode
    if rc != 0 or not os.path.exists(utcu):
        raise SystemExit("treecns failed")
    ut0 = [[int(ln.split()[0]), int(ln.split()[1]), float(ln.split()[2]), float(ln.split()[3])]
           for ln in open(utcu, encoding="utf-8") if len(ln.split()) == 4]
    root_taxa = [ids[nm] for nm in open(os.path.join(d, "rootfile.txt"), encoding="utf-8").read().split()]
    utfile = os.path.join(d, "eval.utreec")
    label = locate_root(ut0, root_taxa)
    with open(utcu, "rb") as fin, open(utfile, "wb") as fout:
        rc = subprocess.run([rert, str(label), str(m)], stdin=fin, stdout=fout).returncode
    if rc != 0:
        raise SystemExit("rert failed")
    ut = [[int(ln.split()[0]), int(ln.split()[1]), float(ln.split()[2]), float(ln.split()[3])]
          for ln in open(utfile, encoding="utf-8") if len(ln.split()) == 4]
    gpfile = os.path.join(d, "eval.gpfile")
    g1, g2 = a.groups.split(",")
    write_gpfile(gpfile, g1, g2)
    start_u = utfile
    start_w = os.path.join(d, "weights.txt")
    start_a = open(os.path.join(d, "alpha.txt")).read().strip()
    start_m = None
    if a.start_from:
        restart = read_param(a.start_from, m, a.classes + 1)
        start_u = os.path.join(d, "start.utreec")
        write_utreec(start_u, restart["utreec"])
        start_w = os.path.join(d, "start.wtc")
        write_weights(start_w, restart["weights"])
        start_a = "%.17g" % restart["alpha"]
        start_m = os.path.join(d, "start.mfile")
        write_mfile(start_m, restart["gammas"])
    param = os.path.join(d, "gfmix.param")
    fit_cmd = [gfmix, "-s", seqfile, "-c", str(a.classes + 1),
               "-f", os.path.join(d, "freq.dat"), "-u", start_u, "-g", gpfile,
               "-w", start_w, "-a", start_a,
               "-e", "9", "-o", str(a.optimize), "-b", str(a.threads), "-p", "1", "-x", param]
    if start_m:
        fit_cmd += ["-m", start_m]
    trace_path = os.path.join(d, "fit_trace.log")
    lines = []
    with open(trace_path, "w", encoding="utf-8") as trace:
        proc = subprocess.Popen(fit_cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            trace.write(line)
            trace.flush()
            lines.append(line)
        rc = proc.wait()
    if rc != 0 or not os.path.exists(param):
        raise SystemExit("GFmixF fit failed: " + "".join(lines)[-500:])
    nums = [x for x in "".join(lines).split() if re.fullmatch(r"[-+0-9.eE]+", x)]
    lnl = float(nums[-1]) if nums else float("nan")
    p = read_param(param, m, a.classes + 1)
    theta = p["theta"]
    with open(os.path.join(d, "theta.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join("%.17g" % x for x in theta) + "\n")
    open(os.path.join(d, "lnL_at_mle.txt"), "w", encoding="utf-8").write("%.17g\n" % lnl)
    open(os.path.join(d, "alpha.txt"), "w", encoding="utf-8").write("%.17g\n" % p["alpha"])
    write_weights(os.path.join(d, "eval.wtc"), p["weights"])
    write_weights(os.path.join(d, "weights_fit.txt"), p["weights"])
    write_mfile(os.path.join(d, "eval.mfile"), p["gammas"])
    open(os.path.join(d, "eval.names"), "w", encoding="utf-8").write("\n".join(names) + "\n")
    write_utreec(utfile, ut)
    write_node_clades(os.path.join(d, "node_clades.txt"), node_clades(ut, m, names))
    with open(os.path.join(d, "gammas.txt"), "w", encoding="utf-8") as fh:
        fh.write("node\tgamma_G\tgamma_F\tlabel\n")
        for node, row in enumerate(p["gammas"]):
            label = names[node] if node < m else ("root" if node == 2 * m - 2 else "node%d" % node)
            fh.write("%d\t%.17g\t%.17g\t%s\n" % (node, row[0], row[1], label))
    print("GFmixF fit: m=%d nclass=%d lnL=%.6f" % (m, a.classes + 1, lnl))


if __name__ == "__main__":
    main()
