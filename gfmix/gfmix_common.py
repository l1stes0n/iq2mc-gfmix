"""Shared helpers for the GFmixF step."""

import os
import re
import subprocess

AA = "ARNDCQEGHILKMFPSTWYV"


def read_fasta(path):
    names, seqs, cur = [], [], -1
    for ln in open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if not ln:
            continue
        if ln.startswith(">"):
            names.append(re.sub(r"^>\s*", "", ln))
            seqs.append("")
            cur += 1
        else:
            seqs[cur] += re.sub(r"\s+", "", ln)
    return names, seqs


def write_integer_seqfile(path, names, seqs, ids):
    """Single-line sequence file with integer taxon names (gfmixf/C engine input)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{len(names)} {len(seqs[0])}\n")
        for nm, s in zip(names, seqs):
            fh.write(f"{ids[nm]:<10}{s}\n")


def write_gpfile(path, garp, fymink):
    cia = [2] * 20
    for ch in garp:
        cia[AA.index(ch)] = 0
    for ch in fymink:
        cia[AA.index(ch)] = 1
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(" ".join(str(x) for x in cia) + "\n")


def write_utreec(path, ut):
    with open(path, "w", encoding="utf-8") as fh:
        for row in ut:
            fh.write("%d %d %.17g %.17g\n" % (int(row[0]), int(row[1]), row[2], row[3]))


def utreec_theta(ut):
    m = len(ut) + 1
    theta = [0.0] * (2 * m - 2)
    for row in ut:
        theta[row[0]] = row[2]
        theta[row[1]] = row[3]
    return theta


def theta_utreec(ut, theta):
    return [[row[0], row[1], theta[row[0]], theta[row[1]]] for row in ut]


def read_param(path, m, nclass):
    vals = [float(x) for x in open(path, encoding="utf-8").read().split()]
    n_ut = 4 * (m - 1)
    if len(vals) != 1 + n_ut + nclass + (2 * m - 1) * 2:
        raise SystemExit("unexpected parameter file layout: " + path)
    ut = [[int(vals[1 + 4 * j]), int(vals[2 + 4 * j]), vals[3 + 4 * j], vals[4 + 4 * j]]
          for j in range(m - 1)]
    p = 1 + n_ut
    weights = vals[p:p + nclass]
    p += nclass
    gammas = [vals[p + 2 * j:p + 2 * j + 2] for j in range(2 * m - 1)]
    return {"alpha": vals[0], "utreec": ut, "weights": weights, "gammas": gammas,
            "m": m, "nclass": nclass, "theta": utreec_theta(ut)}


def lnl_eval(bin_path, seqfile, freq, gpfile, nclass, utreec_path, alpha=None,
             wtfile=None, mfile=None):
    cmd = [bin_path, "-s", seqfile, "-c", str(nclass), "-f", freq, "-g", gpfile,
           "-e", "9", "-n", "0", "-u", utreec_path]
    if alpha is not None:
        cmd += ["-a", f"{alpha:.17g}"]
    if wtfile:
        cmd += ["-w", wtfile]
    if mfile:
        cmd += ["-m", mfile]
    out = subprocess.run(cmd, capture_output=True, text=True)
    vals = [x for x in out.stdout.split() if re.fullmatch(r"[-+0-9.eE]+", x)]
    if out.returncode != 0 or not vals:
        raise SystemExit("gfmix_lnl_align failed: " + (out.stdout + out.stderr)[-500:])
    return float(vals[-1])


def write_weights(path, weights):
    with open(path, "w", encoding="utf-8") as fh:
        for w in weights:
            fh.write(f"{w:.17g}\n")


def write_mfile(path, gammas):
    with open(path, "w", encoding="utf-8") as fh:
        for row in gammas:
            fh.write("%.17g %.17g\n" % (row[0], row[1]))


def node_parent(ut, m):
    parent = {}
    for j, row in enumerate(ut):
        p = m + j
        parent[row[0]] = p
        parent[row[1]] = p
    return parent


def node_clades(ut, m, names):
    parent = node_parent(ut, m)
    dl = {i: [i] for i in range(m)}
    for j, row in enumerate(ut):
        dl[m + j] = dl[row[0]] + dl[row[1]]
    out = []
    for node in range(2 * m - 1):
        taxa = " ".join(names[t] for t in dl[node])
        out.append((node, parent.get(node, -1), 1 if node < m else 0, taxa))
    return out


def write_node_clades(path, nc):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("node\tparent\ttip\tclade\n")
        for row in nc:
            fh.write("%d\t%d\t%d\t%s\n" % row)


def read_node_clades(path):
    rows = {}
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows[int(f[0])] = {"parent": int(f[1]), "tip": int(f[2]), "clade": f[3].split()}
    return rows
