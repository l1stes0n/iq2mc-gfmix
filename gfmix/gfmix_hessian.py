"""Branch-length Hessian of the GFmixF log-likelihood at the ML fit."""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gfmix_common import read_param, write_mfile, write_utreec, write_weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--gfmixf-dir", required=True)
    a = ap.parse_args()
    d = a.fit_dir
    nclass = sum(1 for ln in open(os.path.join(d, "freq.dat"), encoding="utf-8") if ln.strip())
    m = len(open(os.path.join(d, "eval.names"), encoding="utf-8").read().split())
    param = read_param(os.path.join(d, "gfmix.param"), m, nclass)
    utfile = os.path.join(d, "eval.utreec")
    wtfile = os.path.join(d, "eval.wtc")
    mfile = os.path.join(d, "eval.mfile")
    write_utreec(utfile, param["utreec"])
    write_weights(wtfile, param["weights"])
    write_mfile(mfile, param["gammas"])
    hfile = os.path.join(d, "hessian_engine.txt")
    t0 = time.time()
    subprocess.run([os.path.join(a.gfmixf_dir, "gfmix_lnl_align"),
                    "-s", os.path.join(d, "eval.seqfile"), "-c", str(nclass),
                    "-f", os.path.join(d, "freq.dat"), "-u", utfile,
                    "-g", os.path.join(d, "eval.gpfile"), "-e", "9", "-n", "0",
                    "-a", "%.17g" % param["alpha"], "-w", wtfile, "-m", mfile,
                    "-H", hfile], check=True)
    rows = [ln.split() for ln in open(hfile, encoding="utf-8") if ln.strip()]
    os.remove(hfile)
    with open(os.path.join(d, "gradient.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows[3]) + "\n")
    with open(os.path.join(d, "hessian.txt"), "w", encoding="utf-8") as fh:
        for row in rows[4:]:
            fh.write(" ".join(row) + "\n")
    lnl = float(rows[1][0])
    open(os.path.join(d, "lnL_at_mle.txt"), "w", encoding="utf-8").write("%.17g\n" % lnl)
    print("GFmixF Hessian: p=%d lnL=%.6f seconds=%.1f" % (len(rows) - 4, lnl, time.time() - t0))


if __name__ == "__main__":
    main()
