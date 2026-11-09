"""Gradient and Hessian of the GFmixF log-likelihood at the ML fit."""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gfmix_common import (lnl_eval, read_param, theta_utreec, write_mfile,
                          write_utreec, write_weights)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--gfmixf-dir", required=True)
    ap.add_argument("--h-rel", type=float, default=1e-4)
    a = ap.parse_args()
    d = a.fit_dir
    nclass = sum(1 for ln in open(os.path.join(d, "freq.dat"), encoding="utf-8") if ln.strip())
    names = open(os.path.join(d, "eval.names"), encoding="utf-8").read().split()
    m = len(names)
    param = read_param(os.path.join(d, "gfmix.param"), m, nclass)
    theta = param["theta"]
    p = len(theta)
    ut = param["utreec"]
    utfile = os.path.join(d, "eval.utreec")
    write_utreec(utfile, ut)
    write_weights(os.path.join(d, "eval.wtc"), param["weights"])
    write_mfile(os.path.join(d, "eval.mfile"), param["gammas"])
    args = dict(bin_path=os.path.join(a.gfmixf_dir, "gfmix_lnl_align"),
                seqfile=os.path.join(d, "eval.seqfile"),
                freq=os.path.join(d, "freq.dat"),
                gpfile=os.path.join(d, "eval.gpfile"),
                nclass=nclass, utreec_path=utfile, alpha=param["alpha"],
                wtfile=os.path.join(d, "eval.wtc"), mfile=os.path.join(d, "eval.mfile"))

    def f(th):
        write_utreec(utfile, theta_utreec(ut, th))
        return lnl_eval(**args)

    t0 = time.time()
    lnl0 = f(theta)
    h = [a.h_rel * abs(x) for x in theta]
    h = [x if x >= 1e-6 else 1e-6 for x in h]
    grad = [0.0] * p
    for i in range(p):
        hi = h[i]
        tp = list(theta); tp[i] += hi
        tm = list(theta); tm[i] -= hi
        grad[i] = (f(tp) - f(tm)) / (2 * hi)
    hess = [[0.0] * p for _ in range(p)]
    for i in range(p):
        hi = h[i]
        tp = list(theta); tp[i] += hi
        tm = list(theta); tm[i] -= hi
        hess[i][i] = (f(tp) - 2 * lnl0 + f(tm)) / (hi * hi)
    for i in range(p):
        for j in range(i + 1, p):
            hi, hj = h[i], h[j]
            th = list(theta)
            th[i] += hi; th[j] += hj; fpp = f(th)
            th = list(theta)
            th[i] += hi; th[j] -= hj; fpm = f(th)
            th = list(theta)
            th[i] -= hi; th[j] += hj; fmp = f(th)
            th = list(theta)
            th[i] -= hi; th[j] -= hj; fmm = f(th)
            hess[i][j] = hess[j][i] = (fpp - fpm - fmp + fmm) / (4 * hi * hj)
    with open(os.path.join(d, "gradient.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join("%.17g" % x for x in grad) + "\n")
    with open(os.path.join(d, "hessian.txt"), "w", encoding="utf-8") as fh:
        for row in hess:
            fh.write(" ".join("%.17g" % x for x in row) + "\n")
    open(os.path.join(d, "lnL_at_mle.txt"), "w", encoding="utf-8").write("%.17g\n" % lnl0)
    evals = 1 + 2 * p + 4 * (p * (p - 1) // 2)
    print("GFmixF Hessian: p=%d evals=%d seconds=%.1f lnL=%.6f" %
          (p, evals, time.time() - t0, lnl0))


if __name__ == "__main__":
    main()
