"""Gradient and Hessian of the GFmixF log-likelihood at the ML fit."""

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gfmix_common import (read_param, theta_utreec, write_mfile, write_utreec,
                          write_weights)


def run_engine(bin_path, seqfile, freq, gpfile, nclass, utfile, alpha, wtfile, mfile, deriv):
    cmd = [bin_path, "-s", seqfile, "-c", str(nclass), "-f", freq, "-g", gpfile,
           "-e", "9", "-u", utfile, "-a", "%.17g" % alpha, "-w", wtfile,
           "-m", mfile, "-n", str(deriv)]
    if deriv >= 1:
        cmd += ["-q", "1"]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        raise SystemExit("gfmix_lnl_align failed: " + (out.stdout + out.stderr)[-500:])
    return out.stdout


def parse_lnl(text):
    vals = [x for x in text.split() if x.replace(".", "").replace("-", "").replace("+", "")
            .replace("e", "").replace("E", "").isdigit()]
    return float(vals[-1])


def parse_gradient(text, p):
    lines = text.splitlines()
    i0 = next((i for i, ln in enumerate(lines) if "Edge Length Derivatives" in ln), None)
    if i0 is None:
        raise SystemExit("no gradient table in engine output")
    grad = [None] * p
    for ln in lines[i0 + 1:]:
        fld = ln.split()
        if len(fld) == 3 and fld[0].lstrip("+-").isdigit():
            e = int(fld[0])
            if 0 <= e < p:
                grad[e] = float(fld[2])
    if any(g is None for g in grad):
        raise SystemExit("incomplete gradient table in engine output")
    return grad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--gfmixf-dir", required=True)
    ap.add_argument("--hessian-mode", choices=["numeric", "analytic"], default="analytic")
    ap.add_argument("--h-rel", type=float, default=1e-4)
    ap.add_argument("--jobs", type=int, default=1)
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
    engine = os.path.join(a.gfmixf_dir, "gfmix_lnl_align")
    seqfile = os.path.join(d, "eval.seqfile")
    freq = os.path.join(d, "freq.dat")
    gpfile = os.path.join(d, "eval.gpfile")
    wtfile = os.path.join(d, "eval.wtc")
    mfile = os.path.join(d, "eval.mfile")

    def f(th):
        write_utreec(utfile, theta_utreec(ut, th))
        text = run_engine(engine, seqfile, freq, gpfile, nclass, utfile, param["alpha"],
                          wtfile, mfile, 0)
        return parse_lnl(text)

    def g(th):
        write_utreec(utfile, theta_utreec(ut, th))
        text = run_engine(engine, seqfile, freq, gpfile, nclass, utfile, param["alpha"],
                          wtfile, mfile, 1)
        return parse_lnl(text), parse_gradient(text, p)

    t0 = time.time()
    h = [a.h_rel * abs(x) for x in theta]
    h = [x if x >= 1e-6 else 1e-6 for x in h]
    if a.hessian_mode == "numeric":
        eval_dir = os.path.join(d, "hessian_eval")
        os.makedirs(eval_dir, exist_ok=True)
        base_path = os.path.join(eval_dir, "ut_base.utreec")

        def f_at(th, path):
            write_utreec(path, theta_utreec(ut, th))
            text = run_engine(engine, seqfile, freq, gpfile, nclass, path, param["alpha"],
                              wtfile, mfile, 0)
            return parse_lnl(text)

        def column_task(i):
            path = os.path.join(eval_dir, "ut_col_%d.utreec" % i)
            hi = h[i]
            tp = list(theta); tp[i] += hi
            tm = list(theta); tm[i] -= hi
            fp = f_at(tp, path)
            fm = f_at(tm, path)
            grad_i = (fp - fm) / (2 * hi)
            diag = (fp - 2 * lnl0 + fm) / (hi * hi)
            row = []
            for j in range(i + 1, p):
                hj = h[j]
                th = list(theta); th[i] += hi; th[j] += hj; fpp = f_at(th, path)
                th = list(theta); th[i] += hi; th[j] -= hj; fpm = f_at(th, path)
                th = list(theta); th[i] -= hi; th[j] += hj; fmp = f_at(th, path)
                th = list(theta); th[i] -= hi; th[j] -= hj; fmm = f_at(th, path)
                row.append((j, (fpp - fpm - fmp + fmm) / (4 * hi * hj)))
            return i, grad_i, diag, row

        lnl0 = f_at(theta, base_path)
        grad = [0.0] * p
        hess = [[0.0] * p for _ in range(p)]
        done = 0
        step = max(1, p // 10)
        if a.jobs > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
                futures = [pool.submit(column_task, i) for i in range(p)]
                for fut in concurrent.futures.as_completed(futures):
                    i, grad_i, diag, row = fut.result()
                    grad[i] = grad_i
                    hess[i][i] = diag
                    for j, value in row:
                        hess[i][j] = hess[j][i] = value
                    done += 1
                    if done % step == 0:
                        print("numeric Hessian %d/%d" % (done, p), file=sys.stderr, flush=True)
        else:
            for i in range(p):
                i, grad_i, diag, row = column_task(i)
                grad[i] = grad_i
                hess[i][i] = diag
                for j, value in row:
                    hess[i][j] = hess[j][i] = value
                done += 1
                if done % step == 0:
                    print("numeric Hessian %d/%d" % (done, p), file=sys.stderr, flush=True)
        shutil.rmtree(eval_dir, ignore_errors=True)
        evals = 1 + 2 * p + 4 * (p * (p - 1) // 2)
    else:
        eval_dir = os.path.join(d, "hessian_eval")
        os.makedirs(eval_dir, exist_ok=True)

        def gradient_task(index, sign, scale):
            th = list(theta)
            th[index] += sign * h[index] * scale
            tag = "%d_%d_%d" % (index, 0 if sign < 0 else 1, 0 if scale < 1.0 else 1)
            path = os.path.join(eval_dir, "ut_%s.utreec" % tag)
            write_utreec(path, theta_utreec(ut, th))
            text = run_engine(engine, seqfile, freq, gpfile, nclass, path, param["alpha"],
                              wtfile, mfile, 1)
            return index, sign, scale, parse_gradient(text, p)

        lnl0, grad = g(theta)
        hess = [[0.0] * p for _ in range(p)]
        tasks = [(i, s, sc) for i in range(p) for s in (-1, 1) for sc in (1.0, 0.5)]
        results = {}
        done = 0
        step = max(1, p // 10)
        if a.jobs > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
                futures = [pool.submit(gradient_task, i, s, sc) for i, s, sc in tasks]
                for fut in concurrent.futures.as_completed(futures):
                    i, s, sc, gv = fut.result()
                    results[(i, s, sc)] = gv
                    done += 1
                    if done % (4 * step) == 0:
                        print("analytic Hessian %d/%d" % (done // 4, p), file=sys.stderr, flush=True)
        else:
            for i, s, sc in tasks:
                _, _, _, gv = gradient_task(i, s, sc)
                results[(i, s, sc)] = gv
                done += 1
                if done % (4 * step) == 0:
                    print("analytic Hessian %d/%d" % (done // 4, p), file=sys.stderr, flush=True)
        for i in range(p):
            hi = h[i]
            gp1 = results[(i, 1, 1.0)]
            gm1 = results[(i, -1, 1.0)]
            gph = results[(i, 1, 0.5)]
            gmh = results[(i, -1, 0.5)]
            for j in range(p):
                ref = (gp1[j] - gm1[j]) / (2 * hi)
                half = (gph[j] - gmh[j]) / hi
                hess[j][i] = (4 * half - ref) / 3
        for i in range(p):
            for j in range(i + 1, p):
                value = 0.5 * (hess[i][j] + hess[j][i])
                hess[i][j] = hess[j][i] = value
        shutil.rmtree(eval_dir, ignore_errors=True)
        evals = 4 * p + 1
    with open(os.path.join(d, "gradient.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join("%.17g" % x for x in grad) + "\n")
    with open(os.path.join(d, "hessian.txt"), "w", encoding="utf-8") as fh:
        for row in hess:
            fh.write(" ".join("%.17g" % x for x in row) + "\n")
    open(os.path.join(d, "lnL_at_mle.txt"), "w", encoding="utf-8").write("%.17g\n" % lnl0)
    print("GFmixF Hessian: mode=%s p=%d evals=%d seconds=%.1f lnL=%.6f" %
          (a.hessian_mode, p, evals, time.time() - t0, lnl0))


if __name__ == "__main__":
    main()
