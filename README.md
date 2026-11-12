# IQ2MC-GFmix

Divergence-time estimation under the latest site-and-branch-heterogeneous GFmixF model
(McCarthy, Susko & Roger 2026, *Systematic Biology*, syag040). The workflow uses IQ-TREE
for maximum-likelihood starting values, fits GFmixF by full ML, computes the branch-length
Hessian, and samples divergence times with the patched MCMCtree under the GFmixF
approximate likelihood.

## Run

Quick start:

```bash
./iq2mc-gfmix run --aln aln.fa --tree calibrated_tree.nwk --threads 8
```

Step by step:

```bash
# GFmixF ML fit and branch-length Hessian
./iq2mc-gfmix fit --aln aln.fa --tree calibrated_tree.nwk --threads 8

# MCMCtree dating and posterior node ages
./iq2mc-gfmix date --aln aln.fa --tree calibrated_tree.nwk
```

`fit` writes the GFmixF fit, Hessian, and generated MCMCtree files to
`results/<prefix>/`; `date` adds `<prefix>.gfmix.dates.csv` with the posterior node ages.
`--outdir` and `--prefix` select another location; `--model` sets the IQ-TREE starting
model (default `LG+C20+F+G4`). `--optimize N` selects the GFmixF parameter set optimized by the engine (`1` edges, `2` +alpha, `3` +weights, `4` gammas, `5` edges+gammas, `6` all; default `6`). During the fit the L-BFGS iterations are appended to
`<prefix>/gfmix/fit_trace.log` for monitoring convergence.

For large trees the branch-length Hessian dominates the runtime. The default
`--hessian-mode analytic` builds it from full-precision analytic gradients with Richardson
extrapolation (about `4p` instead of `2p^2` evaluations); `--threads N` parallelizes the
evaluations, and `--hessian-mode numeric` is kept for verification with likelihood differences.
`--resume` continues the GFmixF fit from the last `gfmix.param` checkpoint.

The alignment must be FASTA. The tree must be rooted and
binary with calibration tags, e.g. `((A,B)'B(1.0,2.0,0.025,0.025)',C);`.

## MCMCtree parameters

`workflow/mcmctree.ctl` is the control-file template (`@PREFIX@`, `@SEQFILE@`, `@TREEFILE@`
and `@HESSIANFILE@` are filled by the workflow); edit it or pass your own with
`--mcmctree-ctl FILE`.

## References

- Demotte P. et al. (2025) IQ2MC. https://doi.org/10.32942/X2CD2X
- McCarthy C.G.P., Susko E., Roger A.J. (2026) Modeling site-and-branch-heterogeneity with
  GFmix. *Systematic Biology* syag040. Software:
  https://www.mathstat.dal.ca/~tsusko/software.html
