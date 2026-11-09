#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IQTREE="$ROOT/third_party/iqtree3/bin/iqtree3"
MCMCTREE_NONREV="$ROOT/third_party/paml-iq2mc/bin/mcmctree-nonrev"
GFMIXF="$ROOT/third_party/gfmixf"
PYTHON="${PYTHON:-python3}"

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
    cat <<EOF
usage: bash workflow/iq2mc-gfmix.sh <fit|date|run> --aln FILE --tree FILE [options]

  fit       GFmixF ML fit + branch-length Hessian
  date      MCMCtree dating under the GFmixF approximate likelihood
  run       fit then date

  --model STRING         IQ-TREE model for starting values (default: LG+C20+F+G4)
  --prefix NAME          output prefix (default: alignment basename)
  --outdir DIR           output directory (default: results)
  --threads N            threads for IQ-TREE and GFmixF (default: 4)
  --seed N               IQ-TREE seed (default: 1)
  --mcmctree-ctl FILE    MCMCtree control template (default: workflow/mcmctree.ctl)
  --force                rerun existing outputs
EOF
}

CMD="run"
if [ $# -ge 1 ]; then CMD="$1"; shift; fi
case "$CMD" in
    fit|date|run) ;;
    -h|--help|help) usage; exit 0 ;;
    *) usage >&2; die "unknown command: $CMD" ;;
esac

ALN=""; TREE=""; MODEL="LG+C20+F+G4"; PREFIX=""; OUTDIR="results"; THREADS=4; SEED=1
FORCE=0; MCMCTREE_CTL=""
while [ $# -gt 0 ]; do
    case "$1" in
        --aln)          ALN="$2"; shift 2 ;;
        --tree)         TREE="$2"; shift 2 ;;
        --model)        MODEL="$2"; shift 2 ;;
        --prefix)       PREFIX="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --seed)         SEED="$2"; shift 2 ;;
        --mcmctree-ctl) MCMCTREE_CTL="$2"; shift 2 ;;
        --force)        FORCE=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) usage >&2; die "unknown option: $1" ;;
    esac
done

[ -n "$ALN" ] && [ -f "$ALN" ] || { usage >&2; die "--aln is required"; }
[ -n "$TREE" ] && [ -f "$TREE" ] || { usage >&2; die "--tree is required"; }
[ -n "$PREFIX" ] || PREFIX="$(basename "${ALN%.*}")"
[ -n "$MCMCTREE_CTL" ] || MCMCTREE_CTL="$ROOT/workflow/mcmctree.ctl"
[ -f "$MCMCTREE_CTL" ] || die "MCMCtree control template not found: $MCMCTREE_CTL"

CLASSES="$(printf '%s' "$MODEL" | sed -n 's/.*C\([0-9][0-9]*\).*/\1/p')"
[ -n "$CLASSES" ] || CLASSES=20

RUN_DIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)/$PREFIX"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/inputs"
ALN_FA="$RUN_DIR/inputs/$PREFIX.aln.fa"
TREE_NWK="$RUN_DIR/inputs/$PREFIX.tree.nwk"
TREE_PLAIN="$RUN_DIR/inputs/$PREFIX.tree.plain.nwk"
G="$RUN_DIR/gfmix"

stage_prepare() {
    "$PYTHON" "$ROOT/analysis/prepare_inputs.py" "$ALN" "$TREE" "$RUN_DIR/inputs" "$PREFIX" \
        > "$RUN_DIR/logs/prepare.log" 2>&1 || { cat "$RUN_DIR/logs/prepare.log" >&2; die "prepare failed"; }
}

stage_tree() {
    if [ -s "$RUN_DIR/$PREFIX.iqtree" ] && [ "$FORCE" = 0 ] \
        && [ "$(cat "$RUN_DIR/logs/model.txt" 2>/dev/null)" = "$MODEL" ]; then return 0; fi
    [ -s "$ALN_FA" ] || stage_prepare
    ( cd "$RUN_DIR" && "$IQTREE" -redo -s "$ALN_FA" -m "$MODEL" -te "$TREE_PLAIN" \
        --prefix "$PREFIX" -T "$THREADS" --seed "$SEED" > logs/iqtree.log 2>&1 ) \
        || { tail -n 30 "$RUN_DIR/logs/iqtree.log" >&2; die "IQ-TREE failed"; }
    echo "$MODEL" > "$RUN_DIR/logs/model.txt"
}

stage_fit() {
    mkdir -p "$G"
    stage_prepare
    stage_tree
    if [ ! -s "$G/hessian.txt" ] || [ "$FORCE" = 1 ]; then
        "$PYTHON" "$ROOT/gfmix/prep_gfmix_inputs.py" --aln "$ALN_FA" --tree "$TREE_NWK" \
            --iqtree "$RUN_DIR/$PREFIX.iqtree" --gfmixf-dir "$GFMIXF" --outdir "$G" \
            --classes "$CLASSES" --groups GARP,FYMINK --seed-tree "$RUN_DIR/$PREFIX.treefile" \
            > "$G/prep.log" 2>&1 || { tail -n 20 "$G/prep.log" >&2; die "GFmixF input preparation failed"; }
        "$PYTHON" "$ROOT/gfmix/gfmix_fit.py" --fit-dir "$G" --gfmixf-dir "$GFMIXF" \
            --classes "$CLASSES" --groups GARP,FYMINK --threads "$THREADS" > "$G/fit.log" 2>&1 \
            || { tail -n 20 "$G/fit.log" >&2; die "GFmixF fit failed"; }
        "$PYTHON" "$ROOT/gfmix/gfmix_hessian.py" --fit-dir "$G" --gfmixf-dir "$GFMIXF" \
            > "$G/hessian.log" 2>&1 || { tail -n 20 "$G/hessian.log" >&2; die "GFmixF Hessian failed"; }
        cat "$G/fit.log" "$G/hessian.log"
    fi
    cp "$TREE_NWK" "$RUN_DIR/$PREFIX.rooted.nwk"
    cp "$RUN_DIR/inputs/$PREFIX.dummy.phy" "$RUN_DIR/$PREFIX.dummy.phy"
    "$PYTHON" "$ROOT/gfmix/gfmix_to_mcmctree.py" --tree "$TREE_NWK" --aln "$ALN_FA" \
        --fit-dir "$G" --ctl-template "$MCMCTREE_CTL" \
        --ctl-out "$RUN_DIR/$PREFIX.gfmix.mcmctree.ctl" \
        --hessian-out "$RUN_DIR/$PREFIX.gfmix.mcmctree.hessian" \
        --out-prefix "$PREFIX.gfmix" --base-prefix "$PREFIX" > "$G/to_mcmctree.log" 2>&1 \
        || { cat "$G/to_mcmctree.log" >&2; die "Hessian conversion failed"; }
}

stage_date() {
    [ -s "$RUN_DIR/$PREFIX.gfmix.mcmctree.hessian" ] || die "run the fit step first"
    if [ -s "$RUN_DIR/$PREFIX.gfmix.dates.csv" ] && [ "$FORCE" = 0 ]; then return 0; fi
    rm -f "$RUN_DIR/$PREFIX.gfmix.mcmctree.ckp" "$RUN_DIR/$PREFIX.gfmix.mcmctree.log" \
          "$RUN_DIR/$PREFIX.gfmix.mcmctree.out" "$RUN_DIR/FigTree.tre"
    ( cd "$RUN_DIR" && "$MCMCTREE_NONREV" "$PREFIX.gfmix.mcmctree.ctl" > "$G/mcmctree.log" 2>&1 ) \
        || { tail -n 20 "$G/mcmctree.log" >&2; die "MCMCtree failed"; }
    "$PYTHON" "$ROOT/analysis/summarize_mcmctree.py" "$RUN_DIR" "$PREFIX.gfmix" \
        > "$G/summary.log" 2>&1 || { tail -n 20 "$G/summary.log" >&2; die "summary failed"; }
}

case "$CMD" in
    fit)  stage_fit ;;
    date) stage_date ;;
    run)  stage_fit; stage_date ;;
esac
