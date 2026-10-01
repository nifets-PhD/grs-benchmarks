#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

MODEL=kronecker-full
K=7
SEED=1
N=200
REF=400
DURATION=172800
SAMPLES=10
CORES=${CORES:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
PLANNED=()
RUNS=runs
SCHEDULES=models/schedules
BUILDS=models/build
OUT=measurements/$(basename "$0" .sh).csv

if [ -n "${SMOKE:-}" ]; then
    TMP=$(mktemp -d); [ -n "${KEEP:-}" ] || trap 'rm -rf "$TMP"' EXIT
    echo "smoke output: $TMP" >&2
    RUNS=$TMP/runs; SCHEDULES=$TMP/schedules; BUILDS=$TMP/build; OUT=$TMP/out.csv
    K=3; N=20; REF=40; DURATION=200
fi

TAUS=(10 3 1 0.3 0.1)
EPS=0.05

ARMS=(
  "grs-ssa                 $REF   1  grs     --method RSSACR"
  "gillespiessa2-ssa       $N    20  g2      --method exact"
  "copasi-ssa              $N     4  copasi  --method stochastic"
  "copasi-hybrid-eps$EPS   $N     8  copasi  --method adaptiveSA --epsilon $EPS"
  "copasi-tauleap-eps$EPS  $N     4  copasi  --method tauLeap --epsilon $EPS"
)

for t in "${TAUS[@]}"; do
    ARMS+=("grs-tauleap-tau$t            $N  1  grs  --method FixedTau --tau $t --policy AlwaysLeap --thin-negative false")
    ARMS+=("grs-hybrid-tau$t             $N  1  grs  --method FixedTau --tau $t --policy CriticalBlend --nc 2")
    ARMS+=("gillespiessa2-tauleap-tau$t  $N  4  g2   --method etl --tau $t")
done

num_genes=$((2 ** K))
rung=$(printf "genes-%05d" $num_genes)
SCHEDULE=$SCHEDULES/$MODEL/$rung.json
BUILD=$BUILDS/$MODEL/$rung

prepare() {
    local k=$1
    mkdir -p $(dirname $SCHEDULE) $(dirname $BUILD)
    [ -f $SCHEDULE ] || julia --project=. models/kronecker.jl --out $SCHEDULE \
      --k $k --variant full \
      --duration $DURATION --samples $SAMPLES --trajectories 1 --unique false --copies 1
    [ -f $BUILD.sbml.xml ]    || { julia --project=. models/export.jl sbml      $SCHEDULE $BUILD.sbml.tmp      $SEED && mv $BUILD.sbml.tmp $BUILD.sbml.xml; }
    [ -f $BUILD.gillespie.R ] || { julia --project=. models/export.jl gillespie $SCHEDULE $BUILD.gillespie.tmp $SEED && mv $BUILD.gillespie.tmp $BUILD.gillespie.R; }
}

run () {
    local out=$1 runner=$2 ntraj=$3 seed=$4; shift 4
    case $runner in
        grs) julia --project=. runners/grs.jl $out $SCHEDULE --seed $seed --trajectories $ntraj "$@" ;;
        g2) Rscript runners/gillespie.R $out $BUILD.gillespie.R --seed $seed --trajectories $ntraj --duration $DURATION --samples $SAMPLES --stop-on-neg FALSE "$@" ;;
        copasi) .venv/bin/python runners/copasi.py $out $BUILD.sbml.xml --seed $seed --trajectories $ntraj --duration $DURATION --samples $SAMPLES "$@" ;;
    esac
}

prepare $K

for spec in "${ARMS[@]}"; do
    set -- $spec
    arm=$1 total=$2 chunks=$3 runner=$4; shift 4
    ntraj=$((total / chunks))
    for c in $(seq 0 $((chunks - 1))); do
        seed=$((SEED + c * ntraj))
        out=$RUNS/$MODEL/$rung/$arm-n$ntraj-s$seed
        PLANNED+=("$out")
        [ -f $out/metrics.csv ] && continue
        while [ $(jobs -r | wc -l) -ge $CORES ]; do sleep 5; done
        echo "=== [$(date -u +%H:%MZ)] $arm n=$ntraj s=$seed"
        run $out $runner $ntraj $seed "$@" &
    done
done
wait || echo "some jobs failed" >&2

julia --project=. analysis/distances.jl $OUT $num_genes --nulls 20 "${PLANNED[@]}"

COSTS=()
for spec in "${ARMS[@]}"; do
    set -- $spec
    arm=$1 runner=$4; shift 4
    out=$RUNS/$MODEL/$rung/$arm-n3-s$SEED
    COSTS+=("$out")
    [ -f $out/metrics.csv ] && continue
    echo "=== [$(date -u +%H:%MZ)] cost $arm"
    run $out $runner 3 $SEED "$@"
done

julia --project=. analysis/costs.jl ${OUT%.csv}-cost.csv "${COSTS[@]}"
