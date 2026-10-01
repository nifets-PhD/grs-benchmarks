#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

MODEL=kronecker
K=9
SEED=1
N=100
REF=200
CHUNKS=25
DURATION=172800
SAMPLES=10
RATES=(10 1 0.1 0.01 0.001)
CORES=${CORES:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
RUNS=runs/promoter-approximation/$MODEL
SCHEDULES=models/schedules/promoter-approximation/$MODEL
OUT=measurements/$(basename "$0" .sh)

if [ -n "${SMOKE:-}" ]; then
    TMP=$(mktemp -d); [ -n "${KEEP:-}" ] || trap 'rm -rf "$TMP"' EXIT
    echo "smoke output: $TMP" >&2
    RUNS=$TMP/runs; SCHEDULES=$TMP/schedules; OUT=$TMP/out
    K=3; N=10; REF=20; CHUNKS=2; DURATION=200; RATES=(1 0.01)
fi

num_genes=$((2 ** K))
rung=$(printf "genes-%05d" $num_genes)
ARMS=("equilibrium 1 grs-ssa-qss $N")
for rate in "${RATES[@]}"; do ARMS+=("full $rate grs-ssa $REF"); done
mkdir -p $OUT

schedule() {
    SCHEDULE=$SCHEDULES/$1-d$2-profile$3/$rung.json
    [ -f $SCHEDULE ] || julia --project=. models/kronecker.jl --out $SCHEDULE \
      --k $K --rates $MODEL --variant $1 --switching-rate $2 \
      --duration $DURATION --samples $SAMPLES --trajectories 1 --unique true \
      --seed $SEED --profile-reactions $3 < /dev/null
}

simulate() {
    [ -f $1/metrics.csv ] && return
    while [ $(jobs -r | wc -l) -ge $CORES ]; do sleep 5; done
    echo "=== [$(date -u +%H:%MZ)] $1"
    julia --project=. runners/grs.jl $1 $SCHEDULE --seed $3 --trajectories $2 --method RSSACR < /dev/null
}

ENSEMBLES=()
PROFILED=()
for spec in "${ARMS[@]}"; do
    set -- $spec
    variant=$1 rate=$2 arm=$3 ntraj=$(($4 / CHUNKS))
    schedule $variant $rate false
    dirs=""
    for c in $(seq 0 $((CHUNKS - 1))); do
        seed=$((SEED + c * ntraj))
        out=$RUNS/$variant-d$rate/$rung/$arm-n$ntraj-s$seed
        dirs="$dirs $out"
        simulate $out $ntraj $seed &
    done
    ENSEMBLES+=("$dirs")
    schedule $variant $rate true
    out=$RUNS/$variant-d$rate-profiled/$rung/$arm-n3-s$SEED
    PROFILED+=($out)
    simulate $out 3 $SEED &
done
wait

for i in "${!RATES[@]}"; do
    julia --project=. analysis/distances.jl $OUT/distances-d${RATES[$i]}.csv $num_genes \
      --nulls 20 ${ENSEMBLES[$((i + 1))]} ${ENSEMBLES[0]}
done
julia --project=. analysis/events.jl $OUT/events.csv "${PROFILED[@]}"
for species in mrnas proteins; do
    julia --project=. analysis/final-counts.jl $OUT/$species.csv .$species \
      ${ENSEMBLES[0]} ${ENSEMBLES[1]} ${ENSEMBLES[${#RATES[@]}]}
done

COSTS=()
for spec in "${ARMS[@]}"; do
    set -- $spec
    schedule $1 $2 false
    out=$RUNS/$1-d$2/$rung/$3-n3-s$SEED
    COSTS+=($out)
    simulate $out 3 $SEED
done
julia --project=. analysis/costs.jl $OUT/costs.csv "${COSTS[@]}"
