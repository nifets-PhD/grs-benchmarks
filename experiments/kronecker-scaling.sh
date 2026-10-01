#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

MODEL=kronecker-full
N=3
SEED=1
KS=(3 4 5 6 7 8 9 10 11 12 13 14 15)
DURATION=172800
SAMPLES=10
PLANNED=()
RUNS=runs
SCHEDULES=models/schedules
BUILDS=models/build
OUT=measurements/$(basename "$0" .sh)${ONLY:+-$ONLY}.csv

if [ -n "${SMOKE:-}" ]; then
    TMP=$(mktemp -d); [ -n "${KEEP:-}" ] || trap 'rm -rf "$TMP"' EXIT
    echo "smoke output: $TMP" >&2
    RUNS=$TMP/runs; SCHEDULES=$TMP/schedules; BUILDS=$TMP/build; OUT=$TMP/out.csv
    KS=(3); DURATION=200
fi

ARMS=(
  "grs-ssa              32768  grs     --method RSSACR"
  "grs-hybrid-tau3      32768  grs     --method FixedTau --tau 3 --policy CriticalBlend --nc 2"
  "grs-tauleap-tau3     32768  grs     --method FixedTau --tau 3 --policy AlwaysLeap --thin-negative false"
  "grs-tauleap-eps0.05  32768  grs     --method AdaptiveTau --epsilon 0.05"
  "gillespiessa2-ssa      128  g2      --method exact"
  "gillespiessa2-tauleap-tau3 32768 g2 --method etl --tau 3"
  "copasi-ssa                   1024  copasi  --method stochastic"
  "copasi-hybrid-eps0.05         256  copasi  --method adaptiveSA --epsilon 0.05"
  "copasi-tauleap-eps0.05       2048  copasi  --method tauLeap --epsilon 0.05"
)

prepare() {
    local k=$1
    mkdir -p $(dirname $SCHEDULE) $(dirname $BUILD)
    [ -f $SCHEDULE ] || julia --project=. models/kronecker.jl --out $SCHEDULE \
      --k $k --variant full \
      --duration $DURATION --samples $SAMPLES --trajectories 1 --unique false --copies 1
    [ -f $BUILD.sbml.xml ]    || { julia --project=. models/export.jl sbml      $SCHEDULE $BUILD.sbml.tmp      $SEED && mv $BUILD.sbml.tmp $BUILD.sbml.xml; }
    [ -f $BUILD.gillespie.R ] || { julia --project=. models/export.jl gillespie $SCHEDULE $BUILD.gillespie.tmp $SEED && mv $BUILD.gillespie.tmp $BUILD.gillespie.R; }
}

run() {
    local out=$1 runner=$2; shift 2
    case $runner in
        grs) julia --project=. runners/grs.jl $out $SCHEDULE --seed $SEED --trajectories $N "$@" ;;
        g2) Rscript runners/gillespie.R $out $BUILD.gillespie.R --seed $SEED --trajectories $N --duration $DURATION --samples $SAMPLES --stop-on-neg FALSE "$@" ;;
        copasi) .venv/bin/python runners/copasi.py $out $BUILD.sbml.xml --seed $SEED --trajectories $N --duration $DURATION --samples $SAMPLES "$@" ;;
    esac
}

for k in "${KS[@]}"; do
    num_genes=$((2 ** k)); rung=$(printf "genes-%05d" $num_genes)
    SCHEDULE=$SCHEDULES/$MODEL/$rung.json
    BUILD=$BUILDS/$MODEL/$rung
    prepare $k

    for spec in "${ARMS[@]}"; do
        set -- $spec
        arm=$1 ceiling=$2 runner=$3; shift 3
        [ -z "${ONLY:-}" ] || [ $arm = $ONLY ] || continue
        [ $num_genes -le $ceiling ] || continue
        out=$RUNS/$MODEL/$rung/$arm-n$N-s$SEED
        PLANNED+=("$out")
        [ -f $out/metrics.csv ] && continue
        echo "=== [$(date -u +%H:%MZ)] $rung $arm"
        run $out $runner "$@"
    done
done

julia --project=. analysis/costs.jl $OUT "${PLANNED[@]}"
