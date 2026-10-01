#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

N=3
SEED=1
GENES=(10000)
TAU=0.0833333333333333
RUNS=runs
MODELS=models/dyngen
OUT=measurements/$(basename "$0" .sh).csv

if [ -n "${SMOKE:-}" ]; then
    TMP=$(mktemp -d); [ -n "${KEEP:-}" ] || trap 'rm -rf "$TMP"' EXIT
    echo "smoke output: $TMP" >&2
    RUNS=$TMP/runs; MODELS=$TMP/models; OUT=$TMP/out.csv
    GENES=(100)
fi

ARMS=(
  "dyngen-tauleap-tau0.0833    dyngen  --method etl --tau $TAU"
  "grs-tauleap-tau0.0833       grs     --method FixedTau --tau $TAU --policy AlwaysLeap --thin-negative false"
  "grs-tauleap-tau0.0833-thin  grs     --method FixedTau --tau $TAU --policy AlwaysLeap --thin-negative true"
  "grs-tauleap-eps0.05         grs     --method AdaptiveTau --epsilon 0.05"
  "grs-hybrid-tau0.0833-tausplitting  grs  --method FixedTau --tau $TAU --policy CriticalBlend --nc 10 --exact TauSplitting"
  "grs-tausplitting            grs     --method TauSplitting"
)

PLANNED=()
for genes in "${GENES[@]}"; do
    rung=$(printf "genes-%05d" $genes)
    MODEL=$MODELS/$rung
    [ -f $MODEL/schedule.json ] || Rscript models/dyngen.R --out $MODEL --genes $genes \
      --trajectories $N --seed $SEED --promoter thermodynamic

    for spec in "${ARMS[@]}"; do
        set -- $spec
        arm=$1 runner=$2; shift 2
        out=$RUNS/dyngen/$rung/$arm-n$N-s$SEED
        PLANNED+=("$out")
        [ -f $out/metrics.csv ] && continue
        echo "=== [$(date -u +%H:%MZ)] $rung $arm"
        case $runner in
            dyngen) Rscript runners/dyngen.R --out $out --model $MODEL/model.rds \
                      --trajectories $N --seed $SEED "$@" ;;
            grs) julia --project=. runners/grs.jl $out $MODEL/schedule.json \
                   --seed $SEED --trajectories $N "$@" ;;
        esac
    done
done

julia --project=. analysis/costs.jl $OUT "${PLANNED[@]}"

