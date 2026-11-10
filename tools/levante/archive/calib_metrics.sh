#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=01:00:00 --mem=40G --job-name=jl-calib -o /work/ba1442/robinson/kryos-gate/calib/out/jl_calib.%j.log
# calib-cf metrics of every calib run (docs/dev/calib-cf.md): sbatch calib_metrics.sh [run folder ...] (default: all of $ROOT).
# Run folder <Domain>-<case><TAG> -> out/<ant|nor>-<case><TAG>; lgmrecon is measured against the LGM references.
G=/work/ba1442/robinson/kryos-gate; ROOT=$G/runs/calib/kx-624f686/S9
export JULIA_DEPOT_PATH=$G/julia-depot; J="$HOME/.juliaup/bin/julia +1.13 --project=$G/p4w7/kryos/tools"
cd $G/calib
runs=("$@"); [ ${#runs[@]} -eq 0 ] && runs=($(cd $ROOT && ls -d Antarctica-* North-* | grep -v _test))
for r in "${runs[@]}"; do
    dom=${r%%-*}; rest=${r#*-}; case=${rest%%-*}; [ $case = lgmrecon ] && mcase=lgm || mcase=$case
    pfx=$([ $dom = Antarctica ] && echo ant || echo nor)-$rest
    echo "== $r -> out/$pfx"
    $J calib_metrics.jl $ROOT/$r $dom $mcase out/$pfx
done
