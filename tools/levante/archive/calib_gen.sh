#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=00:15:00 --mem=4G --job-name=calib-gen -o /work/ba1442/robinson/kryos-gate/calib/gen/gen.%j.log
# calib-cf pars (docs/dev/calib-cf.md) on Levante: copy tools/bin/kryos-calib-par.jl and par/kryos_NH-ANT{,_vilma2}.nml of the checkout to $G/calib/gen
# (same layout), sbatch calib_gen.sh, copy back par/kryos_NH-ANT*_calib.nml and par/calib/*.nml.
cd /work/ba1442/robinson/kryos-gate/calib/gen; J="$HOME/.juliaup/bin/julia +1.13 --startup-file=no"
$J tools/bin/kryos-calib-par.jl sci par/kryos_NH-ANT.nml par/kryos_NH-ANT_calib.nml
$J tools/bin/kryos-calib-par.jl sci par/kryos_NH-ANT_vilma2.nml par/kryos_NH-ANT_vilma2_calib.nml
for d in Antarctica North; do for c in pd lgm; do $J tools/bin/kryos-calib-par.jl run par/kryos_NH-ANT_calib.nml par/calib/$d-$c.nml $d $c; done; done
$J tools/bin/kryos-calib-par.jl run par/kryos_NH-ANT_calib.nml par/calib/North-lgmrecon.nml North lgmrecon
ls -l par par/calib
