#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=00:15:00 --mem=4G --job-name=g5-gen -o /work/ba1442/robinson/kryos-gate/s11/g5gen/gen.%j.log
# G5 pars of the calib-cf science pars (docs/dev/calib-cf.md) on Levante: copy tools/bin/kryos-ant-par.jl and par/kryos_NH-ANT{,_vilma2}_calib.nml
# of the checkout to $G/s11/g5gen (same layout), sbatch g5_gen.sh, copy back par/gate/g5-{vilma2,fastiso}-calib.nml.
cd /work/ba1442/robinson/kryos-gate/s11/g5gen; J="$HOME/.juliaup/bin/julia +1.13 --startup-file=no"
mkdir -p par/gate
$J tools/bin/kryos-ant-par.jl par/kryos_NH-ANT_vilma2_calib.nml par/gate/g5-vilma2-calib.nml
$J tools/bin/kryos-ant-par.jl par/kryos_NH-ANT_calib.nml par/gate/g5-fastiso-calib.nml
ls -l par/gate
