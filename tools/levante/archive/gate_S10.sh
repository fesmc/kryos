#!/bin/bash -l
# The S10 runs on stack S10 (setup_S10.sh), kryos at the S10 checkouts. Prints the job ids. KSHA: the commit $S/kryos is at.
#   gate_S10.sh g0    G0: make check and make check-long of $S/kryos (vilma=1), check of $S/kryos-v0 (check_kryos.sh)
#   gate_S10.sh g1    G1: K(S10, vilma=1) cont for grl, ant, ant-full, nor (maps S5) and the two-domain sci, sci-const cont
#   gate_S10.sh g4    G4: sci-vilma2 (cont, cont2, chunk at -79800, lossless; g34_vilma2.sh) and the sci chunk (fastiso, 0 -> 500 -> 1000)
#   gate_S10.sh bc    backward compatibility (after g4): c01 of sci-vilma2 (-79800 -> -79600) restored by kryos-x (+ diagnostic) from
#                     the S9 bundle and by kryos-d (- diagnostic) from the S10 bundle; the run folders of tools/gate/sci-vilma2-bc.toml
#                     are assembled in $G/runs/sci-vilma2-bc/kx-$KSHA/S10 (cont, bc-*/c00: links to the S10 runs).
# Compare with compare_S10.sh g1|g4|bc.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S10; P=$G/p4w7
KSHA=$(git -C $S/kryos rev-parse --short=7 HEAD)
RK=${RK:-$P/levante/run_kryos.sh}; RKK=${RKK:-$P/levante/run_kryos_kk.sh}
RV=$G/rst10/run_vilma2.sh   # the repo's (SINGLE=1 runs, WP7c); $P/levante has the older one
M=${MAPS:-$G/maps/S5}
OV=(clock.time_end=-79600 out_2D.dt=50 out_2Dsm.dt=50 out_rst.dt=200 'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"')
case ${1:-} in
  g0)
    sbatch --parsable --export=ALL,STK=S10,V=1,TGT="check check-long",LOGTAG=vilma1 $G/rst10/check_kryos.sh
    sbatch --parsable --export=ALL,STK=S10,V=0,LOGTAG=vilma0 $G/rst10/check_kryos.sh ;;
  g1)
    for g in grl ant ant-full nor; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S10 MAPS=$M $RK $g cont
    done
    for g in sci sci-const; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S10 MAPS=$M TOOLS=$P/kryos $RKK $g cont
    done ;;
  g4)
    STK=S10 K=$S/kryos $P/levante/g34_vilma2.sh
    K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S10 MAPS=$M TOOLS=$P/kryos $RKK sci chunk ;;
  bc)
    B9=$G/runs/sci-vilma2/kx-624f686/S9/cont/restart--79800
    B10=$G/runs/sci-vilma2/kx-$KSHA/S10/cont/restart--79800
    for b in $B9 $B10; do [ -f $b/manifest.nml ] || { echo "no bundle $b"; exit 1; }; done
    BC=(clock.time_init=-79800 clock.time_end=-79600)
    STK=S10 K=$S/kryos-x JOBTIME=00:45:00 SINGLE=1 $RV sci-vilma2 bc-s9x "${OV[@]}" "${BC[@]}" "kryos.restart=\"$B9\""
    STK=S10 K=$S/kryos-d JOBTIME=00:45:00 SINGLE=1 $RV sci-vilma2 bc-drop "${OV[@]}" "${BC[@]}" "kryos.restart=\"$B10\""
    A=$G/runs/sci-vilma2-bc/kx-$KSHA/S10; mkdir -p $A
    ln -sfn $G/runs/sci-vilma2/kx-$KSHA/S10/cont $A/cont
    for r in bc-s9x:kryos-x bc-drop:kryos-d; do
        n=${r%%:*}; sx=$(git -C $S/${r##*:} rev-parse --short=7 HEAD)
        mkdir -p $A/$n
        ln -sfn $G/runs/sci-vilma2/kx-$KSHA/S10/chunk/c00 $A/$n/c00
        ln -sfn $G/runs/sci-vilma2/kx-$sx/S10/$n $A/$n/c01
    done ;;
  *) echo "usage: gate_S10.sh g0|g1|g4|bc"; exit 1 ;;
esac
