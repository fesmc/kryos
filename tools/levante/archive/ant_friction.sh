#!/bin/bash -l
# ant-friction (docs/dev/ant-friction.md): ANT PD runs of the friction sweep on stack S9 via calib.sh (calib-cf pd case, 30 kyr cold start).
#   ant_friction.sh <cf_ref> <z0> <z1> ...      (triples; TAG -af<cf_ref without "0.">-z<-z0>_<-z1>, e.g. -af05-z400_0)
# Fixed on top of $G/calib/par/Antarctica-pd.nml: enh_shear 3, scale_zb 2 (exp, as the file), vm-m16 calving with tau_ice_flt 85 kPa,
# no grounded calving, H_min 0 (calib-cf round 4 `-cf05-vm`).
G=/work/ba1442/robinson/kryos-gate
while [ $# -ge 3 ]; do
  cf=$1; z0=$2; z1=$3; shift 3
  TAG=-af${cf#0.}-z${z0#-}_${z1#-} $G/calib/calib.sh Antarctica pd ymat_ant.enh_shear=3.0 \
    ytill_ant.cf_ref=$cf ytill_ant.z0=$z0 ytill_ant.z1=$z1 \
    'ycalv_ant.calv_flt_method="vm-m16"' 'ycalv_ant.calv_grnd_method="zero"' ycalv_ant.H_min_grnd=0.0 ycalv_ant.H_min_flt=0.0 \
    ycalv_ant.tau_ice_flt=85.0e3 ycalv_ant.tau_ice_grnd=85.0e3
done
