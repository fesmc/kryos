#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=ba1442
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name=pm
# PaleoMIST test of the standalone VILMA2 (vilma.x of the S11 stack, vilma2 523b660) under the PaleoMIST ice history.
#   sbatch --export=ALL,RUN=B1_l256,LMAX=256,EARTH=1d,FORCE=B pm_run.sh   (the runs: -p shared -c 128 --mem=120G, 64 threads on 64 cores)
#   B1_L80: RUN=B1_L80,LMAX=170,EARTH=l80,FORCE=B,TMAX=2.0e6; B1_bm: RUN=B1_bm,LMAX=170,EARTH=bm,FORCE=B,TMAX=2.0e6 (after pm_visc1d.sh)
# Env: RUN (the folder name), LMAX, EARTH (1d: Gowan's three-layer Earth on PREM | l80: Gowan's viscosities, 80 km lithosphere |
#      bm: the science setting collapsed to 1-D | 3d: the science setting, Bagge 2021; see below),
#      FORCE (A: 22.5 ka ice held | B: 80 -> 0 ka steps | S: smoke, 25 -> 20 ka), TMAX (spin-up cap [yr], 5e5),
#      RATE_TOL (1e-5 m/yr), NTHR (OpenMP threads, 64).
#      pm-diag: FORCE also Bg | Bl | Bm | Bm910 | Bgm910 | Bm80 | Bgm80_910 (pm_prep.jl), TMAX=0 (no spin-up: the Earth relaxed to the present-day reference,
#      the 80 ka load switched on at 80 ka), EARTH=t3 (SELEN's Earth of PaleoMIST), ROT (.true.), DEG1 (cm | cf), RHOI (931), RHOW (1000),
#      REF (the reference file in $P/data, pm_ref | pm_ref80: SELEN's initial state, with TMAX=0), SUBGRID (sle_subgrid: .true. the
#      sloping-coast ocean load s = C rsl - z0 (C - C0) | .false. the binary C rsl, an SLE without the topography term),
#      ROTKS (rotation_k_s: 0 = the model's own fluid Love number; 0.94349 = the value from the observed flattening that vilma.x
#      prints, the revised rotation theory), SCHEME (the memory scheme, fe | trap), CFL (the Maxwell-number ceiling of the fe sub-steps, 1.0), ROOT (the parent folder of the
#      run, $P). The vilma2 clone of S11 is at c17d513, whose tree is that of 523b660 ($P/provenance.txt).
# Inputs from pm_prep.jl in $P/data; run folder $P/$RUN; the remap weights are built per run in $RUN/maps.
set -uo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/runs/paleomist; S=$G/stacks/S11
VX=$S/vilma2/bin/vilma.x
RUN=${RUN:?}; LMAX=${LMAX:?}; EARTH=${EARTH:?1d|l80|bm|3d|t3}; FORCE=${FORCE:?A|B|S|Bg|Bl|Bm|Bm910|Bgm910|Bm80|Bgm80_910}
TMAX=${TMAX:-5.0e5}; RATE_TOL=${RATE_TOL:-1.0e-5}
ROT=${ROT:-.true.}; DEG1=${DEG1:-cm}; RHOI=${RHOI:-931.0}; RHOW=${RHOW:-1000.0}; ROOT=${ROOT:-$P}
SCHEME=${SCHEME:-fe}; CFL=${CFL:-1.0}; REF=${REF:-pm_ref}; SUBGRID=${SUBGRID:-.true.}; ROTKS=${ROTKS:-0.0}
[ -f $P/data/pm_force_$FORCE.nc ] || { echo "missing $P/data/pm_force_$FORCE.nc (pm_prep.jl)"; exit 1; }
# The Earths: layers (Gowan's 12 | the science setting's 11, PREM rho and mu as r^2 shell means), eta per layer, lateral structure.
#   1d : Gowan et al. (2021): 120 km elastic lithosphere, upper mantle 4e20 Pa s to 670 km, lower mantle 4e22 Pa s; the layers split on the
#        PREM boundaries 80, 220, 400, 600, 670 km and the lower mantle in four
#   l80: Gowan's viscosities on the science layers (80 km elastic lithosphere): l80 - 1d = the lithosphere
#   bm : the science setting collapsed to 1-D: the profile of its pre_spinup_1d (lateral mean of log10 eta at the FE nodes, pm_visc1d.sh)
#        as a laterally uniform field ($P/data/pm_visc_bm.nc): bm - l80 = the radial profile, 3d - bm = the lateral structure
#   3d : the science setting (par/kryos_NH-ANT_vilma2.nml, as the G2b standalone runs): PREM, 80 km lithosphere, Bagge 2021 structure 6
#   t3 : SELEN's Earth of PaleoMIST (calc_sl/data.inc of the PANGAEA bundle: TABOO NV = 2, CDE = 0; model "ehgr": 120 km, UM 3.9811e20,
#        LM 3.9811e22 Pa s): three homogeneous layers 0-120, 120-670 km, 670 km-CMB and the core, rho and mu the volume means of TABOO's
#        200 s PREM table (pm_prep.jl taboo); ETA_UM, ETA_LM override the viscosities (the other Earths of Gowan's Stokes
#        coefficients, e.g. ehgG: LM 3.9811e21)
LAYERS_GOWAN=" n_layer = 12
 r_top = 6371.0e3, 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
 r_bot = 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
 rheology = 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2"
LAYERS_SCI=" n_layer = 11
 r_bot = 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
 r_top = 6371.0e3, 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
 rheology = 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2"
ETA_GOWAN_SCI=" eta = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 0.0"
ETA_SCI=" eta = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 0.0"
visc3d() { printf ' l_visc_3d = .true.\n visc_3d_file = "%s"\n visc_3d_structure = %s\n pre_spinup_1d = %s\n' "$1" "$2" "$3"; }
NO3D=" l_visc_3d = .false.
 pre_spinup_1d = .false."
case $EARTH in
  1d)  EARTH_NML="$LAYERS_GOWAN
 eta = 1.0e40, 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 0.0
$NO3D" ;;
  l80) EARTH_NML="$LAYERS_SCI
$ETA_GOWAN_SCI
$NO3D" ;;
  bm)  [ -f $P/data/pm_visc_bm.nc ] || { echo "missing $P/data/pm_visc_bm.nc (pm_visc1d.sh)"; exit 1; }
       EARTH_NML="$LAYERS_SCI
$ETA_SCI
$(visc3d $P/data/pm_visc_bm.nc 1 .false.)" ;;
  3d)  EARTH_NML="$LAYERS_SCI
$ETA_SCI
$(visc3d fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc 6 .true.)" ;;
  t3)  EARTH_NML=" n_layer = 4
 r_top = 6371.0e3, 6251.0e3, 5701.0e3, 3480.0e3
 r_bot = 6251.0e3, 5701.0e3, 3480.0e3, 0.0
 rho = 3233.60, 3630.57, 4877.91, 10931.73
 mu = 6.12924e10, 8.83454e10, 2.19481e11, 0.0
 eta = 1.0e40, ${ETA_UM:-3.9811e20}, ${ETA_LM:-3.9811e22}, 0.0
 rheology = 0, 1, 1, 2
$NO3D" ;;
  *)   echo "EARTH: 1d | l80 | bm | 3d | t3"; exit 1 ;;
esac

OUT=$ROOT/$RUN
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
mkdir -p $OUT/input $OUT/maps; cd $OUT
for f in $S/vilma2/input/*; do ln -s $f input/; done
ln -s /work/ba1442/fesmdata fesmdata

cat > pm.nml <<EOF
&ctl
 file_forcing = "$P/data/pm_force_$FORCE.nc"
 name_ice = "h_ice"
 name_time = "time"
 file_ref = ""
 name_zbed_eq = "z_bed_eq"
 name_hice_ref = "h_ice_ref"
 file_out = "vilma_out.nc"
 file_hor = ""
 write_lossless = .true.
 time_init = -1.0e30
 time_end = 1.0e30
 remap_input = .true.
 name_lon = "lon"
 name_lat = "lat"
 i_eq = 1
 z_bed_ref_file = "$P/data/$REF.nc"
 h_ice_ref_file = "$P/data/$REF.nc"
 name_z_bed_ref = "z_bed"
 name_h_ice_ref = "H_ice"
 z_bed_eq_file = ""
 h_ice_eq_file = ""
 rsl_restart_file = ""
 name_rsl = "rsl"
 restart_in_file = ""
/
&vilma
 lmax = $LMAX
 scheme = "$SCHEME"
 sle_subgrid = $SUBGRID
 cfl = $CFL
 rotation = $ROT
 rotation_k_s = $ROTKS
 deg1_frame = "$DEG1"
 rho_ice = $RHOI
 rho_water = $RHOW
 earth = "$([ $EARTH = t3 ] && echo custom || echo PREM)"
$EARTH_NML
 equil_rate_tol = $RATE_TOL
 equil_time_max = $TMAX
/
EOF

export OMP_NUM_THREADS=${NTHR:-64} OMP_PROC_BIND=spread OMP_PLACES=cores OMP_STACKSIZE=512M
ulimit -s unlimited
echo "== vilma.x $RUN lmax $LMAX earth $EARTH force $FORCE tmax $TMAX tol $RATE_TOL rot $ROT deg1 $DEG1 rho $RHOI/$RHOW scheme $SCHEME cfl $CFL ref $REF subgrid $SUBGRID k_s $ROTKS $(date -Is) on $(hostname)"
/usr/bin/time -v $VX pm.nml > run.log 2>&1; echo "vilma.x exit $?"
grep -E "spin-up|WARNING|visc3d split|PROFILE setup|n_accept|Maximum resident|Elapsed" run.log
echo "== done $(date -Is)"
