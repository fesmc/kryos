#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=ba1442
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name=pm
# PaleoMIST test of the standalone VILMA2 (vilma.x of the S11 stack, vilma2 523b660) under the PaleoMIST ice history.
#   sbatch --export=ALL,RUN=B1_l256,LMAX=256,EARTH=1d,FORCE=B pm_run.sh   (the runs: -p shared -c 128 --mem=120G, 64 threads on 64 cores)
# Env: RUN (the folder name), LMAX, EARTH (1d: Gowan's three-layer Earth on PREM | 3d: the science setting, Bagge 2021),
#      FORCE (A: 22.5 ka ice held | B: 80 -> 0 ka steps | S: smoke, 25 -> 20 ka), TMAX (spin-up cap [yr], 5e5),
#      RATE_TOL (1e-5 m/yr), NTHR (OpenMP threads, 64). The vilma2 clone of S11 is at c17d513, whose tree is that of 523b660 ($P/provenance.txt).
# Inputs from pm_prep.jl in $P/data; run folder $P/$RUN; the remap weights are built per run in $RUN/maps.
set -uo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/runs/paleomist; S=$G/stacks/S11
VX=$S/vilma2/bin/vilma.x
RUN=${RUN:?}; LMAX=${LMAX:?}; EARTH=${EARTH:?1d|3d}; FORCE=${FORCE:?A|B|S}
TMAX=${TMAX:-5.0e5}; RATE_TOL=${RATE_TOL:-1.0e-5}
OUT=$P/$RUN
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
mkdir -p $OUT/input $OUT/maps; cd $OUT
for f in $S/vilma2/input/*; do ln -s $f input/; done
ln -s /work/ba1442/fesmdata fesmdata

if [ $EARTH = 1d ]; then
  # Gowan et al. (2021): 120 km elastic lithosphere, upper mantle 4e20 Pa s to 670 km, lower mantle 4e22 Pa s; PREM rho and mu
  # (r^2 shell means), the layers split on the PREM boundaries 80, 220, 400, 600, 670 km and the lower mantle in four.
  EARTH_NML=$(cat <<EOF
 earth = "PREM"
 n_layer = 12
 r_top = 6371.0e3, 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
 r_bot = 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
 eta = 1.0e40, 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 0.0
 rheology = 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
 l_visc_3d = .false.
 pre_spinup_1d = .false.
EOF
)
else
  # the science setting (par/kryos_NH-ANT_vilma2.nml, as the G2b standalone runs): PREM, 80 km lithosphere, Bagge 2021 structure 6
  EARTH_NML=$(cat <<EOF
 earth = "PREM"
 n_layer = 11
 r_bot = 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
 r_top = 6371.0e3, 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
 eta = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 0.0
 rheology = 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
 l_visc_3d = .true.
 visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"
 visc_3d_structure = 6
 pre_spinup_1d = .true.
EOF
)
fi

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
 z_bed_ref_file = "$P/data/pm_ref.nc"
 h_ice_ref_file = "$P/data/pm_ref.nc"
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
 rotation = .true.
$EARTH_NML
 equil_rate_tol = $RATE_TOL
 equil_time_max = $TMAX
/
EOF

export OMP_NUM_THREADS=${NTHR:-64} OMP_PROC_BIND=spread OMP_PLACES=cores OMP_STACKSIZE=512M
ulimit -s unlimited
echo "== vilma.x $RUN lmax $LMAX earth $EARTH force $FORCE tmax $TMAX tol $RATE_TOL $(date -Is) on $(hostname)"
/usr/bin/time -v $VX pm.nml > run.log 2>&1; echo "vilma.x exit $?"
grep -E "spin-up|WARNING|visc3d split|PROFILE setup|n_accept|Maximum resident|Elapsed" run.log
echo "== done $(date -Is)"
