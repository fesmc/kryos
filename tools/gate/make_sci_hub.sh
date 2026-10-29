#!/bin/bash
# par/gate/sci-hub.nml = par/gate/sci.nml (FastIsostasy, two domains) + the global geometry hub as an observer (gate G1b):
# [globe_ice] on GLOBAL-0.5DEG with the PaleoMIST history, target grid GAUSS-L170 (the grid VILMA2 would get), globe_ice.nc at out_2D.
#   bash tools/gate/make_sci_hub.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
{
echo '! par/gate/sci.nml with the global geometry hub as an observer (gate G1b, tools/gate/make_sci_hub.sh): [globe_ice] and its datasrc groups appended.'
echo '! The hub gathers the domains'"'"' ice at every step and writes globe_ice.nc at out_2D; it must not change any result of sci.nml.'
echo
cat par/gate/sci.nml
cat <<'BLOCK'

&globe_ice
    grid      = "GLOBAL-0.5DEG"      ! background grid
    z_bed_ref = "X_z_bed"
    H_ice_ref = "X_H_ice"
    target    = "GAUSS-L170"         ! the global grid of the host terms (observer: no VILMA2)
    ext       = "file"               ! X + the history's anomaly (the run is after the history's last record: it is held, anomaly 0)
    ext_field = "X_dH_ice"
    ext_tref  = 0.0
    out_2D    = "out_2D"             ! globe_ice.nc every 100 yr
/

&X_z_bed
    filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"
    name     = "z_bed"
    mode     = "transient"
/

&X_H_ice
    filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"
    name     = "H_ice"
    mode     = "transient"
/

&X_dH_ice
    filename     = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"
    name         = "H_ice"
    mode         = "transient"
    anomaly_ref  = "slice"
    ref_time     = 0.0
    out_of_range = "hold"
/
BLOCK
} > par/gate/sci-hub.nml
