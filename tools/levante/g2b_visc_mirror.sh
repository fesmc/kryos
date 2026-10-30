#!/bin/bash
# G2b diagnosis: a fesmdata view whose Bagge 3D viscosity is mirrored north-south (data reversed along lat, coordinates kept),
# everything else linked to /work/ba1442/fesmdata. A kryos run with FESM=<this folder> (g2b.sh) puts the viscosity under the
# load where VILMA2 sees it when the host fields enter unflipped (ascending host lat into the north-first Gauss rows).
#   bash g2b_visc_mirror.sh [out folder ($G/p4g2b/diag/fesmdata-mirror)]
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; D=${1:-$G/p4g2b/diag/fesmdata-mirror}; S=/work/ba1442/fesmdata
V=GLOBAL-T170_SolidEarth-Bagge2021.nc
mkdir -p $D/Global/GLOBAL-T170
for e in $S/*; do [ "$(basename $e)" = Global ] || ln -sfn $e $D/; done
for e in $S/Global/*; do [ "$(basename $e)" = GLOBAL-T170 ] || ln -sfn $e $D/Global/; done
for e in $S/Global/GLOBAL-T170/*; do [ "$(basename $e)" = $V ] || ln -sfn $e $D/Global/GLOBAL-T170/; done
F=$D/Global/GLOBAL-T170/$V
cp $S/Global/GLOBAL-T170/$V $F; chmod u+w $F
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot
julia +1.13 --project=$G/p4w7/kryos/tools -e "using NCDatasets; ds = NCDataset(\"$F\", \"a\"); v = ds[\"log10_visc\"]; a = v.var[:, :, :, :]; v.var[:, :, :, :] = a[:, end:-1:1, :, :]; f = ds[\"f_valid\"]; b = f.var[:, :]; f.var[:, :] = b[:, end:-1:1]; ds.attrib[\"g2b_note\"] = \"data reversed along lat (north-south mirror) for the G2b diagnosis\"; close(ds); println(\"mirrored \", size(a))"
