# WP0 data check (Q17): files each new base opens against the gate's links

Date 2026-10-28. Gate root `$G = /work/ba1442/robinson/kryos-gate` on Levante. kryos main `6c731ea`; yelmox `topo-prep` `5ac89c8` (R2': `[topo]` keys in all bases).

## Result

- **Every file a base reads exists** (0 missing among the keys that are read; GRP, LIS, NOR, PYR, SRG, eANT, eGRL, rGRL, and the shipped ANT extras).
- **No link was made.** `fesmdata`, `ice_data` and `isostasy_data` of the gate (stack roots and `$G/kryos`) are symlinks to `/work/ba1442/<name>`, so everything in the shared store is already visible; nothing is missing there. The rule "never write into a link folder that is itself a symlink into a shared store" leaves nothing to do. `maps` (`$G/maps/S`) is frozen read-only and belongs to the mapgen agent; `input` is the overlay of the stack's yelmox `input/` and kryos `input/`.
- Missing **maps** (not generated here): per base below. Only 5 maps exist in `$G/maps/S`; all grid descriptions needed (`grid_*.txt`) are there.
- The only gaps are in `input/` of the **kryos repo** (3 series files that only the yelmox `input/` has) and some quirks listed under "For Alex".

| base | data files | maps missing (src -> dst) | input/ series |
|---|---|---|---|
| GRP | ready | none | `sealevel_spratt2016.dat` (bsl `file`): stack only |
| LIS | ready | GLOBAL-1DEG -> LIS-32KM (bilin) | all in kryos `input/` |
| NOR | ready | GLOBAL-1DEG -> NH-32KM (bilin) | all in kryos `input/` |
| PYR | ready | PYR-250M <-> PYR-2KM (bilin, con, nn 250M->2KM) | `alpha_orbital_t2018.dat`: stack only |
| SRG | ready | SRG-250M <-> SRG-16KM (isos grid; yelmox needs con 250M->16KM, see maps) | none |
| eANT | ready (spinup reference set) | none (all ANT-8KM) | `esm/esm_ant.nml` is stack-only (the groups go into the kryos par file) |
| eGRL | ready (spinup reference set) | GRL-16KM -> GRL-8KM (bilin, esm ocean reference) | `esm/esm_grl_ismip7.nml`: as eANT |
| rGRL | ready | none (all GRL-16KM) | `timeout_ramp_100kyr.txt` (`tm_1D.method = "file"`): stack only |
| ANT-full | ready (no new files vs phase 1) | none | none |

## Method

1. **Resolved configuration.** `KryosTools.Convert.resolve_yelmox` (FesmUtils.jl `Namelist.resolve`, the converter's reference table `tools/data/yelmox2kryos.toml`) on the yelmox worktree `~/models/yelmox/.claude/worktrees/topo-prep`, one run per base with the phase-1 `-p` (`ctrl.time_end=1000 ctrl.write_lossless=True tm_2D.dt=100 tm_2Dsm.dt=10 tm_1D.dt=10 tm_rst.dt=500`) and the gate extras of phase2-plan §3.1 (PYR `ctrl.time_end=-119e3`; SRG `ctrl.time_end=200 tm_rst.dt=100`; eANT/eGRL `spinup.time_end=1000`, because `ctrl.time_end` does not exist with `run_step = "spinup"`; ANT-full without `sim.relax=False sim.opt=False`). Feature gates (grl `sim.use_negis`, `scale_glacial_smb`, `init_marine_H`, `init_kill_shelves`, `init_time_thrm`; grp `ice_bnd.active`) were resolved too: only grp `ice_bnd` reads an extra file (PaleoMIST, present); the others read none.
2. **Every path/file key** of the resolved groups was collected (string values with `/` or a file extension), then each was classified by **reading the Fortran** (yelmox libs at `5ac89c8`, Yelmo and FastIsostasy of stack S3, fesm-utils `regions`/`varslice`, rembo1 v1.3 of S3) for the condition under which it is opened and for the grid its `{grid_name}` takes (hub: topo, regions; ice: Yelmo data, sediments, geothermal, recon; isos: rheology/mask; clim: snapshots and esm; surface: smb_simple mask; shelf: marine_shelf). `{domain}`, `{set}`, `{grid_global}`, `{exp}`, `{gcm}`, `{experiment}` were substituted as yelmox does.
3. **By hand (the resolver does not follow them):** groups named by keys outside the converter's reference rules: `[bias] dtas/fpr` (GRP), `[ice_bnd] H_ice`, the esm `par_file` (`input/esm/*.nml`, groups `gcm_*`; the gate runs the `spinup` clock with `use_esm = use_var = False`, so only the reference groups are read: `ts_ref`, `pr_ref` or `smb_ref` (eGRL, `surface_method = climate`), `zs_ref`, `to_ref`, `so_ref`), and rembo1's own par file. The resolver also keeps groups that the par file defines but nothing references (e.g. `[coef_alpha]` in bases with `coefs = ""`); those were dropped by evaluating `forcing.coefs`, `terms`, `reference`, `snap0/snap1`.
4. **Existence** through the gate's own links (`$G/kryos/{fesmdata,ice_data,isostasy_data}`); `input/` against the stack's yelmox `input/` (`$G/stacks/S3/yelmox/input`, what the gate overlay links) and against kryos main `input/`. Then `ncdump -h` of every existing netCDF file (89): the configured variable names exist (topo, regions layer, yelmo_data, varslice `name`, sed, ghf, isos rheology, esm, rembo, recon `dz`, smb_simple `mask`; 0 missing), the `grid_name` attribute equals the expected grid, and files without the attribute have the `xc`/`yc` size of the expected grid (all do). Region names used in masks exist in the hub `REGIONS.nc` flag_meanings (Greenland, Canadian_Arctic_Archipelago, North_America, Hudson_Bay_and_Strait, Antarctic_Peninsula, East_Antarctica, West_Antarctica, Antarctica; zones open_ocean, shelf_break_buffer).
5. **Maps** by reading the `remap` calls (`libs/kryos.f90`, `kryos_coupling.f90`, `kryos_startup.f90`, `regions_init_nml`, varslice `grid_src`), then compared with what the yelmox mapgen runs have written so far in `$G/maps/gen/*-maps` (read only, 23:59 CET; the mapgen agent is still running).

**Blind spots.** (a) Only the yelmox readers were followed; a kryos backend that is not written yet (esm, rembo, smb_simple, `data`) may open the same files differently or choose another grid for `{grid_name}`. (b) The time coverage of each varslice (`time = ...` of the group against the file's time axis), units and fill values were not checked; only existence, variable names, dimensions and the grid attribute. (c) eANT/eGRL `hist`/`proj`/`esm_ref` groups are not in the gate (clock `spinup`); they were only looked up by path (see below). (d) The map list is derived from code reading and the mapgen output so far, not from a run. (e) Conditions rest on the resolved values (e.g. `marine_shelf.use_obs = False` in all bases, so `bmb` observation files are not read); a `-p` that flips such a key changes the list.

## Per-base tables

"expected" = where the plan says the product lives (grid/topo in fesmdata, esm data in ice_data). Deviations: the legacy products (`SED-M23`, `GEO-P22`, `ICE-6G_C`, Batchelor masks, `BIAS-MARv3.14-ERA5`, RACMO23, ORAS5 for the esm ocean) and the ISMIP7 files are in `ice_data` in every case, as yelmox expects. "action none" = nothing to link or copy.

### GRP (GreenlandPaleo)

`yelmox_GreenlandPaleo.nml`, domain `GreenlandPaleo`; grids: hub GRL-PAL-16KM, ice GRL-PAL-16KM, isos GRL-PAL-16KM, clim GRL-PAL-16KM, surface GRL-PAL-16KM, shelf GRL-PAL-16KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_TOPO-BedMachine-v6.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_basins` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BASINS-Zwally2012.nc` | fesmdata | yes | none | basin set Zwally2012 |
| `regions.path_basins` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BASINS-NEGIS.nc` | fesmdata | yes | none | basin set NEGIS |
| `yelmo_data.pd_age_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Stratigraphy-MacGregor2025.nc` | fesmdata | yes | none | pd_age_load = True |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = file |
| `barysealevel.sl_path` | `input/sealevel_spratt2016.dat` | input | yes (stack `input/`) | none (gate overlay supplies it) | method = file |
| `ghf.obs_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_GHF-FoxMaule2005.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `era5_1981_2010_tas.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) tas |
| `era5_1981_2010_pr.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) pr |
| `era5_1981_2010_zs.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) zs |
| `oras5_1981_2010_to.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) to |
| `oras5_1981_2010_so.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) so |
| `bias_mar_era5_1981_2010_dtas.filename` | `ice_data/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BIAS-MARv3.14-ERA5-1981-2010.nc` | ice_data | yes | none | `[bias] active = True`; group of the library, not followed by the converter's resolver |
| `bias_mar_era5_1981_2010_fpr.filename` | `ice_data/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BIAS-MARv3.14-ERA5-1981-2010.nc` | ice_data | yes | none | `[bias] active = True`; group of the library, not followed by the converter's resolver |
| `paleomist_H_ice.filename` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_IceHistory-PaleoMIST.nc` | fesmdata | yes | none | feature gate `ice_bnd.active=True` only (WP5) |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the present-day topography |
| `yelmo_data.pd_tsrf_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | no | pd_tsrf_load = False |
| `yelmo_data.pd_smb_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | no | pd_smb_load = False |
| `yelmo_data.pd_vel_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Velocity-Joughin2018.nc` | no | pd_vel_load = False |
| `sed.obs_path` | `fesmdata/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_Sediment-CRUST1.nc` | yes | sed.method = 0 |
| `marine_shelf.obs_path` | `ice_data/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BMELT-R13.nc` | no | use_obs = False |

### LIS

`yelmox_LIS.nml`, domain `Laurentide`; grids: hub LIS-32KM, ice LIS-32KM, isos LIS-32KM, clim LIS-32KM, surface LIS-32KM, shelf LIS-32KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_TOPO-GEBCO2026.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `isos.rheology_file` | `ice_data/Laurentide/LIS-32KM/LIS-32KM_GEO-P22.nc` | ice_data | yes | none | lithosphere/mantle = rheology_file |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `sed.obs_path` | `ice_data/Laurentide/LIS-32KM/LIS-32KM_SED-M23.nc` | ice_data | yes | none | sed.method = 1 |
| `ghf.obs_path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_GHF-Shapiro2004.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `sim.recon_path` | `ice_data/Laurentide/LIS-32KM/LIS-32KM_TOPO-ICE-6G_C.nc` | ice_data | yes | none | init_method = recon |
| `smb_simple.mask_file` | `ice_data/Laurentide/LIS-32KM/LIS-32KM_Batchelor2019_ice_masks.nc` | ice_data | yes | none | surface_method = smb_simple |
| `era5_1981_2010_tas.filename` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) tas |
| `era5_1981_2010_pr.filename` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) pr |
| `era5_1981_2010_zs.filename` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) zs |
| `oras5_1981_2010_to.filename` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) to |
| `oras5_1981_2010_so.filename` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) so |
| `pmip3_mean_tas_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_jja.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_pr_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol pr; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_zs.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol zs; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_jja.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_pr_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm pr; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_zs.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm zs; grid_src=GLOBAL-1DEG remap=bilin |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `yelmo_data.pd_tsrf_path` | `fesmdata/Laurentide/LIS-32KM/LIS-32KM_Atmosphere-ERA5-1981-2010.nc` | yes | pd_tsrf_load = False |
| `marine_shelf.obs_path` | `ice_data/Laurentide/LIS-32KM/LIS-32KM_BMELT-R13.nc` | no | use_obs = False |
| `barysealevel.sl_path` | `input/sl_deglaciation.dat` | yes (kryos `input/`) | method = fastiso |

### NOR (North)

`yelmox_North.nml`, domain `North`; grids: hub NH-32KM, ice NH-32KM, isos NH-32KM, clim NH-32KM, surface NH-32KM, shelf NH-32KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/North/NH-32KM/NH-32KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `sed.obs_path` | `fesmdata/North/NH-32KM/NH-32KM_Sediment-CRUST1.nc` | fesmdata | yes | none | sed.method = 1 |
| `ghf.obs_path` | `fesmdata/North/NH-32KM/NH-32KM_GHF-Shapiro2004.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `sim.recon_path` | `ice_data/North/NH-32KM/NH-32KM_TOPO-ICE-6G_C.nc` | ice_data | yes | none | init_method = recon |
| `era5_1981_2010_tas.filename` | `fesmdata/North/NH-32KM/NH-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) tas |
| `era5_1981_2010_pr.filename` | `fesmdata/North/NH-32KM/NH-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) pr |
| `era5_1981_2010_zs.filename` | `fesmdata/North/NH-32KM/NH-32KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) zs |
| `oras5_1981_2010_to.filename` | `fesmdata/North/NH-32KM/NH-32KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) to |
| `oras5_1981_2010_so.filename` | `fesmdata/North/NH-32KM/NH-32KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) so |
| `pmip3_mean_tas_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_jja.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_pr_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol pr; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_zs.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-piControl-mean.nc` | fesmdata | yes | none | term_lgm snapshot_picontrol zs; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_tas_jja.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm tas; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_pr_ann.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm pr; grid_src=GLOBAL-1DEG remap=bilin |
| `pmip3_mean_zs.filename` | `fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-PMIP3-lgm-mean.nc` | fesmdata | yes | none | term_lgm snapshot_lgm zs; grid_src=GLOBAL-1DEG remap=bilin |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/North/NH-32KM/NH-32KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/North/NH-32KM/NH-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/North/NH-32KM/NH-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `yelmo_data.pd_tsrf_path` | `fesmdata/North/NH-32KM/NH-32KM_Atmosphere-ERA5-1981-2010.nc` | yes | pd_tsrf_load = False |
| `marine_shelf.obs_path` | `ice_data/North/NH-32KM/NH-32KM_BMELT-R13.nc` | no | use_obs = False |
| `barysealevel.sl_path` | `input/sl_deglaciation.dat` | yes (kryos `input/`) | method = fastiso |

### PYR (Pyrenees)

`yelmox_Pyrenees.nml`, domain `Pyrenees`; grids: hub PYR-250M, ice PYR-250M, isos PYR-2KM, clim PYR-2KM, surface PYR-250M, shelf PYR-250M.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/Pyrenees/PYR-250M/PYR-250M_TOPO-GLO90-GEBCO2026.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Pyrenees/PYR-250M/PYR-250M_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `ghf.obs_path` | `fesmdata/Pyrenees/PYR-2KM/PYR-2KM_GHF-Lucazeau2019.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `coef_alpha.series_file` | `input/alpha_orbital_t2018.dat` | input | yes (stack `input/`) | none (gate overlay supplies it) | coef_alpha method = series |
| `era5_1981_2010_tas.filename` | `fesmdata/Pyrenees/PYR-2KM/PYR-2KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) tas |
| `era5_1981_2010_pr.filename` | `fesmdata/Pyrenees/PYR-2KM/PYR-2KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) pr |
| `era5_1981_2010_zs.filename` | `fesmdata/Pyrenees/PYR-2KM/PYR-2KM_Atmosphere-ERA5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) zs |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Pyrenees/PYR-250M/PYR-250M_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Pyrenees/PYR-250M/PYR-250M_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Pyrenees/PYR-250M/PYR-250M_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `barysealevel.sl_path` | `input/sl_deglaciation.dat` | yes (kryos `input/`) | method = fastiso |

### SRG

`yelmox_SRG.nml`, domain `SRG`; grids: hub SRG-250M, ice SRG-250M, isos SRG-16KM, clim SRG-250M, surface SRG-250M, shelf SRG-250M.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/SRG/SRG-250M/SRG-250M_TOPO-GLO90-GEBCO2026.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/SRG/SRG-250M/SRG-250M_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_basins` | `fesmdata/SRG/SRG-250M/SRG-250M_BASINS-SanRafael.nc` | fesmdata | yes | none | basin set SanRafael |
| `srg_mar314_1960_1990_tas.filename` | `ice_data/SRG/SRG-250M/SRG-250M_CLIM_PD-1960-1990_MAR3V14.nc` | ice_data | yes | none | snapshot_pd (reference) tas |
| `srg_mar314_1960_1990_pr.filename` | `ice_data/SRG/SRG-250M/SRG-250M_CLIM_PD-1960-1990_MAR3V14.nc` | ice_data | yes | none | snapshot_pd (reference) pr |
| `srg_mar314_1960_1990_zs.filename` | `ice_data/SRG/SRG-250M/SRG-250M_CLIM_PD-1960-1990_MAR3V14.nc` | ice_data | yes | none | snapshot_pd (reference) zs |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/SRG/SRG-250M/SRG-250M_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/SRG/SRG-250M/SRG-250M_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/SRG/SRG-250M/SRG-250M_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `yelmo_data.pd_tsrf_path` | `ice_data/SRG/SRG-250M/SRG-250M_CLIM_PD-1960-1990_MAR3V14.nc` | yes | pd_tsrf_load = False |
| `sed.obs_path` | `ice_data/SRG/SRG-250M/SRG-250M_SED-L97.nc` | no | sed.method = 0 |
| `ghf.obs_path` | `ice_data/SRG/SRG-250M/SRG-250M_GHF-S04.nc` | no | ghf.use_obs = False |
| `marine_shelf.obs_path` | `ice_data/SRG/SRG-250M/SRG-250M_BMELT-R13.nc` | no | use_obs = False |

### eANT (esm_Antarctica)

`yelmox_esm_Antarctica.nml`, domain `Antarctica`; grids: hub ANT-8KM, ice ANT-8KM, isos ANT-8KM, clim ANT-8KM, surface ANT-8KM, shelf ANT-8KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_TOPO-BedMachine-v4.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_basins` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_BASINS-Zwally2012.nc` | fesmdata | yes | none | basin set Zwally2012 |
| `yelmo_data.pd_tsrf_path` | `ice_data/Antarctica/ANT-8KM/RACMO2.3/ANT-8KM_RACMO23-VW23.nc` | ice_data | yes | none | pd_tsrf_load = True |
| `yelmo_data.pd_smb_path` | `ice_data/Antarctica/ANT-8KM/RACMO2.3/ANT-8KM_RACMO23-VW23.nc` | ice_data | yes | none | pd_smb_load = True |
| `yelmo_data.pd_vel_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_Velocity-Rignot2011.nc` | fesmdata | yes | none | pd_vel_load = True |
| `isos.rheology_file` | `isostasy_data/earth_structure/yelmo/ANT-8KM_GIA_HR24.nc` | isostasy_data | yes | none | lithosphere/mantle = rheology_file |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `ghf.obs_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_GHF-HazzardRichards2024.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `gcm_ts_ref.filename` | `ice_data/Antarctica/ANT-8KM/RACMO2.3/ANT-8KM_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc` | ice_data | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_pr_ref.filename` | `ice_data/Antarctica/ANT-8KM/RACMO2.3/ANT-8KM_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc` | ice_data | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_zs_ref.filename` | `ice_data/Antarctica/ANT-8KM/RACMO2.3/ANT-8KM_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc` | ice_data | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_to_ref.filename` | `ice_data/Antarctica/ANT-8KM/ERA-INT-ORAS5/ANT-8KM_ORAS5_1979-2024_thetao.nc` | ice_data | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_so_ref.filename` | `ice_data/Antarctica/ANT-8KM/ERA-INT-ORAS5/ANT-8KM_ORAS5_1979-2024_so.nc` | ice_data | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `sed.obs_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_Sediment-CRUST1.nc` | yes | sed.method = 0 |
| `marine_shelf.obs_path` | `fesmdata/Antarctica/ANT-8KM/ANT-8KM_BasalMelt-Rignot2013.nc` | yes | use_obs = False |
| `marine_shelf.tf_path` | `ice_data/ISMIP6/Antarctica/ANT-8KM/Ocean/ANT-8KM_OCEAN_ISMIP6_J20.nc` | no | tf_correction = False |

### eGRL (esm_Greenland)

`yelmox_esm_Greenland.nml`, domain `Greenland`; grids: hub GRL-8KM, ice GRL-8KM, isos GRL-8KM, clim GRL-8KM, surface GRL-8KM, shelf GRL-8KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `ice_data/ISMIP7/Greenland/GRL-8KM/obs/GRL-8KM_ObsISMIP7-v1.3.nc` | ice_data | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_ismip7` | `ice_data/ISMIP7/Greenland/GRL-8KM/obs/GRL-8KM_ObsISMIP7-v1.3.nc` | ice_data | yes | none | regions layer ismip7 (var mouginot_basins) |
| `yelmo_data.pd_tsrf_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | fesmdata | yes | none | pd_tsrf_load = True |
| `yelmo_data.pd_smb_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | fesmdata | yes | none | pd_smb_load = True |
| `yelmo_data.pd_vel_path` | `ice_data/ISMIP7/Greenland/GRL-8KM/obs/GRL-8KM_ObsISMIP7-v1.3.nc` | ice_data | yes | none | pd_vel_load = True |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `sed.obs_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_Sediment-CRUST1.nc` | fesmdata | yes | none | sed.method = 1 |
| `ghf.obs_path` | `ice_data/ISMIP7/Greenland/GRL-8KM/obs/GRL-8KM_ObsISMIP7-v1.3.nc` | ice_data | yes | none | ghf.use_obs = True |
| `gcm_ts_ref.filename` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_SMB-MARv3.14-1961-1990.nc` | fesmdata | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_smb_ref.filename` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_SMB-MARv3.14-1961-1990.nc` | fesmdata | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_zs_ref.filename` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_SMB-MARv3.14-1961-1990.nc` | fesmdata | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False) |
| `gcm_to_ref.filename` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False); grid_src=GRL-16KM remap=bilin |
| `gcm_so_ref.filename` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | esm spinup reference (run_type = spinup; use_esm = use_var = False); grid_src=GRL-16KM remap=bilin |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the present-day topography |
| `yelmo_data.pd_age_path` | `fesmdata/Greenland/GRL-8KM/GRL-8KM_Stratigraphy-MacGregor2025.nc` | yes | pd_age_load = False |
| `marine_shelf.obs_path` | `ice_data/Greenland/GRL-8KM/GRL-8KM_BMELT-R13.nc` | no | use_obs = False |
| `marine_shelf.tf_path` | `ice_data/ISMIP6/Greenland/GRL-8KM/Ocean/GRL-8KM_OCEAN_ISMIP6_J20.nc` | no | tf_correction = False |
| `barysealevel.sl_path` | `input/paleo300ka_hybrid_aicc2012.nc` | yes (stack `input/`) | method = fastiso |

### rGRL (rembo_Greenland)

`yelmox_rembo_Greenland.nml (+ rembo1 `rembo_Greenland.nml`)`, domain `Greenland`; grids: hub GRL-16KM, ice GRL-16KM, isos GRL-16KM, clim GRL-16KM, surface GRL-16KM, shelf GRL-16KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_TOPO-BedMachine-v6.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_basins` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_BASINS-Zwally2012.nc` | fesmdata | yes | none | basin set Zwally2012 |
| `yelmo_data.pd_tsrf_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | fesmdata | yes | none | pd_tsrf_load = True |
| `yelmo_data.pd_smb_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | fesmdata | yes | none | pd_smb_load = True |
| `yelmo_data.pd_vel_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Velocity-Joughin2018.nc` | fesmdata | yes | none | pd_vel_load = True |
| `yelmo_data.pd_age_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Stratigraphy-MacGregor2025.nc` | fesmdata | yes | none | pd_age_load = True |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `sed.obs_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Sediment-CRUST1.nc` | fesmdata | yes | none | sed.method = 1 |
| `ghf.obs_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_GHF-Shapiro2004.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `tm_1D.file` | `input/timeout_ramp_100kyr.txt` | input | yes (stack `input/`) | none (gate overlay supplies it) | tm_1D.method = file |
| `oras5_1981_2010_to.filename` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) to |
| `oras5_1981_2010_so.filename` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_Ocean-ORAS5-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) so |
| `topo_file` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_TOPO-BedMachine-v6.nc` | fesmdata | yes | none | rembo1 topo (z_srf,H_ice,z_bed) |
| `clim_file` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc` | fesmdata | yes | none | rembo1 precip (rf, sf, annual mean) |
| `precip_mon_file` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1961-1990.nc` | fesmdata | yes | none | rembo1 monthly precip reference (ppcorr) |
| `emb_clim_file` | `ice_data/Greenland/GRL-80KM/ERA-40/GRL-80KM_ERA-40_195801-200112.nc` | ice_data | yes | none | rembo1 EMB boundary climate (ERA-40, 80 km), always read |
| `forcing_file` | `in_rembo/forcing/mis11.climber_paleo.dat` | in_rembo | n/a | none | rembo1 forcing_file: read as a parameter only (cparam); not opened in rembo1 v1.3 sources (checked: no nc_read/open of it) |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Greenland/GRL-16KM/GRL-16KM_TOPO-BedMachine-v6.nc` | yes | the hub supplies the present-day topography |
| `marine_shelf.obs_path` | `ice_data/Greenland/GRL-16KM/GRL-16KM_BMELT-R13.nc` | no | use_obs = False |
| `barysealevel.sl_path` | `input/sl_deglaciation.dat` | yes (kryos `input/`) | method = fastiso |

### ANT-full (shipped Antarctica, relax = opt = True; extras to phase 1 only)

`yelmox_Antarctica.nml`, domain `Antarctica`; grids: hub ANT-16KM, ice ANT-32KM, isos ANT-32KM, clim ANT-32KM, surface ANT-32KM, shelf ANT-16KM.

| key | file | expected | found | action | read when |
|---|---|---|---|---|---|
| `topo.path` | `fesmdata/Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine-v4.nc` | fesmdata | yes | none | hub topography (htopo_init) |
| `regions.path_regions` | `fesmdata/Antarctica/ANT-16KM/ANT-16KM_REGIONS.nc` | fesmdata | yes | none | regions/zones (read at hub, remapped nn to ice/clim/shelf grids) |
| `regions.path_basins` | `fesmdata/Antarctica/ANT-16KM/ANT-16KM_BASINS-Zwally2012.nc` | fesmdata | yes | none | basin set Zwally2012 |
| `yelmo_data.pd_tsrf_path` | `ice_data/Antarctica/ANT-32KM/RACMO2.3/ANT-32KM_RACMO23-ERAINT-HYBRID_1981-2010.nc` | ice_data | yes | none | pd_tsrf_load = True |
| `yelmo_data.pd_smb_path` | `ice_data/Antarctica/ANT-32KM/RACMO2.3/ANT-32KM_RACMO23-ERAINT-HYBRID_1981-2010.nc` | ice_data | yes | none | pd_smb_load = True |
| `yelmo_data.pd_vel_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_Velocity-Rignot2011.nc` | fesmdata | yes | none | pd_vel_load = True |
| `isos.rheology_file` | `isostasy_data/earth_structure/yelmo/ANT-32KM_GIA_HR24.nc` | isostasy_data | yes | none | lithosphere/mantle = rheology_file |
| `barysealevel.A_ocean_path` | `isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc` | isostasy_data | yes | none | method = fastiso |
| `ghf.obs_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_GHF-Shapiro2004.nc` | fesmdata | yes | none | ghf.use_obs = True |
| `racmo23p2_1981_2010_tas.filename` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) tas |
| `racmo23p2_1981_2010_pr.filename` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) pr |
| `racmo23p2_1981_2010_zs.filename` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010.nc` | fesmdata | yes | none | snapshot_pd (reference) zs |
| `ismip7_to.filename` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_Ocean-ISMIP7.nc` | fesmdata | yes | none | snapshot_pd (reference) to |
| `ismip7_so.filename` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_Ocean-ISMIP7.nc` | fesmdata | yes | none | snapshot_pd (reference) so |

Keys present in the resolved configuration but not read (existence checked anyway):

| key | file | found | why not read |
|---|---|---|---|
| `yelmo.grid_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_grid.nc` | yes | the driver passes the grid (maps/grid_<name>.txt) |
| `yelmo_init_topo.init_topo_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the initial topography |
| `yelmo_data.pd_topo_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine-v6.nc` | no | the hub supplies the present-day topography |
| `sed.obs_path` | `fesmdata/Antarctica/ANT-32KM/ANT-32KM_Sediment-CRUST1.nc` | yes | sed.method = 0 |
| `marine_shelf.obs_path` | `fesmdata/Antarctica/ANT-16KM/ANT-16KM_BasalMelt-Rignot2013.nc` | yes | use_obs = False |
| `marine_shelf.tf_path` | `ice_data/Antarctica/ANT-16KM/ANT-16KM_OCEAN_ISMIP6_J20.nc` | yes | tf_correction = False |
## Keys that are present but not read, and are partly absent

Not read in yelmox R2', so not read by kryos either; absent files of this kind, collected from the tables above:

- `yelmo_init_topo.init_topo_path` and `yelmo_data.pd_topo_path` default to `..._TOPO-BedMachine-v6.nc`; the file is absent for ANT-8KM, ANT-32KM, LIS-32KM, NH-32KM, PYR-250M and SRG-250M (present for GRL grids). Never opened: the hub supplies both (already so in v2.4: "the initial topography comes from the hub"), which is why ANT in phase 1 runs without `ANT-32KM_TOPO-BedMachine-v6.nc`. The kryos `--check-par` and the yelmo adapter must not stat them.
- `marine_shelf.obs_path` (`BMELT-R13`): absent for GRL-16KM, GRL-8KM, GRL-PAL-16KM, LIS-32KM, NH-32KM, SRG-250M; `use_obs = False` in every base. Same for the ISMIP6 `tf_path` of ANT-8KM and GRL-8KM (`tf_correction = False`).
- SRG `sed.obs_path` (`SED-L97`) and `ghf.obs_path` (`GHF-S04`): absent; `sed.method = 0`, `ghf.use_obs = False` (constant `ghf_const`).
- GRP `pd_smb_path` and `pd_vel_path` (GRL-PAL-16KM MAR annual mean, Joughin2018): absent; `pd_smb_load = pd_vel_load = False`.
- eGRL/eANT transient set (`gcm_*_esm_ref/hist/proj`): not in the spinup gate. The ISMIP7 Greenland `historical` files exist (`ice_data/ISMIP7/Greenland/GRL-8KM/CESM2-WACCM/historical/...`); the `ctrl` projection files do not exist, and are not read for `experiment = "ctrl"`. For eANT the groups of `esm_ant.nml` contain `{scenario}`, which `esm_forcing.f90` never substitutes, `esm.esm_name = "UKESM"` has no folder in `ice_data/ISMIP6-2300/Antarctica/ANT-8KM/{Atmosphere,Ocean}` (there are `CCSM4_RCP85`, `CESM2-WACCM_ssp585`, ...), and `gcm_so_hist` points to `ice_data/cmip-2300/...`, which does not exist. Irrelevant for the gate; relevant for the transient eANT run and WP9 item 4.

## Maps

Present in `$G/maps/S` (frozen, 5 maps): `scrip-bilin_ANT-16KM_ANT-32KM`, `scrip-bilin_ANT-32KM_ANT-16KM`, `scrip-con_ANT-16KM_ANT-32KM`, `scrip-nn_ANT-16KM_ANT-32KM`, `scrip-bilin_GLOBAL-1DEG_GRL-16KM`. Grid descriptions present for all grids of the 8 bases (`grid_ANT-8KM`, `GRL-8KM`, `GRL-PAL-16KM`, `LIS-32KM`, `NH-32KM`, `PYR-250M`, `PYR-2KM`, `SRG-16KM`, `SRG-250M`); no `grid_GLOBAL-1DEG.txt` (as for the existing GLOBAL map).

Missing, per base (derivation: hub <-> ice maps if grids differ; regions nearest neighbour hub -> clim/shelf/ice; isos load/displacement ice <-> isos with the smooth rule (bilin refine, con coarsen); hub -> clim/surface/shelf; clim -> surface bilin; file grid -> target grid in varslice/sediments/geothermal):

| base | pair (src -> dst) | method | needed for | seen in mapgen (23:59) |
|---|---|---|---|---|
| GRP, eANT, rGRL | none (all components one grid; esm/rembo/era5 files are on that grid) | | | no new map |
| LIS | GLOBAL-1DEG -> LIS-32KM | bilin | PMIP3 snapshots (`grid_src`) | yes |
| NOR | GLOBAL-1DEG -> NH-32KM | bilin | PMIP3 snapshots | yes |
| PYR | PYR-250M -> PYR-2KM | bilin, con, nn | hub -> clim/isos (lapse transform, regions, isos load) | yes (all three) |
| PYR | PYR-2KM -> PYR-250M | bilin (con for kryos' clim -> ice if used) | isos displacement, clim -> surface, `ghf` file on PYR-2KM | bilin yes |
| SRG | SRG-250M -> SRG-16KM | con | isos load at init | yes |
| SRG | SRG-16KM -> SRG-250M | bilin (superset) | isos displacement; yelmox does not need it with `with_isostasy = False` (not generated), kryos with `geo = data const` should not | no |
| eGRL | GRL-16KM -> GRL-8KM | bilin | esm `to_ref`, `so_ref` (`grid_src = GRL-16KM`) | yes |

The kryos wiring may ask for more pairs than yelmox (it fails on the read-only folder, as intended); `con`/`nn` in the reverse direction of the PYR pair are the likely candidates.

## Links made

None. Commands run on Levante were read-only (`ls`, `stat`, `ncdump -h`, `cat`); nothing was written under `$G` or `/work/ba1442/{fesmdata,ice_data,isostasy_data}`.

## For Alex

1. **kryos `input/` is not self-contained for three bases.** `alpha_orbital_t2018.dat` (PYR), `sealevel_spratt2016.dat` (GRP) and `timeout_ramp_100kyr.txt` (rGRL) are only in yelmox `input/`. The gate's `input_overlay` and `run_kryos.sh` (which links the stack's `input/*` first) supply them, so the gates run; copy them into kryos `input/` (committed) if the overlay is to go. Not a fesmdata item, so not linked. `LA2004.*` and `sl_deglaciation.dat` are in kryos `input/`. (The clone `$G/kryos` on Levante is at the `af8f1b9` era; its `input/` lacks the kryos-owned series, which the overlay also covers.)
2. **rGRL `tm_1D.method = "file"`** (the other bases use `const`): the phase-1 `-p tm_1D.dt=10` is inert and the 1D output follows `input/timeout_ramp_100kyr.txt` (0:10:200, 200:20:300, ...). Add `tm_1D.method=const` to the rGRL `-p`, or accept it; the same on both sides, so it only changes the output times.
3. **rGRL, rembo1:** rembo1 reads its parameters from `<run folder>/rembo_Greenland.nml` (yelmox `yelmox/rembo_Greenland.nml`, staged by runme) and `./input/LA2004.INSOLN.5Ma.txt`, `LA2004.INSOLP.1Ma.txt`; its `forcing_file` (`in_rembo/forcing/...`) is only parsed, never opened; it opens `topo_file`, `clim_file`, `precip_mon_file` (GRL-16KM, MAR3.14, present) and `emb_clim_file` (ERA-40 80 km, present). The kryos run folder needs that par file staged; rembo1 is not linked in S3 (`rembo=0`): stack S4.
4. **eGRL hub topography is `ice_data/ISMIP7/Greenland/GRL-8KM/obs/GRL-8KM_ObsISMIP7-v1.3.nc`** (`topo.names = bed thickness surface_grimp`), not a fesmdata product, and has no `grid_name` attribute (size 211 x 361 = GRL-8KM). Same file for `pd_vel`, `ghf` and the `mouginot_basins` layer. Plan expected grid/topo in fesmdata.
5. **eGRL ocean reference** `fesmdata/Greenland/GRL-16KM/GRL-16KM_Ocean-ORAS5-1981-2010.nc` is remapped onto GRL-8KM (the new bilin map above).
6. **esm / rembo / smb_simple have no kryos templates yet** (`input/*_defaults.nml`): WP work, not data.
7. The mapgen agent's folders (`$G/maps/gen/*-maps`) show what yelmox needs (table above); kryos may need more.
