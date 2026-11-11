---
title: "VILMA2"
---

The geo backend `vilma2` couples the global viscoelastic solid-Earth and sea-level
model VILMA2 (library `libvilma`) to kryos. It has two parts: one **global** instance
that owns the Earth on its own Gauss grid and solves the sea-level equation on the
load of the whole globe, and one stateless **regional port** per domain that receives
the bed and sea-surface change on the domain's grid. The ice of every domain enters
the global load at every Earth step, so regional ice sheets and the global sea level
see each other (glacial isostatic adjustment, GIA, with rotation and a migrating
coastline).

- **Build:** `vilma=1`. `config/common.mk` expects `vilma2/libvilma/include`
  (`libvilma.a`), and the SHTns that fesm-utils builds (`make shtns`); VILMA2 must be
  built with the same `openmp=` as kryos and fesm-utils. Without it the stub stops the
  run.
- **Code:** [fesmc/vilma2](https://github.com/fesmc/vilma2), documentation
  <https://fesmc.github.io/vilma2/>.
- **Select:** `[globe] geo = "<group>"` with `backend = "vilma2"` (the global part),
  and `geo = "<group>"` with `backend = "vilma2"` in each domain (the port). The keys
  are the template `[geo_vilma2]` of
  [`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml);
  the physics are the group `&vilma`, overrides of
  [`input/vilma_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/vilma_defaults.nml)
  (a synchronised copy of VILMA2's own file).
- **Source:** [`src/backends/geo/kryos_geo_vilma2.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/geo/kryos_geo_vilma2.f90);
  the gather and the reference in `src/core/kryos_run_gia.f90` and
  `src/core/kryos_globe_ice.f90`.
- **Par files:** `par/kryos_NH-ANT_plastic_vilma2_b1.nml` (Gowan's 1-D Earth),
  `par/kryos_NH-ANT_plastic_vilma2.nml` (the 3-D Earth), the Yelmo science case
  `par/kryos_NH-ANT_vilma2.nml`; see [Par files](reference-par-files.md).

## Global instance and regional ports

The global instance (`geo_glob` in the par files) is global only, one per run. It
builds the Gauss grid `GAUSS-L<lmax>` that VILMA2 itself uses: `nlat = 2 lmax + 2`,
`nphi = 4 lmax` unless `&vilma` sets `nlat`, `nphi` (> 0), longitudes from 0, a sphere
of VILMA2's radius `r_earth`, so that integrals over the grid are VILMA2's own
quadrature. Its grid description is written to `<map_fldr>/grid_<name>.txt` for the
couplers (an existing file must agree). The par files leave `grid = ""`. VILMA2 itself
defaults `nlat = 0` to 2 lmax + 2 and `nphi = 0` to 4 lmax (`vilma_drive.f90`,
`vilma_coupling.f90`), as kryos does.

The core turns the domains' ice into VILMA2's load (see [Global hub](global-hub.qmd)):
each domain's ice (lag 1: the ice of the previous step) becomes the host terms
`host_w` (coverage of the host regions), `host_Ig` (grounded ice), `host_C` (ocean
function) and `host_F` (flotation term) on the Gauss grid, and the ice outside the
domains is the background of `&globe_ice` (key `ice_bg`) plus, with
`ext = "file"`, the anomaly of an ice history (`H_ice_ext`). VILMA2 ramps from the
load of the previous call to this one over the step and returns the bed and
sea-surface change; the core scatters `dz_bed`, `dz_ss`, `dphi` and `bsl` into the
ports.

A regional port (`geo_nor`, `geo_ant`, ...) lies on the grid the chain gives it (the
ice grid by default), holds only the scattered anomalies and passes the ice's static
reference bed through. Its cadence must be that of the global instance (a different
`dt` stops the run). It writes no restart file.

The densities are the run's (`phys_const`: `rho_ice`, `rho_sw`), not keys: the backend
replaces VILMA2's `rho_ice` and `rho_water` and logs both. Gravity comes from VILMA2's
Earth model. `solver` must be `"v2"` and `sle_fixed_ocean` false (the host terms need
the native solver and a migrating coastline).

## Keys

Template `[geo_vilma2]` (global instance), besides the core keys of `[instance]`:

| Key | Default | |
|---|---|---|
| `vilma` | `"vilma"` | the `&vilma` group; must be named `vilma` (VILMA2 reads its defaults under that name only) |
| `ice_bg` | `"globe_ice"` | the global ice field group |
| `allow_no_spinup` | `False` | a run without VILMA2's spin-up (neither `pre_spinup_1d` nor `equil_time_max > 0`) stops unless True (tests, present-day starts) |
| `out_gia` | `"none"` | schedule of `<instance>_gia.nc`: the global ice field VILMA2 got |
| `out_vilma` | `"none"` | schedule of `<instance>_vilma.nc`: VILMA2's own diagnostics |

The `&vilma` keys used by the runs (defaults from `input/vilma_defaults.nml`):

| Key | Default | |
|---|---|---|
| `lmax` | 64 | maximum spherical-harmonic degree (the runs: 170) |
| `nlat`, `nphi` | 0, 0 | Gauss grid size (0: from lmax) |
| `earth` | `"M3-L70-V01"` | a built-in model, `"PREM"` or `"custom"`; `"PREM"` takes `n_layer`, `r_bot`, `r_top`, `eta`, `rheology` and fills `rho`, `mu` from the r²-weighted shell means of incompressible PREM |
| `n_layer`, `r_top`, `r_bot` | 5, ... | layers, surface first [m] |
| `eta` | ... | [Pa s] viscosity per layer |
| `rheology` | ... | per layer: 0 elastic, 1 Maxwell, 2 fluid |
| `rotation` | `.true.` | rotational feedback (true polar wander) |
| `l_visc_3d` | `.false.` | load a lateral log10(η) field |
| `visc_3d_file` | `"fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"` | the 3-D field |
| `visc_3d_structure` | 6 | ≥ 1: the file is in the FesmData layout and this is the structure (6 = `v_1.0_s16`) |
| `visc_log10_min`, `visc_log10_max` | 19.5, 30.0 | [dex] floor and ceiling of log10(η) |
| `pre_spinup_1d` | `.false.` | a cheap 1-D pre-equilibration before the full model |
| `equil_time_max` | 0.0 | [yr] > 0: cap of the full-model spin-up |
| `equil_rate_tol` | 1.0e-3 | [m yr-1] mean bed speed over a pass that ends the spin-up |
| `scheme` | `"fe"` | viscoelastic memory scheme |
| `deg1_frame` | `"cm"` | degree-1 reference frame (centre of mass) |

## Earth structure: 1-D and 3-D

Two Earths are in use; both are `earth = "PREM"` at `lmax = 170` with rotation.

- **Gowan's 1-D Earth** (`_b1` par files, as the standalone run B1 of
  `docs/dev/vilma2-paleomist.md`): 12 layers, elastic 0-80-120 km, Maxwell
  4×10²⁰ Pa s to 670 km (split at 220, 400, 600 km), 4×10²² Pa s to the core-mantle
  boundary, fluid core. The radial profile of Gowan et al. (2021) (120 km lithosphere,
  upper mantle 4×10²⁰, lower mantle 4×10²² Pa s), on which PaleoMIST was built.
- **The 3-D science Earth** (as `par/kryos_NH-ANT_vilma2.nml`, the CLIMBER-X set):
  11 layers, 80 km elastic lithosphere, the 1-D radial reference (upper mantle
  4×10²⁰, lower mantle 10²² Pa s) and the lateral viscosity of Bagge et al. (2021)
  from FesmData (`GLOBAL-T170_SolidEarth-Bagge2021.nc`, structure 6 = `v_1.0_s16`,
  floored at 10^19.5 Pa s). With `l_visc_3d = .true.` the file's log10 η (bilinear
  onto the Gauss grid, linear in radius onto the FE nodes, clamped to
  `visc_log10_min`, `visc_log10_max`) replaces the radial viscosity of every Maxwell
  layer; the elastic lithosphere (80 km) and the fluid core are kept. A layer whose
  lateral spread is below `visc3d_tol` runs 1-D on its lateral mean.

The `&vilma` groups of the two par files, `par/kryos_NH-ANT_plastic_vilma2_b1.nml`
(1-D) and `par/kryos_NH-ANT_plastic_vilma2.nml` (3-D):

```fortran
&vilma                          ! 1-D: Gowan et al. (2021), as B1
    lmax              = 170
    rotation          = .true.
    earth             = "PREM"
    n_layer           = 12
    r_top             = 6371.0e3, 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
    r_bot             = 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
    eta               = 1.0e40, 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 0.0
    rheology          = 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
    l_visc_3d         = .false.
    pre_spinup_1d     = .false.
    equil_time_max    = 2.0e6     ! [yr] relaxation to the 80 ka load (B1)
    equil_rate_tol    = 1.0e-5    ! [m/yr]
/
```

```fortran
&vilma                          ! 3-D: PREM 11 layers + Bagge et al. (2021)
    lmax              = 170
    earth             = "PREM"
    n_layer           = 11
    r_bot             = 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0
    r_top             = 6371.0e3, 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3
    eta               = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 0.0
    rheology          = 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
    l_visc_3d         = .true.
    visc_3d_file      = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"
    visc_3d_structure = 6
    pre_spinup_1d     = .true.    ! the 1-D pre-equilibration
    equil_time_max    = 2.0e4     ! [yr] cap of the full-model relaxation after it
    equil_rate_tol    = 1.0e-5    ! [m/yr] mean bed rate that ends a pass (both phases)
/
```

The 3-D file leaves `rotation` out (default `.true.`). Besides `&vilma`, the two files
differ only in the reference (below): the ice, the margins, τ_b and the clock are the
same.

| Group | 1-D (`_b1`) | 3-D |
|---|---|---|
| `[X_z_bed]`, `[X_H_ice]` (background X of `&globe_ice`) | `fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc` (`z_bed`, `H_ice`), `mode = "transient"`, read at `ext_tref` = 0: PaleoMIST's 0 ka bed and ice (B1's reference) | `kryos_data/refcorr/PaleoMIST-sci/{grid_name}/{grid_name}_REFCORR.nc`, `mode = "static"`: PaleoMIST's 0 ka bed + r, its 0 ka ice |
| `[topo_nor]` `path` | `{grid_name}_TOPO-GEBCO2026.nc`: today's topography | `{grid_name}_TOPO-GEBCO2026-refcorr.nc`: today's + r |
| `[topo_ant]` `path` | `{grid_name}_TOPO-BedMachine-v4.nc` | `{grid_name}_TOPO-BedMachine-v4-refcorr.nc` |

r is the present-day reference correction: the rebound still to come at 0 ka after the
glacial cycle, computed once per Earth and ice history by iterating standalone VILMA2
under the history (`tools/levante/gia_refcorr.sh`), so that the coupled run ends on
today's bed instead of one 100-300 m too deep under the former ice sheets. The 1-D
runs take the present-day bed without correction as the relaxed reference, as B1.
In both, the reference ice is present-day: inside the domains the ice backend's
`H_ice_pd` (see [PlasticIceSheet](backend-plastic.md#exports)), outside them X.

In the standalone comparison under PaleoMIST (`docs/dev/vilma2-paleomist.md`), the
3-D Earth changes the deglacial RSL by about 50 m RMS over the formerly glaciated
areas (up to 230 m), mostly through its radial profile; lmax 170 differs from 256 by
0.4 m RMS. A 3-D run costs far more: a standalone run at lmax 170 (spin-up and the
80 kyr transient) took 1 h 44 on 64 threads in 3-D against 3 min in 1-D.

## Spin-up

A cold start holds the load of the first time (here the 80 ka ice) and relaxes the
Earth under it, so that the transient starts with viscous memory, then seeds the
state with a dt = 0 solve under the same load. The reference stays the datum.

- `equil_time_max > 0`: the full-model relaxation, in passes (in the standalone runs:
  200 yr, growing by a factor 1.6), until the mean bed speed over a pass is below
  `equil_rate_tol`, or the cap is reached (warning).
- `pre_spinup_1d = .true.`: first a cheap 1-D pre-equilibration on the radial
  geometric mean of the 3-D viscosity (the unweighted mean of log10 η over the Gauss
  grid), then the full model carrying the spun-up memory.

The 1-D runs relax up to 2×10⁶ yr at 10⁻⁵ m yr-1 (B1: converged after 93 kyr of
relaxation); the 3-D runs pre-spin in 1-D and cap the 3-D relaxation at 2×10⁴ yr
(about 14 kyr to converge at 10⁻⁵ m yr-1, per the par file). The 1-D pre-spin-up has
no time cap: it relaxes to the same `equil_rate_tol`, in passes of 200 yr growing by
1.6, up to an internal limit of 60 passes (no warning if it is reached);
`equil_time_max` caps only the full-model relaxation that follows.

## Rotation

`rotation = .true.` (the default) adds the rotational feedback (true polar wander) to
the sea-level equation: `rotation_k_s` (≤ 0: the model's own fluid limit) and
`rotation_c_minus_a` (C − A [kg m²]). Gowan et al. (2021) included rotation as well;
`.false.` is for the non-rotating benchmarks.

## Output

- **`out_2D`, `out_1D`:** the exports on the Gauss grid (lon, lat) and their series
  (`geo_glob.nc`, `geo_glob_ts.nc`), table `geo`: `dz_bed` (bed displacement),
  `dz_ss` (sea-surface change relative to the uniform shift `dphi`), `C_ocn`, `bsl`,
  `bsl_ive` (ice-volume equivalent), `bsl_ive_host`, `A_ocn`, `mass_resid`.
- **`out_gia`:** `<instance>_gia.nc`, the global ice field VILMA2 got (`H_ice_ext`,
  `host_w`, `host_Ig`, `host_C`, `host_F`; table
  [`input/kryos-variables-geo-vilma2.md`](https://github.com/fesmc/kryos/blob/main/input/kryos-variables-geo-vilma2.md)).
- **`out_vilma`:** `<instance>_vilma.nc` (`geo_glob_vilma.nc`), VILMA2's own
  diagnostics file, rounded by its table
  [`input/vilma-variables.md`](https://github.com/fesmc/kryos/blob/main/input/vilma-variables.md):
  `h_ice` (the grounded-ice load of the solve: with the host terms of ice instances
  I_g + (1 − w) h_ice_ext, the instances' ice inside their weight w and the external ice
  outside), `h_ice_ext` (the external ice of the 1 − w part, e.g. PaleoMIST outside the
  plastic grids; = `h_ice` without host terms; runs of before 2026-11-11 have only it, as
  `h_ice`), `rsl` (relative sea-level change, full field), `z_bed` (= z_bed_eq − rsl),
  `C_ocean`, `dz_bed`, `dz_ss` (rsl = dz_ss − dz_bed), and the series `bsl`, `dphi`,
  `bsl_ive`, `bsl_ive_host`, `fill_int`, `fill_int_own`. Not the viscoelastic memory.
- **Restart:** VILMA2's restart file, labelled with its own clock, always lossless.

The plastic par files write `geo_glob_vilma.nc` every 2.5 kyr (`out_vilma = "out_2D"`,
the PaleoMIST slice times).

### GAPSLIP scoring

`geo_glob_vilma.nc` is what the sea-level scoring reads. `tools/levante/pm_compare.jl`
(`read_run`, in place) takes `rsl`, `dz_bed`, `bsl` and `bsl_ive`; the model RSL at an
indicator is rsl(t) − rsl(0 ka) (relative to present, as Gowan's), bilinear in space
and linear in time between the records, against GAPSLIP 2.1 (Gowan 2023; FesmData
`Records/RSL/RSL-GAPSLIP_global.nc`). Index points give z = (model − obs)/σ, limiting
data are violated or not, by group (Laurentide, Greenland, Fennoscandia/Barents,
Antarctica, far field). The same reading drives `tools/levante/plastic_score.sh`, the
comparison across Earths (`plastic_earth_compare.jl`) and the τ_b calibration
(`plastic_esmda.jl`), whose members write only `geo_glob_vilma.nc` and the 1-D series.
See [Calibration and ensembles](experiment-calibration.md) and
[Analysis and figures](experiment-analysis.md).

## Configuration

The Earth part of `par/kryos_NH-ANT_plastic_vilma2_b1.nml`:

```fortran
&globe
    geo = "geo_glob"
/

&geo_glob
    backend         = "vilma2"
    grid            = ""                 ! the Gauss grid GAUSS-L<lmax>, built by the backend
    dt              = 100.0              ! [yr] cadence of the Earth (science case: 10)
    out_2D          = "out_2D"
    out_1D          = "out_1D"
    vilma           = "vilma"            ! the &vilma group
    ice_bg          = "globe_ice"        ! the global ice field group
    allow_no_spinup = False
    out_gia         = "none"
    out_vilma       = "out_2D"           ! geo_glob_vilma.nc: rsl, bsl, bsl_ive (the GAPSLIP scores)
/

&globe_ice
    grid      = "GLOBAL-0.5DEG"          ! background grid
    z_bed_ref = "X_z_bed"                ! background bed X at ext_tref
    H_ice_ref = "X_H_ice"                ! background ice X at ext_tref
    ext       = "file"                   ! ice outside the footprints: X + the history's anomaly
    ext_field = "X_dH_ice"
    ext_tref  = 0.0
/

&geo_nor
    backend     = "vilma2"
    grid        = "NH-32KM"              ! the ice grid
    dt          = 100.0                  ! the cadence of geo_glob (a different value stops)
/
```

Outside the plastic domains the ice is PaleoMIST's history (`X_dH_ice`: `H_ice` of
`IceHistory-PaleoMIST` as the anomaly against its 0 ka slice).
