---
title: "PaleoEmulator (climate)"
---

The climate backend `paleoemu` couples the CESM-3Ma climate emulator of
[PaleoEmulator](https://github.com/fesmc/PaleoEmulator) to kryos. The emulator
returns monthly 2-m temperature and precipitation anomalies against pre-industrial,
and the annual ocean temperature change at 8 levels between 4 and 1475 m.
Its inputs are the orbit, CO2 and the Northern Hemisphere ice sheets: the regional
grounded ice volumes, their slow filters and the surface-elevation change. It was
trained on the 3-Myr transient CESM simulation of Yun et al. (2023), and its
uncertainty is calibrated.

Like VILMA2, it has two parts:

- **The global part.** One instance in the globe steps the emulator on the ice of
  the whole Northern Hemisphere: the domains' ice merged into a background ice
  history (`globe_ice`).
- **The regional ports.** Each domain's climate is a port that receives the result
  on the domain's grid.

The surface (chion) takes `tas` and `pr` from the port through the lapse rows, and
computes its insolation on the same orbit. The surface's `smb` and `T_srf` then
drive the ice. A port with a present-day ocean also gives the shelf its ocean
profiles (set A), the present day plus the emulator's change ([The ocean](#the-ocean)).

- **Build:** `paleoemu=1`.
  - `config/common.mk` expects `PaleoEmulator/fortran/build` (`libpaleoemu.a`,
    `paleoemu.mod`), a sibling checkout or symlink `PaleoEmulator`.
  - Build it with the compiler and the fesm-utils of kryos: `make paleoemu-static`
    (it is serial and uses ncio).
  - Without it, the stub stops the run.
- **Code:** [fesmc/PaleoEmulator](https://github.com/fesmc/PaleoEmulator).
  - The Fortran module is `fortran/paleoemu.f90`.
  - The API, the emulator file and its formulas are described in the README.
- **Select:**
  - The global part: `[globe] climate = "<group>"` with `backend = "paleoemu"`.
  - The port: `climate = "<group>"` with `backend = "paleoemu"` in each domain.
  - The keys are in the template `[climate_paleoemu]` of
    [`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml).
- **Source:**
  - [`src/backends/climate/kryos_climate_paleoemu.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/climate/kryos_climate_paleoemu.f90).
  - The ice geometry is in `globe_ice_geometry` and `globe_ice_external`
    (`src/core/kryos_globe_ice.f90`).
  - The globe side is in `globe_climate_ice`, `globe_climate_reference` and
    `globe_forcing_to_climate` (`src/core/kryos_globe.f90`).
- **Par files:** `par/kryos_North_plastic_paleoemu.nml`: the PaleoMIST plastic ice of
  North with chion, from 21 to 20 ka; `par/kryos_NH-ANT_data_paleoemu.nml` (North and
  Antarctica, 80 ka to 0) and `par/kryos_GRL-ANT_data_paleoemu_pd.nml` (Greenland and
  Antarctica at present day): PaleoMIST's or today's ice imposed (ice backend `data`)
  with chion, and the shelf diagnostics ([Shelf diagnostics](#shelf-diagnostics)).
- **Inputs:** `kryos_data/paleoemu/` (its README): the emulator file, the reference
  climate on GLOBAL-0.5DEG and the CO2 table, in `$KRYOS_DATA` on Levante and albedo.

## The global part

The global instance (`clim_glob` in the par file) is global only, one per run.

**Grid.** It runs on its own lon-lat grid, the key `grid`. This is the background grid
of `[globe_ice]` (`GLOBAL-0.5DEG`), a `grid_<name>.txt` of the map folder. The grid
must be at least as fine as the emulator's ice grid (1.5° × 0.75°): each ice-grid cell
takes the area-weighted mean of the instance cells whose centres it holds.

**Imports.** The core fills all of them except the coefficient:

| Import | From | When |
|---|---|---|
| `H_ice`, `z_bed`, `z_srf` | `globe_ice_geometry`: the ice of the last gather of the domains (lag 1), merged into the background history (owner weights; the ice extensive, the bed and the surface intensive) | every update |
| `H_ice_ref`, `z_bed_ref`, `z_srf_ref` | the spliced present day of `globe_ice`: the domains' reference in the background X | once, after the reference |
| `coef_co2` | the forcing (`[forcing] coefs = "co2"`), through the globe's own rows `forcing -> climate_g` (lag 0) | every update |
| `orbit_ecc`, `orbit_obl`, `orbit_per` (optional) | the forcing (`[forcing] orbit = "orbit"`), the same rows | every update |

The surface where `globe_ice` has no surface of its own is `z_bed + H` where grounded,
and the floating surface `z_sl + (1 − ρ_i/ρ_sw) H` elsewhere.

**Update.** Each update makes one emulator step (`paleoemu_update`):

- **time:** the instance's evaluation time, which is the forcing time with
  `clock = "forcing"`;
- **CO2:** the coefficient `coef_co2`;
- **orbit:** the forcing's elements, which are also the surface's; without them, the
  emulator's own table (Laskar et al. 2004) at that time, as when PaleoEmulator runs
  standalone;
- **ice:** the geometry on the emulator's ice grid, against the reference.

**Adjustments.** The temperature can be adjusted before it is applied
(`paleoemu_downscale`), not its standard deviation:

- `scale_dtas` (default 1) multiplies the emulator's anomaly: a colder or warmer
  climate with CESM-3Ma's pattern, within the emulator's training states (unlike a CO2
  outside the range CESM-3Ma saw at that ice, 185-236 ppm at LGM ice).
- `bias` (default "", none) adds a bias correction, `bias_scale` (default 1) times the
  NH ice-volume weight `w` = (dV North America + Eurasia) / `dV_ref`: it follows the
  glacial state and vanishes at present. PaleoEmulator's
  `scripts/prep_bias_pmip4.jl` makes it from PMIP4: per model, the LGM minus piControl
  temperature at the reference surface (the lapse effect of its ice sheets removed)
  minus the emulator's 21 ka state, in `kryos_data/paleoemu/bias/`
  (`bias_pmip4_<model>.nc`, and `bias_pmip4_mean.nc` with the spread).
- `beta_adj` (default 0.05, 5 %/K) carries both into the precipitation, exp(beta_adj ΔT).

**Filter start.** The emulator also sees CO2 and the regional ice volumes through
exponential filters (3 and 10 kyr), its memory; they advance over the elapsed time.
At a cold start (the first update, PASS_BOUNDARY) the key `filter_start` sets them:

- `"training"` (default): the run's input plus the lag CESM-3Ma's trajectory had at
  that time, f(t0) = x(t0) + f_train(t0) − x_train(t0), from the emulator file's
  training inputs (PaleoEmulator's `paleoemu_init`). At 0 ka the 10-kyr filters still
  hold the deglaciation (North America +21.2 m SLE, Eurasia +3.9 m SLE, log(CO2/280)
  −0.11): a run at the reference state (the reference ice, CESM-3Ma's 280.4 ppm, the
  CO2 series' value at 0) reproduces the reference climate, `dtas`, `dlogpr`, `dto` 0.
  Outside 0-3 Ma the filters start in equilibrium.
- `"equilibrium"`: the run's input, as if it had been constant forever. At present day
  this lacks the deglacial memory: `dtas` −0.56 K over Greenland's ice and −0.58 K over
  Antarctica's (annual means; JJA +0.23 and −0.45 K), `dlogpr` +0.8 %, `dto` −2.8 to
  +1.5 K.

A restart restores the filters of the bundle, whatever `filter_start`.

Its anomalies are applied to the pre-industrial
reference climate on the grid (key `reference`), with the same formulas as
`paleoemu_downscale`.

**Exports:**

| Export | Units | |
|---|---|---|
| `tas` | K | at the reference surface `z_srf_clim` |
| `pr` | kg m-2 yr-1 | the emulator's mm d-1 × the days of the run's year |
| `z_srf_clim` | m | the surface of the reference climate (ERA5) |
| `lapse` | K m-1 | Γ = 6.5e-3, the emulator's lapse rate |
| `beta_p` | K-1 | key `beta_p`, default 0: precipitation does not change with the lapse rows' elevation, as in the emulator's own downscaling. Its `dlogpr` already holds CESM's response to the ice sheets; only `dtas` is at the reference surface |
| `dto_lev` | K | the annual ocean temperature change on the emulator's levels (axis `level`), filled and interpolated ([The ocean](#the-ocean)) |

**Reference climate.** It is a file on the instance grid holding `tas_ref` [K] and
`pr_ref` [mm d-1] by month, `z_srf_ref` [m] and the global attribute `grid_name`.
PaleoEmulator writes it from ERA5 1981-2010 minus the HadCRUT 5.1 warming since
1850-1900, corrected over the ice sheets by regional climate models of the same
period (MAR v3.14 for Greenland, RACMO2.3p2 for Antarctica), remapped by FesmData
onto the same grid. ERA5 alone is 1-2 K too cold in summer over Greenland's ablation
zone (600-1600 m) and too snowy at its margins. Each RCM replaces ERA5 over its domain
(its temperature taken to `z_srf_ref` by Γ), and its correction decays to 0 over 200 km
beyond it (`taper_km`):

```bash
scripts/prep_era5_clim.sh <dir> <kryos>/maps/grid_GLOBAL-0.5DEG.txt GLOBAL-0.5DEG
julia fesmdata.jl remap $FESMDATA_PREPARED/Fettweis2017_mar314/MARv3.14_ERA5_1981-2010.nc GLOBAL-0.5DEG --vars=t2m,pr,z_srf --name=SMB-MARv3.14-1981-2010
julia fesmdata.jl remap $FESMDATA_PREPARED/vanWessem2023_racmo23/RACMO2.3p2_1981-2010.nc GLOBAL-0.5DEG --vars=t2m,pr,z_srf --name=SMB-RACMO2.3p2-1981-2010
julia --project=scripts scripts/prep_kryos_reference.jl <dir> GLOBAL-0.5DEG <out.nc> \
    GLOBAL-0.5DEG_SMB-MARv3.14-1981-2010.nc GLOBAL-0.5DEG_SMB-RACMO2.3p2-1981-2010.nc
```

The par files use `kryos_data/paleoemu/paleoemu_reference_rcm_GLOBAL-0.5DEG.nc` (the
remapped RCMs in `kryos_data/paleoemu/rcm/`); `paleoemu_reference_GLOBAL-0.5DEG.nc`
is ERA5 alone. The present-day run (`par/kryos_GRL-ANT_data_paleoemu_pd.nml`) uses
`paleoemu_reference_rcm_1981-2010_GLOBAL-0.5DEG.nc`, made with `shift_pi=false`:
1981-2010 itself, to compare with MAR and RACMO of the same period.

**Restart.** The file `<name>_restart.nc` holds the emulator's memory: the filter
states and the time of the last update. On restore, the reference is set again from
the bundle's imports, and then the filters are restored.

**Output.**
- Streams `2D` and `1D`: the default writer, with the exports on `lon` and `lat`;
  `dto_lev` on `lon`, `lat` and the axis `level`.
- Stream `emu` (`out_emu`, `<name>_emu.nc`): the emulator's state on the CESM grid,
  as written by `paleoemu_write_step`. It holds the inputs, the regional volume
  changes, the elevation change, `dtas`, `dlogpr` and the ocean's `dto` (unfilled,
  0 outside the ocean), and their standard deviations.

## The regional port

The port is stateless (no backend restart file) and steps on the cadence and clock of
the global part. The rows `climate_g -> climate` remap `tas`, `pr`, `z_srf_clim`,
`lapse` and `beta_p` onto its grid: bilinear for the intensive fields, conservative
for `pr`. The port passes them on with their validity (and, with an ocean, makes the
profiles: [The ocean](#the-ocean)). The rows `climate -> surface`
then take `tas` and `pr` to the surface's elevation:

- `tas − lapse · (z_srf − z_srf_clim)`;
- `pr · exp(−beta_p · lapse · (z_srf − z_srf_clim))`.

## The ocean

The emulator predicts the annual ocean potential temperature change `dto` against
pre-industrial at 8 levels (`ocean_depth`, 4-1475 m) on the CESM grid, with its
standard deviation `dto_sd`. Its values are 0 on land and below the sea floor
(`ocean_mask`, fixed: CESM-3Ma's present-day bathymetry).

**The global part** always predicts and exports it, on its own axis `level` (the
emulator's depths, m):

1. `d = scale_dto · dto + dto_nsd · dto_sd` (keys, defaults 1 and 0);
2. `ocean_fill`: on each level, the cells outside the ocean are filled from the ocean,
   a Poisson fill periodic in longitude whose first guess is the zonal mean of the ocean
   (fesm-utils' `fill_poisson`, method 1); the ocean cells keep their values; a level
   without ocean takes the level above;
3. `dto_lev`: each level interpolated bilinearly onto the grid
   (`paleoemu_target_bilinear`). At `PASS_INIT` it is 0.

The fill matters on Antarctica's shelves: the CESM grid (~3.75°) sees little ocean
there, and the shelves of a regional grid lie mostly in CESM land cells.

At present day `dto` is 0 with the default filter start (above): the emulator's
reference is CESM-3Ma at 0 ka, whose filters hold the deglaciation, and a run started
there inherits that memory. A transient that reaches 0 ka carries its own memory
instead, and its `dto` there is the emulator's response to it.

**The port** takes the ocean only with a present-day reference on its grid, the keys
`to_ref` and `so_ref`: var_defs groups read as the snapshots read a snapshot's `to`
and `so` (units, scaling, `fill`, the annual mean of a monthly climatology), e.g.
`ismip7_to`, `ismip7_so` of the ISMIP7 climatology. Without them (the default) the
port has no ocean, as before. With them it sets `caps%has_ocn`, imports `dto_lev`
on the axis `level` (the domains see the global part's axes: `domain_create`'s
`axes_g`) through the row `climate_g -> climate dto_lev` (bilinear, level by level),
and exports the profiles of set A on the file's depths (axis `depth`, |z|):

- `dto_ann = W · dto_lev`, `W` linear in depth between the levels, constant above
  the top level and below the bottom one;
- `to_ann = to_ref + dto_ann`, `so_ann = so_ref`.

The shelf (e.g. marine_shelf, set A) interpolates them to its shelf base as with the
snapshots climate.

**Basin means.** With `dto_basins` (a basin set of the domain's regions, e.g.
`Zwally2012` as the shelf corrections use), `dto_lev` is replaced in each basin by
its mean, level by level, over the basin's sample `dto_sample` (a regions selection,
default `zone:shelf_break_buffer`: FesmData's zone between the continental shelf and
the open ocean). The basin sets end at the continental shelf, so the basins are first
extended over the whole domain by nearest neighbour (`fill_nearest`): the shelf-break
buffer, and the shelf the ice may advance onto, belong to the nearest basin. The shelf
then sees the change of the ocean it exchanges with, uniform within the basin, rather
than the filled values on the shelf itself. A basin without sampled cells (e.g. a basin
whose shelf-break buffer is beyond the grid) samples its ocean cells instead: those of
the reference topography of the domain on the port's grid (the port imports no geometry):
bed below sea level and not grounded, `rho_ice H_ice < −rho_sw z_bed`. A basin without
ocean either keeps its cells' own values. Both are logged at init. The sample is static
(the present-day zones and topography), as is the emulator's ocean mask.

**Caveats.**

- `dto_lev` is the change against the emulator's reference (CESM-3Ma at 0 ka:
  pre-industrial CO2, today's ice), while ISMIP7 (and ORAS5) are present-day
  climatologies: the profiles add that change to a present-day ocean, so the warming
  since pre-industrial stays in the reference. Accepted: the shelf's `tf_corr` are
  calibrated at present day on that same reference.
- Drafts below 1475 m take the change at 1475 m (constant below the bottom level).
- No sea-level datum shift: the depths are those of the present-day reference, not
  shifted with the glacial sea level.
- Salinity stays at present day (`so_ref`).

**Output.** The port's stream `2D` (the default writer) holds `to_ann`, `so_ann` and
`dto_ann` on `xc`, `yc` and its axis `depth`; the global part's holds `dto_lev`.

## Shelf diagnostics

The shelf's response to the emulator's ocean runs on imposed geometry only, with no ice
model and no GIA: the ice backend `data` (PaleoMIST's `H_ice` and `z_bed`, the bed
relative to the sea level of the time, `z_sl = 0`; or the reference topography) and a
shelf instance (`marine_shelf`) per domain, which imports `H_ice`, `z_bed`, `f_grnd`,
`z_sl`, `z_srf` from the hub and set A from the port. The data ice ignores `bmb_shlf`
(its import is optional): the shelf is a diagnostic. Its stream `2D` holds `bmb_shlf`
[m/yr], `T_shlf`, `tf_shlf` and `dT_shlf` (the port's `dto_ann` at the shelf base). The
emulator's levels are depths below the sea surface of each period, the frame of the
data ice. `par/kryos_NH-ANT_data_paleoemu.nml` and
`par/kryos_GRL-ANT_data_paleoemu_pd.nml` hold them: Antarctica as
`par/kryos_NH-ANT.nml`'s shelf (quad-nl, `tf_method` 1, the Zwally2012 tf corrections)
on ISMIP7 with Zwally2012 basin means; North and Greenland `anom` (`tf_method` 2) on
ORAS5 1981-2010 (the annual mean of the monthly climatology).

## Sequence

**Cold start.** The global part runs before the domains' climates, at `PASS_BOUNDARY`.
The domains' ice is not known yet, so its geometry is the background alone at t0
(`globe_ice_external`). The ports then take its atmosphere, and the surfaces spin up
on it.

**Step.** At the start of a step:
1. the gather of the previous step's ice, made whenever the global climate (or a
   global geo) is due;
2. the forcing;
3. the global climate, on that gather's merged ice and this step's CO2;
4. in each domain, after the ice and the hub: the rows into the port, the port, then
   the surface.

## One orbit

The forcing resolves the run's orbit once, from the group named by `[forcing] orbit`
(`kryos_orbit`, insol's tables), at the forcing time or fixed. It exports the elements
`orbit_ecc`, `orbit_obl` and `orbit_per`, and every component that takes the orbit reads
them through the rows `forcing -> climate_g / climate / surface`:

- the global climate passes them to `paleoemu_update`;
- chion's insolation (itm, bessi) computes its daily S from them on its own latitudes,
  with an insol state without tables.

The climate and the insolation therefore always see one orbit. The emulator's own
perihelion convention is Berger's plus 180°: `paleoemu_update` converts the elements.
The forcing's stream `insol` (`out_insol`, `<instance>_insol.nc`) writes the elements and
the monthly mean top-of-atmosphere insolation on 1° latitudes, as a diagnostic.

## Limits

- **Antarctica.** Its ice was fixed in CESM-3Ma, so Antarctic ice changes do not feed
  back. The emulator's ice grid covers 21-85.5°N only.
- **The training ice.** The emulator learnt the response to the ice of Willeit et al.
  (2019), whose 0 ka state holds a Scandinavian ice cap (3.2 m SLE).
  - kryos sets the reference to its own present day: the spliced present day of
    `globe_ice`.
  - The volume and elevation changes are therefore measured against kryos's present
    day.
- **Cold LIG and mid-Holocene.** The emulator inherits CESM-3Ma's: the run's 5×
  acceleration, the 1.5× GHG scaling (see PaleoEmulator's PMIP4 comparison).
