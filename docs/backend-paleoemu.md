---
title: "PaleoEmulator (climate)"
---

The climate backend `paleoemu` couples the CESM-3Ma climate emulator of
[PaleoEmulator](https://github.com/fesmc/PaleoEmulator) to kryos. The emulator
returns monthly 2-m temperature and precipitation anomalies against pre-industrial.
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
drive the ice.

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
- **Par file:** `par/kryos_North_plastic_paleoemu.nml`: the PaleoMIST plastic ice of
  North with chion, from 21 to 20 ka.
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

**Scaled anomalies.** The keys `scale_dtas` and `scale_dlogpr` (default 1) multiply the
emulator's anomalies before they are applied (`paleoemu_downscale(..., scale)`), not
their standard deviations. Equal values give a colder or warmer climate with CESM-3Ma's
pattern, a first-order sensitivity dimension of an ensemble. Unlike a CO2 outside the
range CESM-3Ma saw at that ice (185-236 ppm at LGM ice), it stays within the
emulator's training states.

The emulator's filters of CO2 and ice volume start in equilibrium at the first update
and advance over the elapsed time. Its anomalies are applied to the pre-industrial
reference climate on the grid (key `reference`), with the same formulas as
`paleoemu_downscale`.

**Exports:**

| Export | Units | |
|---|---|---|
| `tas` | K | at the reference surface `z_srf_clim` |
| `pr` | kg m-2 yr-1 | the emulator's mm d-1 × the days of the run's year |
| `z_srf_clim` | m | the surface of the reference climate (ERA5) |
| `lapse` | K m-1 | Γ = 6.5e-3, the emulator's lapse rate |
| `beta_p` | K-1 | key `beta_p`, default 0: precipitation does not change with elevation, as in the emulator's own downscaling |

**Reference climate.** It is a file on the instance grid holding `tas_ref` [K] and
`pr_ref` [mm d-1] by month, `z_srf_ref` [m] and the global attribute `grid_name`.
PaleoEmulator writes it from ERA5 1981-2010 minus the HadCRUT 5.1 warming since
1850-1900:

```bash
scripts/prep_era5_clim.sh <dir> <kryos>/maps/grid_GLOBAL-0.5DEG.txt GLOBAL-0.5DEG
julia --project=scripts scripts/prep_kryos_reference.jl <dir> GLOBAL-0.5DEG <out.nc>
```

**Restart.** The file `<name>_restart.nc` holds the emulator's memory: the filter
states and the time of the last update. On restore, the reference is set again from
the bundle's imports, and then the filters are restored.

**Output.**
- Streams `2D` and `1D`: the default writer, with the exports on `lon` and `lat`.
- Stream `emu` (`out_emu`, `<name>_emu.nc`): the emulator's state on the CESM grid,
  as written by `paleoemu_write_step`. It holds the inputs, the regional volume
  changes, the elevation change, `dtas` and `dlogpr`, and their standard deviations.

## The regional port

The port is stateless (no backend restart file) and steps on the cadence and clock of
the global part. The rows `climate_g -> climate` remap `tas`, `pr`, `z_srf_clim`,
`lapse` and `beta_p` onto its grid: bilinear for the intensive fields, conservative
for `pr`. The port passes them on with their validity. The rows `climate -> surface`
then take `tas` and `pr` to the surface's elevation:

- `tas − lapse · (z_srf − z_srf_clim)`;
- `pr · exp(−beta_p · lapse · (z_srf − z_srf_clim))`.

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
