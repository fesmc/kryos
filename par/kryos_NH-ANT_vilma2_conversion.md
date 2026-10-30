# par/kryos_NH-ANT_vilma2.nml: the science case with geo = vilma2

Generated from `par/kryos_NH-ANT.nml` (the merged NH + ANT file, FastIsostasy; regenerated for this: the committed file predated the shelf move (1d6305b) and lacked `[shelf_corr_ant]`, so `shelf_ant` stopped at init) by `tools/bin/kryos-vilma2-par.jl`
(a text edit by group; nothing else changes, every other key and comment is the base's). Phase 4 WP8 (part, started in WP7:
the file the G3 and G4 runs use). Every choice below is open to review; the ones marked **(choice)** are mine, not decided
in a plan.

## What the file is

NH (`yelmox_North`, ice grid NH-32KM) and ANT (`yelmox_Antarctica`, ice grid ANT-32KM, hub ANT-16KM, relax and opt off as in
the base) in one run, coupled to the global VILMA2 Earth at lmax 170 through two regional ports. Climate and surface as in
the base: snapshots + tsgen + chion. No `obm`.

## Changes against `kryos_NH-ANT.nml`

| Group | Change | Why |
|---|---|---|
| `[clock]` | `method = "rel"`, `time_init = -80e3`, `time_end = 0`, `time_const` removed | the transient clock: forcing time = model time; 80 ka is the length of PaleoMIST (plan Q4). The base is `method = "const"`, `time_const = 0` (frozen forcing). A shorter run: `clock.time_end`. |
| `[coef_alpha]` | `method = "series"`, `series_file = "input/alpha_combined_125kyr_interp.dat"`; `norm_t0`, `norm_t1` removed **(choice)** | a transient clock with the base's constant alpha = 1 (the LGM) would be an LGM climate for 80 kyr. The file is the glacial index of the yelmox transient setups (0 at 0 ka, about 1 at the LGM, 0.50 at -80 ka), already normalised: with `norm_t0`/`norm_t1` present in the group the series is normalised again and the run stops (both 0: "reference values are equal"). Back to the base: `coef_alpha.method="const"`. |
| `[geo_glob]` | `backend = "vilma2"`, `dt = 10`, `clock = "model"`, `vilma = "vilma"`, `ice_bg = "globe_ice"`, `allow_no_spinup = False`, `out_2D = "out_2D"`, `out_gia = out_vilma = "none"` | plan 3.4; `dt = 10` is Q6 (revisit after G7). The 2D output of the Earth on the Gauss grid is on at the schedule of the run; the streams `gia` (what VILMA2 got) and `vilma` (its own diagnostics) are off, the gates switch them on. |
| `[recon_nor]`, `[recon_ant]` | PaleoMIST at -80 ka, whole domain, `recon_zb_min = -1e5`; `[Antarctica] experiments = "recon_ant"` | decided (Alex), choice 1 below. |
| `[geo_nor]`, `[geo_ant]` | `backend = "vilma2"`, `dt = 10`, `grid` = the ice grid, the FastIsostasy keys removed | the ports (plan D16: on the ice grid; a `dt` that differs from the global's stops). |
| `[barysealevel]`, `[isos_nor]`, `[isos_ant]` | removed | FastIsostasy only. |
| `[vilma]` (new) | the CLIMBER-X set: `earth = "PREM"`, 11 layers (`r_bot`, `r_top`, `eta`, `rheology` of `input/vilma_defaults.nml`), `lmax = 170`, `l_visc_3d`, `visc_3d_file` = FesmData `GLOBAL-T170_SolidEarth-Bagge2021.nc`, `visc_3d_structure = 6`, `pre_spinup_1d = .true.`, `equil_time_max = 0.0` | the plan's Earth. Densities and `r_earth` are the run's (`phys_const`), not keys. **`equil_time_max = 0` is a placeholder**: the full-model relaxation after the 1-D pre-spin-up is G2c's question (Alex). With `pre_spinup_1d` alone the run has a spin-up (the `gia_check` condition) and holds the t0 load. |
| `[globe_ice]` + `[X_z_bed]`, `[X_H_ice]`, `[X_dH_ice]` (new) | background `GLOBAL-0.5DEG`, X = PaleoMIST at 0 ka (`z_bed`, `H_ice`), `ext = "file"` with the history's anomaly against its 0 ka slice, `ext_tref = 0`, `out_of_range = "error"` (the gather holds the history after `ext_tref` itself; the template requires it), `out_2D = "none"` | plan 3.4, Q3, Q4. The file is `fesmdata/Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc` (33 records, 80 to 0 ka, "years ago"). The description `grid_GLOBAL-0.5DEG.txt` is copied from FesmData into the map folder `maps/S6`. |

## Choices to review

1. **Start state at -80 ka (decided, Alex):** both domains start, each over the whole domain, from PaleoMIST at -80 ka: `[recon_nor]` and
   `[recon_ant]` read `fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc`, `recon_var = "H_ice"`, `recon_time = -80e3`
   (the first record of the "years ago" axis), `recon_regions = "all"`, `recon_zb_min = -1e5` (no bed filter: marine ice is kept); `[Antarctica]`
   `experiments = "recon_ant"` (the base had none). The VILMA2 reference stays present day (the ice's `z_bed_ref`, `H_ice_pd`). The base's NH recon
   (ICE-6G_C at -26 ka in North America) is replaced. The t0 load of each footprint then equals PaleoMIST at -80 ka (checked in the first short
   science run). `recon_smooth_sigma = 0` and `recon_H_float_min = 0` in both groups (the PaleoMIST start is already smoothed).
2. **Alpha series (choice):** above.
3. **`dt = 10` and the cadence of the ports** equal the base's `dt_isos`.
4. **Spin-up:** `pre_spinup_1d` only (above).
5. **Output:** `out_2D` of the Earth at the run's 2D schedule (the base's 200 yr; the gates use 50); `out_gia`, `out_vilma` off.

## The `smb_simple` variant (`-p`, not run yet)

Surface from `smb_simple` instead of chion, both domains: `surface_nor.backend="smb_simple"`, `surface_ant.backend="smb_simple"`
(the template `[surface_smb_simple]` gives `smb_simple = "smb_simple"`, the module defaults of `[smb_simple]`: not tuned).
The `surface_chion`, `chion`, `orbit` keys of the two groups stay (unused). G6 runs it; whether the unused keys are accepted is
checked there.

## Regenerate

```
julia tools/bin/kryos-vilma2-par.jl par/kryos_NH-ANT.nml par/kryos_NH-ANT_vilma2.nml
```
