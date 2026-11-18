---
title: "Par files"
---

A kryos run reads one par file, a Fortran namelist file passed as the only argument of `kryos.x`. The shipped
files are in [par/](https://github.com/fesmc/kryos/tree/main/par); the gate files in `par/gate/`, the friction
calibration in `par/calib/`. Every file starts with a comment header that says what it sets up, what it was
derived from and the build it needs.

## The plastic and VILMA2 runs

Hand-written; keys left out take the templates of `input/kryos_defaults.nml`. All run 80 ka to present
(`time_init = -80e3`, `time_end = 0`, `dtt = 100` yr), with the PaleoMIST margins (`f_ice > 0.5`, signed-distance
interpolation) and Gowan's time-dependent tau_b domains of FesmData IceHistory. Build: `plastic=1`, and `vilma=1`
where geo is `vilma2` (`yelmo=0` is enough). See [Online ice](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/howto/online-ice.md) and
[PlasticIceSheet](backend-plastic.md). The standalone VILMA2 runs of imposed ice ([Imposed ice](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/howto/imposed-ice.md))
read only the `&vilma` group of these files: their Earth. Status as of 2026-11-16.

| Par file | Domains (hub) | Ice | Geo | What it sets up |
|---|---|---|---|---|
| [kryos_North_plastic.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic.nml) | North (NH-32KM) | plastic | none | Smoke run on a fixed bed: today's bed (TOPO-GEBCO2026), `z_sl = 0`, no globe, no climate |
| [kryos_North_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2.nml) | North (NH-32KM) | plastic | vilma2 | The global VILMA2 Earth in 3-D (Bagge structure 6), reference bed today's + r (`kryos_data/refcorr/PaleoMIST-sci`), the ice outside NH from the PaleoMIST history; Earth cadence 100 yr |
| [kryos_North_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2_b1.nml) | North (NH-32KM) | plastic | vilma2 | B1's Earth and reference: Gowan's 1-D Earth (120 km lithosphere, 4e20 / 4e22 Pa s, 12 layers, lmax 170), the relaxed present-day reference (no r) |
| [kryos_NH-ANT_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2_b1.nml) | North (NH-32KM), Antarctica (ANT-16KM) | plastic, plastic | vilma2 | As `_b1` with a second plastic instance on ANT-16KM (BedMachine Antarctica v4, its own geo port). Its `&vilma` group is `E1D-PMIST` (I1, I5) and the 1-D template of the Earth search |
| [kryos_NH-ANT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2.nml) | North, Antarctica | plastic, plastic | vilma2 | Both hemispheres on Bagge's 3-D Earth: PREM 11 layers, 80 km lithosphere, Bagge et al. (2021) 3-D viscosity (structure 6), lmax 170, 1-D pre-spin-up, 2e4 yr cap on the 3-D relaxation, refcorr reference |
| [kryos_NH-ANT-PAT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2.nml) | North, Antarctica, PatagoniaPaleo (PAT-PAL-8KM) | plastic x 3 | vilma2 | As the previous with a third plastic instance on PAT-PAL-8KM (Gowan's region 4 domains, IceHistory v2.4.0, GEBCO2026 + r); the maps folder needs `grid_PAT-PAL-8KM.txt`. The par file of ens-v1 and ens-v2, `E3D-B21-6` of the imposed-ice runs, and the 3-D template of the Earth search |
| [kryos_NH-ANT-PAT_plastic_vilma2_nh16.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_nh16.nml) | North (NH-16KM), Antarctica, PatagoniaPaleo | plastic x 3 | vilma2 | As the previous with the Northern Hemisphere on NH-16KM (D3; needs `grid_NH-16KM.txt` and the case's NH-16KM reference topography) |
| [kryos_NH-ANT-PAT_plastic_vilma2_bagge1d.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_bagge1d.nml), [`_bagge1d_layers.nml`](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_bagge1d_layers.nml) | North, Antarctica, PatagoniaPaleo | plastic x 3 | vilma2 | The radial mean of Bagge's structure 6 as a 1-D Earth: VILMA2's 1Db (`visc_3d_tol = 99`), and the same as plain layer means (`E1D-B21-6`, `E1D-B21-6-layers`) |
| [kryos_NH-ANT-PAT_plastic_vilma2_inv.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml) | North, Antarctica, PatagoniaPaleo | plastic x 3 | vilma2 | The implied mass balance by joint inversion (`mass_balance = "inverted"`, `slide_law = "coulomb_n"`) with a surface instance `smb_simple` per domain as the prior (no climate: `tsl_const`); see [Implied mass balance](implied-mb.md) |

In the `_b1` files and the NH-ANT(-PAT) files, `taub_domain_factor_file` (the present-day factors per domain,
Stage 1) is `"kryos_data/plastic/taub_domain_factor_v3.nc"` on every instance (the released factors; v3 is v2 plus
Patagonia, identical on all NH and ANT domains), and `taub_scale_file` (the multipliers per scale group of an
ensemble member) is `""`. `kryos_North_plastic.nml` and `kryos_North_plastic_vilma2.nml` list neither key (the default
`""`: Gowan's tau_b as it is). The ensemble tools set them per member ([Ensembles](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/howto/ensembles.md)); for one run,
`-p ice_nor.taub_scale_file=<file> ...`.

**Generated par files of the Earth search.** `plastic_ensemble.jl earth-search-sample`, `earth-search-sample3d` and
`earth-search-sample3d-box` write one par file per sampled Earth from the two templates (`_b1` for 1-D, `NH-ANT-PAT` for
3-D): `earths/E<nnn>.nml` (wave A, 1-D), `earthsB/E<nnn>.nml` with its field `E<nnn>.nc` and `E<nnn>_1d.nml` (wave B, 3-D and
its radial profile), and the same for wave B2 in `waveB2/earths/` (with `visc_log10_min = 19.0`). They live in the search's
work folder, not in `par/` ([The Earth search](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/howto/earth-search.md)). `ice-earth-recon/experiments/refs-ice.sh` writes the Bagge structures
1-8 likewise (`$KRYOS_WORK/experiments/par/E3D-B21-<k>.nml`: the 3-D par file with `visc_3d_structure = k`).

## The Yelmo configurations

Converted from yelmox or derived from the converted files; generated, so do not edit them: change the source and
generate again.

| Par file | Domains | Generated by | What it sets up |
|---|---|---|---|
| [kryos_Greenland.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_Greenland.nml) | Greenland (GRL-16KM) | `tools/bin/yelmox2kryos.jl` from `yelmox_Greenland.nml` | yelmox's Greenland: Yelmo, chion, FastIsostasy, 15 kyr, the yelmox output schedules |
| [kryos_Antarctica.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_Antarctica.nml) | Antarctica (hub ANT-16KM, ice ANT-32KM) | `yelmox2kryos.jl` from `yelmox_Antarctica.nml`, `-p sim.relax=False sim.opt=False` | yelmox's Antarctica without relaxation and optimisation, 15 kyr |
| [kryos_North.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North.nml) | North (NH-32KM) | `yelmox2kryos.jl` from `yelmox_North.nml` | yelmox's North: 100 kyr at a constant LGM forcing time (`time_const = -21e3`) |
| [kryos_NH-ANT.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT.nml) | North, Antarctica | `tools/bin/kryos-merge.jl` from `kryos_North.nml` + `kryos_Antarctica.nml` | The two domains in one run with FastIsostasy (15 kyr, forcing time 0); merge report `kryos_NH-ANT_merge.md` |
| [kryos_NH-ANT_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_vilma2.nml) | North, Antarctica | `tools/bin/kryos-vilma2-par.jl` from `kryos_NH-ANT.nml` | The science case: Yelmo in both domains coupled to the global VILMA2 Earth (lmax 170, `dt = 10`), 80 ka transient from PaleoMIST, refcorr reference, bundles every 5 kyr; every choice in `kryos_NH-ANT_vilma2_conversion.md` |
| [kryos_NH-ANT_calib.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_calib.nml) | North, Antarctica | `tools/bin/kryos-calib-par.jl sci` from `kryos_NH-ANT.nml` | The calib-cf variant: `[ytill_ant] scale_zb` 0 -> 2, `cf_ref` 0.1 -> 0.05, threshold calving in `[ycalv_ant]` |
| [kryos_NH-ANT_vilma2_calib.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_vilma2_calib.nml) | North, Antarctica | `kryos-calib-par.jl sci` from `kryos_NH-ANT_vilma2.nml` | The same changes on the science case |

Each converted file has its report next to it (`par/kryos_<domain>_conversion.md`): every yelmox key mapped,
default-equal or dropped, with the reason. The converter needs the Julia project of `tools/`; see
[Tools](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/reference-tools.md).

## par/gate/

The files of the regression gates (yelmox against kryos, phase by phase; [Developer notes](developer-notes.md)).
The yelmox conversions are 1000-yr runs, lossless (`write_lossless = True`), with dense output schedules.

| Par file | What it sets up |
|---|---|
| `grl.nml`, `nor.nml` | `yelmox_Greenland`, `yelmox_North`, 1000 yr |
| `ant.nml` | `yelmox_Antarctica` without relaxation and optimisation, 1000 yr |
| `ant-cf.nml` | `yelmox_Antarctica` with `opt.opt_tf = none` (the interim gate of phase 2) |
| `ant-full.nml` | `yelmox_Antarctica` with relaxation and optimisation (the phase 2 gate) |
| `sci.nml` | `nor.nml` + `ant.nml` merged (`kryos-merge.jl`) |
| `sci-hub.nml` | `sci.nml` with the global geometry hub as an observer (`[globe_ice]` appended, gate G1b, `tools/gate/make_sci_hub.sh`); must not change any result of `sci.nml` |
| `g5-vilma2.nml`, `g5-fastiso.nml` | Gate G5: Antarctica the only active domain, from PaleoMIST at -20e3 yr, paleo forcing (`lgm` + `ocean`), VILMA2 or FastIsostasy (`tools/bin/kryos-ant-par.jl`) |
| `g5-vilma2-calib.nml`, `g5-fastiso-calib.nml` | the same from the calib-cf variants |
| `g6-smb.nml` | Gate G6: `kryos_NH-ANT_vilma2.nml` with the surface of both domains from `smb_simple` (`tools/bin/kryos-smb-par.jl`) |

## par/calib/

The single-domain runs of the basal-friction (cf_ref) calibration, generated by `tools/bin/kryos-calib-par.jl`
from `kryos_NH-ANT_calib.nml`, with FastIsostasy and frozen forcing
([docs/dev/calib-cf.md](https://github.com/fesmc/kryos/blob/main/docs/dev/calib-cf.md); submitted with
`tools/levante/archive/calib.sh`):

| Par file | Case |
|---|---|
| `Antarctica-pd.nml`, `North-pd.nml` | pd: 0 -> 30 kyr, present-day forcing, cold start from the hub reference |
| `Antarctica-lgm.nml`, `North-lgm.nml` | lgm: 30 -> 60 kyr, restarted from pd, LGM forcing (alpha 1, sea level -120 m, orbit -21 ka) |
| `North-lgmrecon.nml` | lgmrecon: LGM forcing, cold start from PaleoMIST at 20 ka, 0 -> 20 kyr |

## How a par file is structured

A par file is a set of namelist groups. The run group `[kryos]` names the rest:

```
&kryos      domains = "North", globe = "globe", clock = "clock", restart_out = "out_rst", out_dir = "./", ...
&clock      time_init, time_end, dtt, method, ...
&out_1D     method = "const", dt = 100.0           schedule groups, named by the out_* keys
&out_2D     ...
&globe      geo = "geo_glob"                       the global instances (forcing, geo)
&geo_glob   backend = "vilma2", dt, out_2D, vilma = "vilma", ice_bg = "globe_ice", ...
&vilma      overrides of input/vilma_defaults.nml
&globe_ice  the global ice field; its field groups &X_z_bed, &X_H_ice, &X_dH_ice
&North      grid_hub = "NH-32KM", topo = "topo_nor", regions = "regions_nor", ice = "ice_nor", geo = "geo_nor"
&ice_nor    backend = "plastic", grid, dt, out_2D, out_1D, ...  and the [ice_plastic] keys
&geo_nor    backend = "vilma2", grid = "NH-32KM", ...  the regional port
&topo_nor   path = "fesmdata/{domain}/{grid_name}/{grid_name}_TOPO-GEBCO2026.nc", ...
&regions_nor
```

- **Domains.** `[kryos] domains` lists the domain groups. A domain group names its hub grid and one instance group
  per component (`ice`, `geo`, `climate`, `surface`, `shelf`, `geothermal`, `sediments`; `""` = absent): a component
  runs if and only if the domain names an instance. Instance names are unique in a run (`ice_nor`, `ice_ant`).
- **Instances.** An instance group holds the core keys of `[instance]` (`backend`, `grid` (`""` = the default chain:
  ice on the hub, geo on the ice, ...), `dt` (0: every step), `phase`, `clock`, `out_2D`, `out_1D`, `out_vars`) and
  the keys of its backend's template `[<component>_<backend>]` (e.g. `[ice_plastic]`, `[geo_vilma2]`,
  `[surface_smb_simple]`).
- **Schedules.** The `out_*` and `restart_out` keys name schedule groups (`[out_2D]`, `[out_rst]`, ...) of the
  template `[schedule]` (`method = none | const | times | file`); `"none"` is off.
- **Margin patches.** A plastic instance's `margin_patches` names groups of the template `[margin_patch]`, each a
  regional margin chronology (file, region map, region ids; [backend-plastic.md](backend-plastic.md#margin-patches)).
- **Paths.** Relative to the run directory; `{domain}` and `{grid_name}` are replaced (e.g. the margin file
  `fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc`).

### Templates and defaults

- The kryos groups (`[kryos]`, the clock, schedules, `[globe]`, domains, instances) are validated against the
  templates in [input/kryos_defaults.nml](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml) (one
  per group; an instance against `[instance]` plus `[<component>_<backend>]`): a key left out takes the template's
  value, and a key marked "required" there must be set.
- Backend library groups (`[yelmo]`, `[ytopo]`, ..., `[chion]`) take their defaults from the library's defaults file
  in `input/` (`yelmo_defaults.nml`, `chion_defaults.nml`). `[vilma]` takes its defaults from
  `input/vilma_defaults.nml`: the group sets overrides only, and its name is fixed (VILMA2 looks the defaults up under
  it). The groups of code kryos hosts (`[smb_simple]`, the chion host's `[orbit]`, `[globe_ice]`) have their schema at
  the end of `kryos_defaults.nml`.
- Data field groups (`[X_z_bed]`, ...) follow `input/datasrc_defaults.nml`; snapshot, term and coefficient groups
  `snapshots_defaults.nml` and `forcing_defaults.nml`.
- Groups whose library has no defaults file (`[topo]`, `[regions]`, `[isos]`, `[barysealevel]`, `[marine_shelf]`)
  have no template: set every key the library reads.
- Named library groups (`snapshot_*`, `term_*`, `coef_*`, data field groups) are looked up in the par file, then in
  the files of `[kryos] library` (kryos ships none).

So a hand-written par file lists only what differs from the defaults (the plastic files). The converted files are
fully expanded instead: every key of every group is set, `library = ""`, and the header records the sha256 of the
input files they were made against, so they stay valid when a default changes.

### Overrides

- **runme**: `-p group.key=value` overrides a key of the par file for one run (`bsl` is an alias of
  `barysealevel`); a comma list makes an ensemble: `-p ydyn.beta_q=0.2,0.5,1.0`. The key must exist in the par file
  or in `input/kryos_defaults.nml`. The clock group is the one `[kryos] clock` names
  ([Getting started](getting-started.md)).
- **[tools/levante/set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py)**: sets existing keys of a
  fully expanded par file in place (`group.key=value`, case-insensitive names; quote strings yourself; a trailing comment
  is kept; a key the group does not have stops it: it never adds keys).

```bash
python3 tools/levante/set_nml.py par.nml clock.time_end=-70e3 'kryos.restart_out="none"' 'ice_nor.taub_scale_file="taub_scale.nc"'
```

- `kryos.x` itself takes no overrides: one argument, the par file.

Do not edit a generated file (the converted, merged and derived files above): change the source or the generator
and generate again (Julia project `tools/`, [Tools](https://github.com/awi-esc/ice-earth-recon/blob/main/docs/reference-tools.md)).
