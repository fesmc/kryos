---
title: "data (ice)"
---

The ice backend `data` reads the ice geometry instead of solving it. The thickness,
the bed and the sea surface come from field groups
([`input/datasrc_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/datasrc_defaults.nml)),
static or transient, e.g. an ice history such as PaleoMIST. A key left empty takes
the domain's reference topography (the sea surface: 0), so an instance without keys
is the present day. Every component downstream sees the fields as an ice model's
exports: the domain's hub, the climate and the surface, and `globe_ice`, which merges
every domain's ice into the background of the globe.

Use it to drive a climate and a surface on an imposed ice sheet, e.g. to tune the
surface mass balance on PaleoMIST without solving the ice.

- **Build:** always compiled; no library.
- **Select:** `backend = "data"` in the ice instance group of a domain; the keys are
  those of the template `[ice_data]` in
  [`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml).
- **Source:** [`src/backends/ice/kryos_ice_data.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/kryos_ice_data.f90).
- **Test:** `tests/test_ice_data.f90`.

## Keys

| Key | | Default |
|---|---|---|
| `H_ice` | field group of the ice thickness [m] | `""`: the domain's reference topography |
| `z_bed` | field group of the bed [m] | `""`: the domain's reference topography |
| `z_sl` | field group of the sea surface [m] | `""`: 0 |

Each key names a field group of the par file (`filename`, `name`, `mode`, ...). Each
domain has its own instance, so domains can read different sources, and they can mix
`data` with `plastic` or `yelmo`. Times on the file follow its time units
(`"years ago"` counts from 1950 through the run clock's reference year), or
`time_units = "clock"` for model years.

```
&ice_grl
    backend = "data"
    H_ice   = "grl_H"
    z_bed   = "grl_zb"
    z_sl    = "grl_sl"
/
&grl_H
    filename = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"
    name     = "H_ice"
    mode     = "transient"
/
```

## Update

At each update the fields are read at the instance's evaluation time (the forcing
time with `clock = "forcing"`), linear in time between the records. The thickness is
floored at 0.

**Exports**, as `plastic`'s:

| Export | |
|---|---|
| `H_ice`, `z_bed`, `z_sl` | the data |
| `z_srf` | `z_bed + H` grounded; `z_sl + (1 − ρ_i/ρ_sw) H` floating; `max(z_bed, z_sl)` without ice |
| `f_ice`, `f_grnd`, `H_grnd` | the ice mask; grounded where ρ_i H ≥ ρ_sw (z_sl − z_bed), with the run's densities |
| `dHidt`, `mb_net` | the change of the data over the step |
| `smb`, `bmb`, `cmb` | 0 |
| `z_bed_ref`, `H_ice_pd`, `footprint` | static: the reference topography; the whole domain |

**Restart.** The backend is stateless: no restart file. The fields are read again at
the restored time.
