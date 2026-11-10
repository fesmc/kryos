---
title: "smb_simple"
---

The surface backend `smb_simple` is the synthetic-elevation SMB scheme `syn` of yelmox
(`libs/smb_simple.f90`): the surface mass balance and surface temperature from a
**target ice mask**, the surface elevation and the sea-level air temperature alone,
with a latitude-dependent snow line. It needs no library and is stateless: every call
computes the fields from the imports, nothing is kept and no restart file is written.
The kernels are the yelmox expressions unchanged (a test compares them with the
yelmox module), computed in single precision.

- **Select:** `backend = "smb_simple"` in the surface instance group of a domain;
  instance key `smb_simple` (template `[surface_smb_simple]`, default
  `"smb_simple"`) names the parameter group, whose schema is the template
  `&smb_simple` of
  [`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml).
- **Source:** [`src/backends/surface/kryos_surface_smb_simple.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/surface/kryos_surface_smb_simple.f90)
  (the kernels `kryos_smb_simple` and the backend `kryos_surface_smb_simple`).
- **Imports:** `z_srf` (hub), `tsl_ann` (climate), `H_ice_ref` (ice), depending on the
  mask and `tsl_const` (below).
- **Exports:** `smb` [kg m-2 yr-1 = mm w.e. yr-1] and `T_srf` [K], means over the
  step; they enter the ice lagged by one step.

## Scheme

With φ = |latitude| (hemisphere-agnostic) and d the signed distance to the boundary of
the target mask (m, positive inside):

1. **Synthetic elevation** z_syn. With `use_plastic = True`, a plastic (Nye/Vialov)
   profile inside, $z = \min(z_\mathrm{max,in},\ C\sqrt{d})$,
   $C = \sqrt{2\tau_0/(\rho_i g)}$, and a linear ramp outside,
   $z = \max(-z_\mathrm{max,out},\ s_\mathrm{out}\, d)$; inside the mask
   z_syn = max(z, z_srf), outside min(z, z_srf). With `use_plastic = False`, a linear
   wedge of `slope` inside. ρ_i and g are the run's constants.
2. **Snow line** $z_{SL} = a_1 + a_2\varphi + a_3\,\mathrm{CO_2} + a_4 (f - f_\mathrm{mean})$,
   with the constants `co2_const`, `f_const`.
3. **SMB** = acc − abl, with accumulation $acc = c_0\, f_\mathrm{lat}(\varphi)\,
   e^{-d_\mathrm{in}/L_\mathrm{acc}}$ (largest at the margin, decaying inland; d_in in
   km, 0 outside), $f_\mathrm{lat} = \mathrm{clamp}(1 - k_\mathrm{acc,lat}
   \max(0, \varphi - \varphi_\mathrm{acc,ref}),\ f_\mathrm{acc,lat,min},\ 1)$, and
   ablation $abl = \beta \max(0, z_{SL} - z_\mathrm{syn})$,
   $\beta = \max(\beta_\mathrm{floor},\ \beta_0 + \beta_1(\varphi_\mathrm{ref} - \varphi))$.
4. Outside the mask the SMB is at most `smb_min`.
5. **Surface temperature** $T_\mathrm{srf} = T_{sl} - \gamma_t \max(z_\mathrm{syn}, 0)$,
   capped at `t_ice_max` on the mask.

So the SMB follows the target mask, not the geometry the ice model makes.

## Target mask

| Setting | Mask | Imports |
|---|---|---|
| default (`mask_file = ""`, `mask_ice = False`) | `H_ice_ref > 0` of the ice, every call (it follows a reconstruction start) | `H_ice_ref` |
| `mask_file`, `mask_var`, `mask_idx` | slice `mask_idx` (third dimension) of `mask_var` in `mask_file` on the surface grid, cells > 0; read once at init (`{domain}`, `{grid_name}` replaced) | none |
| `mask_ice = True` | the current ice, `H_ice > 0` of the hub, every call | `H_ice` (hub) |

`mask_file` and `mask_ice` are exclusive. `mask_ice` is for an ice backend whose
`H_ice_ref` is a fixed reference while its ice follows a reconstruction: the plastic
backend exports today's ice as `H_ice_ref` (see [PlasticIceSheet](backend-plastic.md)),
so the default mask would be the present-day ice at every time.

## Without a climate

`tsl_const > 0` [K] sets the sea-level temperature to that constant and `tsl_ann` is
not imported, so a domain needs no climate component. The SMB does not use the
temperature; only `T_srf` does. `tsl_const ≤ 0` (default -1) imports `tsl_ann`.

## Prior of the inversion

In `par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml` each of the three domains (North,
Antarctica, PatagoniaPaleo) has a surface instance `smb_simple` whose `smb` is the
prior of the plastic backend's joint inversion (`mass_balance = "inverted"`, see
[Implied mass balance](experiment-implied-mb.md)): a rough first estimate, so the
inversion gives it a large error (`inv_sigma_smb_rel = 0.5`, `inv_sigma_smb_min = 100`
kg m-2 yr-1). The three instances share one parameter group:

```fortran
&surface_nor
    backend    = "smb_simple"
    smb_simple = "smb_simple_plastic"
    grid       = ""             ! the hub grid
    dt         = 0.0            ! every step
    out_2D     = "out_2D"
/

&smb_simple_plastic
    ! the prior of the inversion: smb_simple's defaults (not tuned), the mask the current ice, no climate
    mask_ice   = True
    tsl_const  = 263.15       ! [K] only T_srf uses it (the plastic ice does not read T_srf)
/
```

The Yelmo variant of the science case, `par/gate/g6-smb.nml`, also uses `smb_simple`
with the untuned module defaults.

## Keys

Template `&smb_simple` of `input/kryos_defaults.nml` (the defaults are the yelmox
module's; a par file of a science configuration lists its tuned keys):

| Key | Default | |
|---|---|---|
| `scheme` | `"syn"` | the only one implemented (`"tg24"` is not) |
| `co2_const`, `f_const` | 280.0, 494.5 | [ppm], [W m-2] (enter the snow line only if `a3`, `a4` ≠ 0) |
| `mask_file`, `mask_var`, `mask_idx` | `""`, `""`, 1 | target mask from a file |
| `mask_ice` | `False` | the mask is the current ice |
| `tsl_const` | -1.0 | [K] > 0: constant sea-level temperature |
| `a1`, `a2`, `a3`, `a4`, `fmean` | 3750.0, -50.0, 0.0, 0.0, 494.5 | snow line [m], [m/°N], [m/ppm], [m/(W m-2)], [W m-2] |
| `beta0`, `beta1`, `phiref`, `beta_floor` | 4.0, 0.05, 60.0, 0.0 | ablation gradient [(mm/yr)/m], per °N, [°N], floor |
| `c0`, `L_acc` | 150.0, 1000.0 | accumulation at the margin [mm w.e./yr], inland e-folding length [km] |
| `k_acc_lat`, `phi_acc_ref`, `facc_lat_min` | 0.0, 0.0, 0.3 | latitude factor of the accumulation (0 = off) |
| `use_plastic`, `slope`, `tau0` | `True`, 6.09e-3, 36.0e3 | profile inside: plastic or linear; slope [m/m]; yield stress [Pa] |
| `slope_out`, `z_max_in`, `z_max_out` | 6.09e-3, 2500.0, 200.0 | ramp outside [m/m]; inside cap, outside floor [m] |
| `gamma_t`, `t_ice_max` | 6.5e-3, 273.15 | lapse rate [K/m]; cap on the mask [K] |
| `smb_min` | -2000.0 | [mm w.e./yr] outside the mask the SMB is at most this (≤ 0) |
