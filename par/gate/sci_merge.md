# Merged par file `sci.nml`

Sources (domain, instance suffix):
- `par/gate/nor.nml`: North (`nor`)
- `par/gate/ant.nml`: Antarctica (`ant`)

Run-wide groups (taken once, the files agree except where resolved): `&kryos`, `&clock`, `&out_1D`, `&out_2D`, `&out_2Dsm`, `&out_rst`, `&globe`, `&forcing`, `&geo_glob`, `&barysealevel`, `&coef_alpha`, `&coef_dTo`.

## Resolved keys

| key | values of the files | taken |
|---|---|---|
| `clock.time_const` | `-21e3`, `0.0` | `0.0` |
| `barysealevel.sl_path` | `"input/sl_deglaciation.dat"`, `"none"` | `"none"` |

## Renamed groups (46)

- North: `&topo` -> `&topo_nor`, `&regions` -> `&regions_nor`, `&isos` -> `&isos_nor`, `&snapshots` -> `&snapshots_nor`, `&orbit` -> `&orbit_nor`, `&surface_chion` -> `&surface_chion_nor`, `&chion` -> `&chion_nor`, `&chion_const` -> `&chion_const_nor`, `&pdd_chion` -> `&pdd_chion_nor`, `&marine_shelf` -> `&marine_shelf_nor`, `&yelmo` -> `&yelmo_nor`, `&ytopo` -> `&ytopo_nor`, `&ycalv` -> `&ycalv_nor`, `&ydyn` -> `&ydyn_nor`, `&ytill` -> `&ytill_nor`, `&yhyd` -> `&yhyd_nor`, `&ymat` -> `&ymat_nor`, `&ytherm` -> `&ytherm_nor`, `&yelmo_masks` -> `&yelmo_masks_nor`, `&yelmo_init_topo` -> `&yelmo_init_topo_nor`, `&yelmo_data` -> `&yelmo_data_nor`, `&ytrc` -> `&ytrc_nor`, `&snapshot_pd` -> `&snapshot_pd_nor`
- Antarctica: `&topo` -> `&topo_ant`, `&regions` -> `&regions_ant`, `&isos` -> `&isos_ant`, `&snapshots` -> `&snapshots_ant`, `&orbit` -> `&orbit_ant`, `&surface_chion` -> `&surface_chion_ant`, `&chion` -> `&chion_ant`, `&chion_const` -> `&chion_const_ant`, `&pdd_chion` -> `&pdd_chion_ant`, `&marine_shelf` -> `&marine_shelf_ant`, `&yelmo` -> `&yelmo_ant`, `&ytopo` -> `&ytopo_ant`, `&ycalv` -> `&ycalv_ant`, `&ydyn` -> `&ydyn_ant`, `&ytill` -> `&ytill_ant`, `&yhyd` -> `&yhyd_ant`, `&ymat` -> `&ymat_ant`, `&ytherm` -> `&ytherm_ant`, `&yelmo_masks` -> `&yelmo_masks_ant`, `&yelmo_init_topo` -> `&yelmo_init_topo_ant`, `&yelmo_data` -> `&yelmo_data_ant`, `&ytrc` -> `&ytrc_ant`, `&snapshot_pd` -> `&snapshot_pd_ant`

## Repointed keys (46)

| group | key | was | now |
|---|---|---|---|
| `&North` | `topo` | `topo` | `topo_nor` |
| `&North` | `regions` | `regions` | `regions_nor` |
| `&ice_nor` | `yelmo` | `yelmo` | `yelmo_nor` |
| `&geo_nor` | `isos` | `isos` | `isos_nor` |
| `&climate_nor` | `snapshots` | `snapshots` | `snapshots_nor` |
| `&climate_nor` | `reference` | `pd` | `pd_nor` |
| `&surface_nor` | `surface_chion` | `surface_chion` | `surface_chion_nor` |
| `&surface_nor` | `chion` | `chion` | `chion_nor` |
| `&surface_nor` | `orbit` | `orbit` | `orbit_nor` |
| `&shelf_nor` | `marine_shelf` | `marine_shelf` | `marine_shelf_nor` |
| `&chion_nor` | `nml_pdd` | `pdd_chion` | `pdd_chion_nor` |
| `&chion_nor` | `nml_const` | `chion_const` | `chion_const_nor` |
| `&yelmo_nor` | `nml_ytopo` | `ytopo` | `ytopo_nor` |
| `&yelmo_nor` | `nml_ycalv` | `ycalv` | `ycalv_nor` |
| `&yelmo_nor` | `nml_ydyn` | `ydyn` | `ydyn_nor` |
| `&yelmo_nor` | `nml_ytill` | `ytill` | `ytill_nor` |
| `&yelmo_nor` | `nml_ymat` | `ymat` | `ymat_nor` |
| `&yelmo_nor` | `nml_ytrc` | `ytrc` | `ytrc_nor` |
| `&yelmo_nor` | `nml_ytherm` | `ytherm` | `ytherm_nor` |
| `&yelmo_nor` | `nml_yhyd` | `yhyd` | `yhyd_nor` |
| `&yelmo_nor` | `nml_masks` | `yelmo_masks` | `yelmo_masks_nor` |
| `&yelmo_nor` | `nml_init_topo` | `yelmo_init_topo` | `yelmo_init_topo_nor` |
| `&yelmo_nor` | `nml_data` | `yelmo_data` | `yelmo_data_nor` |
| `&Antarctica` | `topo` | `topo` | `topo_ant` |
| `&Antarctica` | `regions` | `regions` | `regions_ant` |
| `&ice_ant` | `yelmo` | `yelmo` | `yelmo_ant` |
| `&geo_ant` | `isos` | `isos` | `isos_ant` |
| `&climate_ant` | `snapshots` | `snapshots` | `snapshots_ant` |
| `&climate_ant` | `reference` | `pd` | `pd_ant` |
| `&surface_ant` | `surface_chion` | `surface_chion` | `surface_chion_ant` |
| `&surface_ant` | `chion` | `chion` | `chion_ant` |
| `&surface_ant` | `orbit` | `orbit` | `orbit_ant` |
| `&shelf_ant` | `marine_shelf` | `marine_shelf` | `marine_shelf_ant` |
| `&chion_ant` | `nml_pdd` | `pdd_chion` | `pdd_chion_ant` |
| `&chion_ant` | `nml_const` | `chion_const` | `chion_const_ant` |
| `&yelmo_ant` | `nml_ytopo` | `ytopo` | `ytopo_ant` |
| `&yelmo_ant` | `nml_ycalv` | `ycalv` | `ycalv_ant` |
| `&yelmo_ant` | `nml_ydyn` | `ydyn` | `ydyn_ant` |
| `&yelmo_ant` | `nml_ytill` | `ytill` | `ytill_ant` |
| `&yelmo_ant` | `nml_ymat` | `ymat` | `ymat_ant` |
| `&yelmo_ant` | `nml_ytrc` | `ytrc` | `ytrc_ant` |
| `&yelmo_ant` | `nml_ytherm` | `ytherm` | `ytherm_ant` |
| `&yelmo_ant` | `nml_yhyd` | `yhyd` | `yhyd_ant` |
| `&yelmo_ant` | `nml_masks` | `yelmo_masks` | `yelmo_masks_ant` |
| `&yelmo_ant` | `nml_init_topo` | `yelmo_init_topo` | `yelmo_init_topo_ant` |
| `&yelmo_ant` | `nml_data` | `yelmo_data` | `yelmo_data_ant` |
