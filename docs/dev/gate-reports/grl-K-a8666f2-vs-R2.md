# kryos-compare: PASS (exit 0)

- ref: `/private/tmp/claude-1319099207/-Users-alrobi001-models-kryos/8261b972-1ee8-4f0c-82dc-aa36928f8b66/scratchpad/gate1/R2`
- run: `/private/tmp/claude-1319099207/-Users-alrobi001-models-kryos/8261b972-1ee8-4f0c-82dc-aa36928f8b66/scratchpad/gate1/K`
- map: `tools/gate/map_phase1.toml` (mapped), gate `grl`
- date: 2026-10-25 11:04
- compared: 202 variables, 0 failing, 0 absent (A), 12 new (N)

## First divergence

None.

## Variables

| ref | run | class | status | cells | differing | max abs | max rel | max ulp | first diff (t, rec) | loc max | note |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `chion.nc:H_snow` | `surface_grl.nc:H_snow` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:PDDs` | `surface_grl.nc:PDDs` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:S` | `surface_grl.nc:S` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:alb_s` | `surface_grl.nc:alb_s` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:mask` | `surface_grl.nc:mask` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:melt` | `surface_grl.nc:melt` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:pr` | `surface_grl.nc:pr` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:refrz` | `surface_grl.nc:refrz` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:runoff` | `surface_grl.nc:runoff` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:sf` | `surface_grl.nc:sf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:smb` | `surface_grl.nc:smb` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:t2m` | `surface_grl.nc:t2m` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:tsrf` | `surface_grl.nc:tsrf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:xc` | `surface_grl.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `chion.nc:yc` | `surface_grl.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:H_ice` | `Greenland_hub.nc:H_ice` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:NEGIS` | `Greenland_hub.nc:NEGIS` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:NEGIS.mask` | `Greenland_hub.nc:NEGIS.mask` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:Zwally2012` | `Greenland_hub.nc:Zwally2012` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:Zwally2012.group` | `Greenland_hub.nc:Zwally2012.group` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:Zwally2012.mask` | `Greenland_hub.nc:Zwally2012.mask` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:basin` | `Greenland_hub.nc:basin` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:f_grnd` | `Greenland_hub.nc:f_grnd` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:region` | `Greenland_hub.nc:region` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:xc` | `Greenland_hub.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:yc` | `Greenland_hub.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:z_bed` | `Greenland_hub.nc:z_bed` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:z_sl` | `Greenland_hub.nc:z_sl` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:z_srf` | `Greenland_hub.nc:z_srf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `htopo.nc:zone` | `Greenland_hub.nc:zone` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:w` | `geo_grl_isos.nc:w` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:w,we` | `geo_grl.nc:dz_bed` | C | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:we` | `geo_grl_isos.nc:we` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:xc` | `geo_grl_isos.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:xc` | `geo_grl.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:yc` | `geo_grl_isos.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:yc` | `geo_grl.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:z_bed` | `geo_grl_isos.nc:z_bed` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `isos.nc:z_ss` | `geo_grl_isos.nc:z_ss` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:bsl` | `geo_grl_isos_ts.nc:bsl` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:w_mean` | `geo_grl_isos_ts.nc:w_mean` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:we_mean` | `geo_grl_isos_ts.nc:we_mean` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:z_bed_max` | `geo_grl_isos_ts.nc:z_bed_max` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:z_bed_mean` | `geo_grl_isos_ts.nc:z_bed_mean` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `isos_ts.nc:z_bed_min` | `geo_grl_isos_ts.nc:z_bed_min` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `mshlf.nc:T_shlf` | `shelf_grl.nc:T_shlf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `mshlf.nc:bmb_shlf` | `shelf_grl.nc:bmb_shlf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `mshlf.nc:tf_shlf` | `shelf_grl.nc:tf_shlf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `mshlf.nc:xc` | `shelf_grl.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `mshlf.nc:yc` | `shelf_grl.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:depth` | `climate_grl.nc:depth` | E | pass | 23 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:month` | `climate_grl.nc:month` | E | pass | 12 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:pr` | `climate_grl.nc:pr` | E | pass | 2532552 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:pr_ann` | `climate_grl.nc:pr_ann` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:so_ann` | `climate_grl.nc:so_ann` | E | pass | 4854058 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:t2m_ann` | `climate_grl.nc:t2m_ann` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:t2m_sum` | `climate_grl.nc:t2m_sum` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:tas` | `climate_grl.nc:tas` | E | pass | 2532552 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:term_lgm` | `climate_grl.nc:term_lgm` | E | pass | 11 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:term_ocean` | `climate_grl.nc:term_ocean` | E | pass | 11 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:to_ann` | `climate_grl.nc:to_ann` | E | pass | 4854058 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:tsl_ann` | `climate_grl.nc:tsl_ann` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:xc` | `climate_grl.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:yc` | `climate_grl.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `snap.nc:z_srf` | `climate_grl.nc:z_srf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:H_ice` | `ice_grl_yelmo.nc:H_ice` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:H_ice` | `ice_grl.nc:H_ice` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:H_ice_ref` | `ice_grl_yelmo.nc:H_ice_ref` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:H_sed` | `ice_grl_yelmo.nc:H_sed` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:N_eff` | `ice_grl_yelmo.nc:N_eff` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:T_prime_b` | `ice_grl_yelmo.nc:T_prime_b` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:area` | `ice_grl_yelmo.nc:area` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:basins` | `ice_grl_yelmo.nc:basins` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:beta` | `ice_grl_yelmo.nc:beta` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:bmb` | `ice_grl_yelmo.nc:bmb` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:bmb` | `ice_grl.nc:bmb` | C | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:bmb_shlf` | `ice_grl_yelmo.nc:bmb_shlf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:cb_ref` | `ice_grl_yelmo.nc:cb_ref` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:cmb` | `ice_grl_yelmo.nc:cmb` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:cmb` | `ice_grl.nc:cmb` | C | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:crs` | `ice_grl_yelmo.nc:crs` | E | pass | 1 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:hyd_W_til` | `ice_grl_yelmo.nc:hyd_W_til` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:lat2D` | `ice_grl_yelmo.nc:lat2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:lon2D` | `ice_grl_yelmo.nc:lon2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:mask_bed` | `ice_grl_yelmo.nc:mask_bed` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:mask_ice` | `ice_grl_yelmo.nc:mask_ice` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:mb_net` | `ice_grl_yelmo.nc:mb_net` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:mb_net` | `ice_grl.nc:mb_net` | C | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:month` | `ice_grl_yelmo.nc:month` | E | pass | 12 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:pc_steps` | `ice_grl_yelmo.nc:pc_steps` | E | pass | 3 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:pd_time_iso` | `ice_grl_yelmo.nc:pd_time_iso` | E | pass | 10 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:regions` | `ice_grl_yelmo.nc:regions` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:smb` | `ice_grl_yelmo.nc:smb` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:smb` | `ice_grl.nc:smb` | C | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:taub` | `ice_grl_yelmo.nc:taub` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:taud` | `ice_grl_yelmo.nc:taud` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:time_iso` | `ice_grl_yelmo.nc:time_iso` | E | pass | 1 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:ux_bar` | `ice_grl_yelmo.nc:ux_bar` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:uxy_b` | `ice_grl_yelmo.nc:uxy_b` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:uxy_bar` | `ice_grl_yelmo.nc:uxy_bar` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:uxy_s` | `ice_grl_yelmo.nc:uxy_s` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:uy_bar` | `ice_grl_yelmo.nc:uy_bar` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:visc_bar` | `ice_grl_yelmo.nc:visc_bar` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:x2D` | `ice_grl_yelmo.nc:x2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:xc` | `ice_grl_yelmo.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:xc` | `ice_grl.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:y2D` | `ice_grl_yelmo.nc:y2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:yc` | `ice_grl_yelmo.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:yc` | `ice_grl.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:z_bed` | `ice_grl_yelmo.nc:z_bed` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:z_bed_sd` | `ice_grl_yelmo.nc:z_bed_sd` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:z_sl` | `ice_grl_yelmo.nc:z_sl` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:z_srf` | `ice_grl_yelmo.nc:z_srf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:z_srf` | `ice_grl.nc:z_srf` | E | pass | 211046 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:zeta` | `ice_grl_yelmo.nc:zeta` | E | pass | 6 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:zeta_ac` | `ice_grl_yelmo.nc:zeta_ac` | E | pass | 7 | 0 | 0 | 0 |  |  |  |  |
| `yelmo.nc:zeta_rock` | `ice_grl_yelmo.nc:zeta_rock` | E | pass | 5 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:H_ice` | `ice_grl_yelmo_sm.nc:H_ice` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:H_sed` | `ice_grl_yelmo_sm.nc:H_sed` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:T_prime_b` | `ice_grl_yelmo_sm.nc:T_prime_b` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:area` | `ice_grl_yelmo_sm.nc:area` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:basins` | `ice_grl_yelmo_sm.nc:basins` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:bmb` | `ice_grl_yelmo_sm.nc:bmb` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:crs` | `ice_grl_yelmo_sm.nc:crs` | E | pass | 1 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:hyd_W_til` | `ice_grl_yelmo_sm.nc:hyd_W_til` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:lat2D` | `ice_grl_yelmo_sm.nc:lat2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:lon2D` | `ice_grl_yelmo_sm.nc:lon2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:mask_bed` | `ice_grl_yelmo_sm.nc:mask_bed` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:mask_ice` | `ice_grl_yelmo_sm.nc:mask_ice` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:month` | `ice_grl_yelmo_sm.nc:month` | E | pass | 12 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:pc_steps` | `ice_grl_yelmo_sm.nc:pc_steps` | E | pass | 3 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:pd_time_iso` | `ice_grl_yelmo_sm.nc:pd_time_iso` | E | pass | 10 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:regions` | `ice_grl_yelmo_sm.nc:regions` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:smb` | `ice_grl_yelmo_sm.nc:smb` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:time_iso` | `ice_grl_yelmo_sm.nc:time_iso` | E | pass | 1 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:uxy_s` | `ice_grl_yelmo_sm.nc:uxy_s` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:x2D` | `ice_grl_yelmo_sm.nc:x2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:xc` | `ice_grl_yelmo_sm.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:y2D` | `ice_grl_yelmo_sm.nc:y2D` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:yc` | `ice_grl_yelmo_sm.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:z_bed` | `ice_grl_yelmo_sm.nc:z_bed` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:z_bed_sd` | `ice_grl_yelmo_sm.nc:z_bed_sd` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:z_sl` | `ice_grl_yelmo_sm.nc:z_sl` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:z_srf` | `ice_grl_yelmo_sm.nc:z_srf` | E | pass | 1937786 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:zeta` | `ice_grl_yelmo_sm.nc:zeta` | E | pass | 6 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:zeta_ac` | `ice_grl_yelmo_sm.nc:zeta_ac` | E | pass | 7 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_sm.nc:zeta_rock` | `ice_grl_yelmo_sm.nc:zeta_rock` | E | pass | 5 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:A_ice` | `ice_grl_yelmo_ts.nc:A_ice` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:A_ice_f` | `ice_grl_yelmo_ts.nc:A_ice_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:A_ice_g` | `ice_grl_yelmo_ts.nc:A_ice_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:H_ice` | `ice_grl_yelmo_ts.nc:H_ice` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:H_ice_f` | `ice_grl_yelmo_ts.nc:H_ice_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:H_ice_g` | `ice_grl_yelmo_ts.nc:H_ice_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:H_ice_max` | `ice_grl_yelmo_ts.nc:H_ice_max` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:T_shlf` | `ice_grl_yelmo_ts.nc:T_shlf` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:T_srf` | `ice_grl_yelmo_ts.nc:T_srf` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:V_ice` | `ice_grl_yelmo_ts.nc:V_ice` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:V_ice_f` | `ice_grl_yelmo_ts.nc:V_ice_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:V_ice_g` | `ice_grl_yelmo_ts.nc:V_ice_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:V_sl` | `ice_grl_yelmo_ts.nc:V_sl` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:V_sle` | `ice_grl_yelmo_ts.nc:V_sle` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:W_til` | `ice_grl_yelmo_ts.nc:W_til` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:bmb` | `ice_grl_yelmo_ts.nc:bmb` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:bmb_g` | `ice_grl_yelmo_ts.nc:bmb_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:bmb_shlf` | `ice_grl_yelmo_ts.nc:bmb_shlf` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:bmb_tot` | `ice_grl_yelmo_ts.nc:bmb_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:cmb` | `ice_grl_yelmo_ts.nc:cmb` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:cmb_flt` | `ice_grl_yelmo_ts.nc:cmb_flt` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:cmb_grnd` | `ice_grl_yelmo_ts.nc:cmb_grnd` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:coef_alpha` | `forcing_ts.nc:coef_alpha` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:coef_alpha_df_dt` | `forcing_ts.nc:coef_alpha_df_dt` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:coef_dTo` | `forcing_ts.nc:coef_dTo` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:dHidt` | `ice_grl_yelmo_ts.nc:dHidt` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:dVidt` | `ice_grl_yelmo_ts.nc:dVidt` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:dmb` | `ice_grl_yelmo_ts.nc:dmb` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:dzsdt` | `ice_grl_yelmo_ts.nc:dzsdt` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:f_pmp` | `ice_grl_yelmo_ts.nc:f_pmp` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:fmb_tot` | `ice_grl_yelmo_ts.nc:fmb_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:fwf` | `ice_grl_yelmo_ts.nc:fwf` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:mask` | `ice_grl_yelmo_ts.nc:mask` | E | pass | 19186 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:mb_border_tot` | `ice_grl_yelmo_ts.nc:mb_border_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:mb_clip_tot` | `ice_grl_yelmo_ts.nc:mb_clip_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:mb_relax_tot` | `ice_grl_yelmo_ts.nc:mb_relax_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:mb_resid_tot` | `ice_grl_yelmo_ts.nc:mb_resid_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:smb` | `ice_grl_yelmo_ts.nc:smb` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:smb_tot` | `ice_grl_yelmo_ts.nc:smb_tot` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_b` | `ice_grl_yelmo_ts.nc:uxy_b` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_b_f` | `ice_grl_yelmo_ts.nc:uxy_b_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_b_g` | `ice_grl_yelmo_ts.nc:uxy_b_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_bar` | `ice_grl_yelmo_ts.nc:uxy_bar` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_bar_f` | `ice_grl_yelmo_ts.nc:uxy_bar_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_bar_g` | `ice_grl_yelmo_ts.nc:uxy_bar_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_s` | `ice_grl_yelmo_ts.nc:uxy_s` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_s_f` | `ice_grl_yelmo_ts.nc:uxy_s_f` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:uxy_s_g` | `ice_grl_yelmo_ts.nc:uxy_s_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:xc` | `ice_grl_yelmo_ts.nc:xc` | E | pass | 106 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:yc` | `ice_grl_yelmo_ts.nc:yc` | E | pass | 181 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:z_bed` | `ice_grl_yelmo_ts.nc:z_bed` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:z_sl` | `ice_grl_yelmo_ts.nc:z_sl` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:z_srf` | `ice_grl_yelmo_ts.nc:z_srf` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:z_srf_g` | `ice_grl_yelmo_ts.nc:z_srf_g` | E | pass | 101 | 0 | 0 | 0 |  |  |  |  |
| `yelmo_ts.nc:zeta` | `ice_grl_yelmo_ts.nc:zeta` | E | pass | 6 | 0 | 0 | 0 |  |  |  |  |

## New (N)

Files: `geo_glob_bsl_ts.nc`, `geo_glob_ts.nc`, `geo_grl_ts.nc`, `ice_grl_ts.nc`

Variables: `Greenland_hub.nc:f_ice`, `Greenland_hub.nc:mask`, `geo_grl.nc:dz_ss`, `geo_grl.nc:z_bed_ref`, `geo_grl.nc:dz_bed_dt`, `ice_grl.nc:z_bed`, `ice_grl.nc:z_sl`, `ice_grl.nc:f_ice`, `ice_grl.nc:f_grnd`, `ice_grl.nc:dHidt`, `ice_grl.nc:z_bed_ref`, `ice_grl.nc:H_ice_ref`

## Deviations

- used: none
- matched nothing (remove them): none
