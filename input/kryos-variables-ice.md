# ice

Exports of the ice component (contracts-design section 2.1), written by the
default writer: dimensions without `time` (empty: 0-D, written to `<instance>_ts.nc`),
`units` as declared by the component, `prec` as variable_io (absolute, `r` relative,
`-` lossless). The `prec` column is the only source of rounding; 0-D rows are
lossless (time series stay lossless) and a unit is required in every row.

| id | variable       | dimensions | units       | long_name                                        | prec       |
|----|----------------|------------|-------------|--------------------------------------------------|------------|
|  1 | H_ice          | xc, yc     | m           | Ice thickness                                    | 1e-3       |
|  2 | z_srf          | xc, yc     | m           | Surface elevation                                | 1e-3       |
|  3 | z_bed          | xc, yc     | m           | Bedrock elevation (as used by the ice)           | 1e-3       |
|  4 | z_sl           | xc, yc     | m           | Sea level elevation (as used by the ice)         | 1e-3       |
|  5 | f_ice          | xc, yc     | 1           | Ice-covered fraction                             | 1e-4       |
|  6 | f_grnd         | xc, yc     | 1           | Grounded fraction                                | 1e-4       |
|  7 | dHidt          | xc, yc     | m yr-1      | Ice thickness rate of change                     | 1e-4       |
|  8 | smb            | xc, yc     | kg m-2 yr-1 | Surface mass balance                             | 1e-1       |
|  9 | bmb            | xc, yc     | kg m-2 yr-1 | Basal mass balance                               | 1e-1       |
| 10 | cmb            | xc, yc     | kg m-2 yr-1 | Calving mass balance                             | 1e-1       |
| 11 | mb_net         | xc, yc     | kg m-2 yr-1 | Net mass balance applied                         | 1e-1       |
| 12 | M_ice          |            | Gt          | Ice mass                                         | -          |
| 13 | ux_s           | xc, yc     | m yr-1      | Surface x-velocity (cell centres)                | 1e-3,1e-4r |
| 14 | uy_s           | xc, yc     | m yr-1      | Surface y-velocity (cell centres)                | 1e-3,1e-4r |
| 15 | dz_bed_dt_corr | xc, yc     | m yr-1      | Bedrock rate correction                          | 1e-4       |
| 16 | z_bed_ref      | xc, yc     | m           | Reference bedrock elevation                      | 1e-3       |
| 17 | H_ice_ref      | xc, yc     | m           | Reference ice thickness                          | 1e-3       |
| 18 | footprint      | xc, yc     | 1           | Where the ice model holds ice (static)           | 1e-4       |
| 19 | H_ice_pd       | xc, yc     | m           | Present-day ice thickness (static)               | 1e-3       |
| 20 | H_grnd         | xc, yc     | m           | Ice thickness above flotation                    | 1e-3       |
| 21 | H_grnd_pd      | xc, yc     | m           | Present-day thickness above flotation (static)   | 1e-3       |
