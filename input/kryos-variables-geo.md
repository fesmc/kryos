# geo

Exports of the geo component, regional and global part (contracts-design
section 2.2). See kryos-variables-ice.md for the columns. The 0-D rows are
lossless (`-`): time series stay lossless.

| id | variable  | dimensions | units  | long_name                                          | prec  |
|----|-----------|------------|--------|----------------------------------------------------|-------|
|  1 | dz_bed    | xc, yc     | m      | Bedrock displacement (viscous + elastic)           | 1e-3  |
|  2 | dz_ss     | xc, yc     | m      | Sea-surface displacement relative to sea level     | 1e-3  |
|  3 | z_bed_ref | xc, yc     | m      | Reference bedrock elevation (consistency check)    | 1e-3  |
|  4 | dz_bed_dt | xc, yc     | m yr-1 | Bedrock displacement rate                          | 1e-4  |
|  5 | bsl       |            | m      | Barystatic sea level                               | -     |
|  6 | A_ocn     |            | m2     | Ocean area                                         | -     |
|  7 | dV_bsl    |            | m3     | Volume change contributing to barystatic sea level | -     |
|  8 | dphi      |            | m      | Uniform sea-surface shift (dz_ss is relative to it) | -     |
|  9 | C_ocn     | xc, yc     | 1      | Ocean function (global part, host-resolved Earth)  | 1e-4  |
| 10 | bsl_ive   |            | m      | Ice-volume equivalent sea level (global part)      | -     |
| 11 | mass_resid |           | 1      | Worst sea-level-equation mass residual of the step | -     |
| 12 | bsl_ive_host |         | m      | bsl_ive's host part: the host columns' grounded-ice change (Yelmo's flotation) | -     |
