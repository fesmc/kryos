# geo

Exports of the geo component, regional and global part (contracts-design
section 2.2). See kryos-variables-ice.md for the columns.

| id | variable  | dimensions | units  | long_name                                          | prec  |
|----|-----------|------------|--------|----------------------------------------------------|-------|
|  1 | dz_bed    | xc, yc     | m      | Bedrock displacement (viscous + elastic)           | 1e-3  |
|  2 | dz_ss     | xc, yc     | m      | Sea-surface displacement relative to sea level     | 1e-3  |
|  3 | z_bed_ref | xc, yc     | m      | Reference bedrock elevation (consistency check)    | 1e-3  |
|  4 | dz_bed_dt | xc, yc     | m yr-1 | Bedrock displacement rate                          | 1e-4  |
|  5 | bsl       |            | m      | Barystatic sea level                               | 1e-3  |
|  6 | A_ocn     |            | m2     | Ocean area                                         | 1e-6r |
|  7 | dV_bsl    |            | m3     | Volume change contributing to barystatic sea level | 1e-6r |
