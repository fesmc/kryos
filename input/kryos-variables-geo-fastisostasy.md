# geo, backend fastisostasy

The 2D stream of the regional FastIsostasy instance (`<instance>.nc`, yelmox
`isos.nc`), written by the backend itself (`kryos_geo_fastisostasy`) with the
precision of the `prec` column. A backend variable the table does not list is
written lossless; the time series (`<instance>_ts.nc`) are lossless. See
kryos-variables-ice.md for the columns. `units` must be those of the data.

| id | variable | dimensions | units | long_name              | prec |
|----|----------|------------|-------|------------------------|------|
|  1 | z_bed    | xc, yc     | m     | Bedrock elevation      | 1e-3 |
|  2 | z_ss     | xc, yc     | m     | Sea-surface height     | 1e-3 |
|  3 | w        | xc, yc     | m     | Viscous displacement   | 1e-3 |
|  4 | we       | xc, yc     | m     | Elastic displacement   | 1e-3 |
