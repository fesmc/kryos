# hub

The variables of the hub output (`<domain>_hub.nc`, yelmox `htopo.nc`) that kryos
writes itself (`domain_hub_write`) with the precision of the `prec` column. The
five fields of fesm-utils `htopo_write_step` (`z_bed`, `H_ice`, `z_srf`, `f_grnd`,
`z_sl`) are written by that library with its own precisions and are not listed; the
integer `mask` is lossless. A variable the table does not list is written lossless.
See kryos-variables-ice.md for the columns. `units` must be those of the data.

| id | variable | dimensions | units | long_name                                                | prec |
|----|----------|------------|-------|----------------------------------------------------------|------|
|  1 | f_ice    | xc, yc     | 1     | Ice-covered fraction                                     | 1e-4 |
