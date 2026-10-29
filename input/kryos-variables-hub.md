# hub

The variables of the hub output (`<domain>_hub.nc`, yelmox `htopo.nc`), written
by `domain_hub_write` with the precision of the `prec` column. The five fields of
fesm-utils `htopo_write_step` (`z_bed`, `H_ice`, `z_srf`, `f_grnd`, `z_sl`) get
their precisions from here (the library's `htopo_write_prec_t`); their long names
and units are the library's. The integer `mask` is not listed and lossless. A
variable the table does not list is written lossless.
See kryos-variables-ice.md for the columns. `units` must be those of the data.

| id | variable | dimensions | units | long_name                                                | prec |
|----|----------|------------|-------|----------------------------------------------------------|------|
|  1 | z_bed    | xc, yc     | m     | Bedrock elevation                                        | 1e-3 |
|  2 | H_ice    | xc, yc     | m     | Ice thickness                                            | 1e-3 |
|  3 | z_srf    | xc, yc     | m     | Surface elevation                                        | 1e-3 |
|  4 | f_grnd   | xc, yc     | 1     | Grounded-ice fraction                                    | 1e-4 |
|  5 | z_sl     | xc, yc     | m     | Sea-surface height                                       | 1e-3 |
|  6 | f_ice    | xc, yc     | 1     | Ice-covered fraction                                     | 1e-4 |
