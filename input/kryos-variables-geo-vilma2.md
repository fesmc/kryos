# geo, backend vilma2

The `gia` stream of the global instance (`<instance>_gia.nc`): the global ice field as VILMA2 got it, on its
Gauss grid (dimensions `lon`, `lat`, `time`), written by the backend itself (`kryos_geo_vilma2`) with the
precision of the `prec` column (lossless under `[kryos] write_lossless`; a variable the table does not list is
written lossless). The weighted fields are per unit Gauss-cell area (phase4-plan 2.3). See
kryos-variables-ice.md for the columns. `units` must be those of the data. The exports of the instance
(`<instance>.nc`, `<instance>_ts.nc`) are the table `geo`; the `vilma` stream is VILMA2's own file, with the
table `vilma-variables.md` of VILMA2.

| id | variable    | dimensions | units | long_name                                                       | prec |
|----|-------------|------------|-------|-----------------------------------------------------------------|------|
|  1 | H_ice_ext   | lon, lat   | m     | External ice, mean over the part outside the host region        | 1e-3 |
|  2 | host_w      | lon, lat   | 1     | Coverage of the host region (sum over the domains)              | 1e-4 |
|  3 | host_Ig     | lon, lat   | m     | Grounded ice of the host region, weighted                       | 1e-3 |
|  4 | host_C      | lon, lat   | 1     | Ocean function of the host region, weighted                     | 1e-4 |
|  5 | host_F      | lon, lat   | m     | Flotation term of the host region, weighted                     | 1e-3 |
