# climate, backend snapshots

The 2D stream of the climate instance with the snapshots backend (`<instance>.nc`,
yelmox `snap.nc`), written by the backend itself (`kryos_climate_snapshots`, the
snapshots library and the climate transforms) with the precision of the `prec`
column. A backend variable the table does not list is written lossless (the
static fields of the transforms, `taper_w`, `bias_dtas_ann`, `bias_fpr_ann`, and
the term coefficients `term_<name>`). See kryos-variables-ice.md for the columns.
`units` must be those of the data.

| id | variable          | dimensions     | units | long_name                                                           | prec  |
|----|-------------------|----------------|-------|---------------------------------------------------------------------|-------|
|  1 | t2m_ann           | xc, yc         | K     | Near-surface air temperature (ann)                                  | 1e-3  |
|  2 | t2m_sum           | xc, yc         | K     | Near-surface air temperature (sum)                                  | 1e-3  |
|  3 | pr_ann            | xc, yc         | mm/a  | Precipitation (ann)                                                 | 1e-4r |
|  4 | to_ann            | xc, yc, depth  | K     | Ocean temperature (ann)                                             | 1e-3  |
|  5 | so_ann            | xc, yc, depth  | PSU   | Ocean salinity (ann)                                                | 1e-4r |
|  6 | tas               | xc, yc, month  | K     | Near-surface air temperature                                        | 1e-3  |
|  7 | pr                | xc, yc, month  | mm/d  | Precipitation                                                       | 1e-4r |
|  8 | tsl_ann           | xc, yc         | K     | Sea-level air temperature (ann)                                     | 1e-3  |
|  9 | z_srf             | xc, yc         | m     | Surface elevation                                                   | 1e-3  |
| 10 | bias_dtas_ann_now | xc, yc         | K     | Bias correction of the air temperature at the forcing time (ann)    | 1e-3  |
