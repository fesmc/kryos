# shelf, backend marine_shelf

The 2D stream of the shelf instance with the marine_shelf backend
(`<instance>.nc`, yelmox `mshlf.nc`), written by the backend itself
(`kryos_shelf_marine_shelf`) with the precision of the `prec` column. Its
`bmb_shlf` is in m yr-1 of ice (the export of the instance table is in
kg m-2 yr-1). A backend variable the table does not list is written lossless.
See kryos-variables-ice.md for the columns. `units` must be those of the data.

| id | variable | dimensions | units | long_name                   | prec |
|----|----------|------------|-------|-----------------------------|------|
|  1 | bmb_shlf | xc, yc     | m/yr  | Shelf basal mass balance    | 1e-4 |
|  2 | T_shlf   | xc, yc     | K     | Shelf temperature           | 1e-3 |
|  3 | tf_shlf  | xc, yc     | K     | Thermal forcing             | 1e-3 |
