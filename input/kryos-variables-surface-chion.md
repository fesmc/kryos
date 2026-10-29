# surface, backend chion

The 2D stream of the surface instance with the chion backend (`<instance>.nc`,
yelmox `chion.nc`): the annual fields of the last annual cycle, written by the
backend itself (`kryos_surface_chion`) with the precision of the `prec` column.
A backend variable the table does not list is written lossless. See
kryos-variables-ice.md for the columns. `units` must be those of the data.

| id | variable | dimensions | units      | long_name                                  | prec  |
|----|----------|------------|------------|--------------------------------------------|-------|
|  1 | smb      | xc, yc     | mm w.e./yr | Surface mass balance (sf+rf-runoff-subl)   | 1e-1  |
|  2 | tsrf     | xc, yc     | K          | Surface temperature                        | 1e-3  |
|  3 | melt     | xc, yc     | mm w.e./yr | Melt                                       | 1e-1  |
|  4 | runoff   | xc, yc     | mm w.e./yr | Runoff                                     | 1e-1  |
|  5 | refrz    | xc, yc     | mm w.e./yr | Refreezing                                 | 1e-1  |
|  6 | pr       | xc, yc     | mm w.e./yr | Precipitation                              | 1e-4r |
|  7 | sf       | xc, yc     | mm w.e./yr | Snowfall                                   | 1e-4r |
|  8 | t2m      | xc, yc     | K          | Near-surface air temperature               | 1e-3  |
|  9 | S        | xc, yc     | W m-2      | TOA insolation                             | 1e-4r |
| 10 | PDDs     | xc, yc     | K d        | Positive degree days                       | 1e-4r |
| 11 | alb_s    | xc, yc     | 1          | Surface albedo                             | 1e-4  |
| 12 | H_snow   | xc, yc     | mm w.e.    | Snowpack thickness (end of year)           | 1e-1  |
| 13 | mask     | xc, yc     | 1          | chion active (land or ice)                 | -     |
