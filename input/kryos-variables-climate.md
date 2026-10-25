# climate

Exports of the climate component (contracts-design section 2.3). Derived
coefficients (0-D, named by the par file) are not listed: the default writer
writes a 0-D export without a row lossless. See kryos-variables-ice.md for the columns.

| id | variable   | dimensions    | units       | long_name                                      | prec  |
|----|------------|---------------|-------------|------------------------------------------------|-------|
|  1 | tas        | xc, yc, month | K           | Near-surface air temperature (monthly)         | 1e-3  |
|  2 | pr         | xc, yc, month | kg m-2 yr-1 | Precipitation (monthly)                        | 1e-4r |
|  3 | z_srf_clim | xc, yc        | m           | Surface elevation of tas and pr                | 1e-3  |
|  4 | lapse      | xc, yc, month | K m-1       | Near-surface lapse rate (monthly)              | 1e-6  |
|  5 | beta_p     | xc, yc        | K-1         | Precipitation scaling per kelvin               | 1e-5r |
|  6 | to_ann     | xc, yc, depth | K           | Ocean temperature profile (annual)             | 1e-3  |
|  7 | so_ann     | xc, yc, depth | psu         | Ocean salinity profile (annual)                | 1e-3  |
|  8 | dto_ann    | xc, yc, depth | K           | Ocean temperature anomaly profile (annual)     | 1e-3  |
|  9 | to_shlf    | xc, yc        | K           | Ocean temperature at the shelf base            | 1e-3  |
| 10 | so_shlf    | xc, yc        | psu         | Ocean salinity at the shelf base               | 1e-3  |
| 11 | dto_shlf   | xc, yc        | K           | Ocean temperature anomaly at the shelf base    | 1e-3  |
| 12 | dso_shlf   | xc, yc        | psu         | Ocean salinity anomaly at the shelf base       | 1e-3  |
| 13 | smb        | xc, yc        | kg m-2 yr-1 | Surface mass balance                           | 1e-1  |
| 14 | T_srf      | xc, yc        | K           | Surface temperature                            | 1e-3  |
| 15 | Qd         | xc, yc        | m3 s-1      | Subglacial discharge per cell                  | 1e-4r |
