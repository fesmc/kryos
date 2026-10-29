# globe_ice

Output of the global geometry hub (`globe_ice.nc`, phase4-plan 5.4, 9.1), on the
background grid (dimensions `lon`, `lat`, without `time`; see kryos-variables-ice.md for
the columns). `z_bed_ref`, `H_ice_ref` and `cov` are static (written once); `H_ice` is
the merged ice of the last gather, one record per output. `prec` as variable_io
(absolute, `r` relative); lossless under `[kryos] write_lossless`.

| id | variable  | dimensions | units | long_name                                                          | prec  |
|----|-----------|------------|-------|--------------------------------------------------------------------|-------|
|  1 | z_bed_ref | lon, lat   | m     | Reference bedrock elevation (background, domains spliced in)       | 1e-3  |
|  2 | H_ice_ref | lon, lat   | m     | Reference ice thickness (background, domains spliced in)           | 1e-3  |
|  3 | cov       | lon, lat   | 1     | Fraction of the cell covered by the domains (sum of the owned parts) | 1e-4  |
|  4 | H_ice     | lon, lat   | m     | Ice thickness of the last gather (domains merged into the background) | 1e-3  |
