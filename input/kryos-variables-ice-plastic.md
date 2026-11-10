# ice, backend plastic

The `plastic` stream of an instance with backend `plastic` (`<instance>_plastic.nc`): the state of the last
solve and its inputs, written by the backend itself (`kryos_ice_plastic`) with the precision of the `prec`
column (lossless under `[kryos] write_lossless`; a variable the table does not list is written lossless:
`converged`, `n_outer_used`, `stress_res_*`, and with mass_balance = "implied" or "inverted" `abl_margin`, `discharge`
[Gt yr-1], with "inverted" `inv_iters`, `inv_J_law`, `inv_J_smb`, `inv_J_cold`, `inv_n_bound`). See kryos-variables-ice.md for the columns. `units` must be those of the data.
The exports of the instance (`<instance>.nc`, `<instance>_ts.nc`) are the table `ice`.

| id | variable    | dimensions | units  | long_name                                                       | prec |
|----|-------------|------------|--------|-----------------------------------------------------------------|------|
|  1 | H_ice       | xc, yc     | m      | Ice thickness                                                   | 1e-1 |
|  2 | z_srf       | xc, yc     | m      | Surface elevation                                               | 1e-1 |
|  3 | z_bed       | xc, yc     | m      | Bedrock elevation (as used by the solve)                        | 1e-1 |
|  4 | z_sl        | xc, yc     | m      | Sea level (as used by the solve)                                | 1e-1 |
|  5 | mask        | xc, yc     | 1      | Ice mask of the margin chronology                               | 1e-4 |
|  6 | sdist       | xc, yc     | m      | Signed distance to the margin (positive inside the ice)         | 1e-1 |
|  7 | tau_b       | xc, yc     | Pa     | Basal shear stress                                              | 1e-1 |
|  8 | taub_domain | xc, yc     | 1      | Basal shear stress domain (0 = none)                            | 1e-4 |
|  9 | div_q       | xc, yc     | m yr-1 | Divergence of the ice flux (ice equivalent)                     | 1e-4 |
| 10 | q_land      | xc, yc     | m yr-1 | Outflow through the margin faces to land (ice equivalent)       | 1e-4 |
| 11 | q_marine    | xc, yc     | m yr-1 | Outflow through the margin faces to the ocean (ice equivalent)  | 1e-4 |
| 12 | u_inv       | xc, yc     | m yr-1 | Depth-averaged speed of the inversion                           | 1e-2 |
| 13 | u_law       | xc, yc     | m yr-1 | Depth-averaged speed of the flow law (the inversion's prior)    | 1e-2 |
