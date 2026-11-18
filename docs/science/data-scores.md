---
title: "Data and scores"
---

What every simulation is scored against, and how. All data come from
[FesmData](https://github.com/fesmc/FesmData) products; the scores from `plastic_ensemble.jl scores`, the same code and
error model for every simulation, published reconstructions included.

## Data

| group | data | n | judges |
|---|---|---|---|
| Laurentide | GAPSLIP 2.1 RSL index points | 892 | yes |
| Greenland | GAPSLIP 2.1 | 276 | yes |
| Fennoscandia | GAPSLIP 2.1 (mainland and the British Isles) | 718 | yes |
| Barents | GAPSLIP 2.1 (Svalbard, Franz Josef Land, Novaya Zemlya, the Kara Sea) | 370 | yes |
| Antarctica | GAPSLIP 2.1 | 47 | yes |
| far field | GAPSLIP 2.1 | 2248 | yes |
| ESL | the Spratt and Lisiecki (2016) sea-level stack, 25-10 ka every 2.5 kyr, against the run's barystatic sea level ([below](#esl)) | 7 | yes |
| GNSS Fennoscandia | present-day uplift rates, BIFROST2015 (Kierulf et al. 2021) | 164 | yes |
| GNSS Hudson Bay | Schumacher et al. (2018), elastic-corrected, and MIDAS stations (Blewitt et al. 2016) | 61 | yes |
| GNSS Greenland | GNET (Berg et al. 2024), less its elastic part ([below](#gnss-greenland-and-antarctica)) | 55 | yes |
| GNSS Antarctica | Martín-Español et al. (2016), elastic-corrected | 67 | yes |
| ice cores (Greenland) | Holocene surface elevation at Camp Century, NGRIP, GRIP and DYE-3 (Vinther et al. 2009), 11-3 ka every 2 kyr ([below](#ice-cores)) | 20 | yes |
| rotation | GIA J̇2 and true polar wander rate and direction, from the literature ([below](#rotation)) | 3 | yes |

The RSL points are those of age ≤ 25 ka: GLAC-1D's own RSL ends at 26 ka, so the published reconstructions can be scored
on the same points. GNSS rates are compared with the run's rate of bed change at its last record (elastic-corrected
data; VILMA2's centre-of-mass frame against ITRF/IGS, a few tenths of a mm/yr apart, not corrected). Until 2026-11-18
GNSS Greenland and Antarctica were reported, not judged.

### GNSS Greenland and Antarctica

Every GNSS group compares GIA alone: our model has no present-day ice loss, so its elastic response must be out of the
data. GNET (FesmData `GPS-GNET`, `rate = "gia"`) is the observed rate less Berg et al.'s elastic response to the ice
loss of 2000-2020 (the ice sheet, its peripheral glaciers and Arctic Canada), 1.5-13.8 mm/yr at the stations against
observed rates of 1.3-17.5 mm/yr; the GIA estimate left is -1.6 to +5.6 mm/yr. Its sd combines the two formal sds
(0.1-0.8 mm/yr). Martín-Español et al. (2016) removed the elastic response to the Antarctic ice loss of 2009-2013
themselves (REAR on the RATES mass trends: up to 16-19 mm/yr in the Amundsen Sea, -2.4 to -0.5 mm/yr in Queen Maud
Land; the Earth model alone 10 % of it), and atmospheric loading.

What remains and no model of ours has: in Greenland the viscous response to the loss since the Little Ice Age on a weak
upper mantle, 3-5 mm/yr on the west and south-east coasts (Adhikari et al. 2021, whose elastic plus GIA models miss the
rates by 3.4 ± 1.9 mm/yr), and the elastic part of the peripheral glaciers (up to a third of the elastic rebound); in
Antarctica the response to the recent loss of the Amundsen Sea and the northern Peninsula on a weak mantle (BERP 19.7,
FONP 13.1, DUPT 10.0 mm/yr; every GIA model of Martín-Español et al. misses the Amundsen Sea by 7.7-9.8 mm/yr WRMS). The
error model carries it: σ_m and L_x fitted as Fennoscandia's and Hudson Bay's, from the 6 best Earths of wave A of the
Earth search (`plastic_ensemble.jl errmodel refs=E1D-ES1-E`; `kryos-work/scores-v3/errfit-gnss`): Greenland 2.0 mm/yr
over 1000 km (a common offset of the coast), Antarctica 4.0 mm/yr over 250 km. Every station is kept (Martín-Español et
al. themselves left out only TOMO, not in the product).

### Ice cores

The data: FesmData `IceCores-Vinther2009`, the surface elevation change relative to the present at Camp Century, NGRIP,
GRIP and DYE-3 (Vinther et al. 2009), from the δ18O of the cores against that of the small ice caps of Agassiz and Renland,
corrected for upstream flow and for the isostatic uplift of the two ice caps (from the relative sea level near them): the
elevation above the sea level of the time, the bed's isostatic change included. 1 sigma 40 m (GRIP, NGRIP), 65 m (Camp
Century, DYE-3). The points: at 11, 9, 7, 5 and 3 ka the mean of the 20-yr values within ±0.25 kyr, every 2 kyr since the
data errors (the uplift correction above all) are correlated in time and the covariance has them as independent. Before
8 ka the correction of Agassiz is uncertain (the Innuitian ice; Lecavalier et al. 2013 revised the curves downward, with
a larger uncertainty at GRIP and NGRIP): the revision is not used here.

The model, alike: the surface above the sea level of the time, `z_srf - z_sl` of the run's `ice_nor.nc` (the plastic ice
on NH-32KM, `z_sl` the sea level the ice sees), bilinear at the site in the grid's projection, linear in time between the
records, relative to 0 ka; a standalone VILMA2 run under a reconstruction's ice `z_bed + h_ice` of its stream; a
published reconstruction its `z_topo` (`pm_compare.jl` `core_at`, `core_product`). Error model: σ_m 75 m over 500 km and 5
kyr, starting values (`[likelihood."ice cores (Greenland)"]`): the plastic surface on 32 km misses the present elevation
of the sites by 30-380 m (RMS 240 m, C7), most of which cancels in the change relative to the present; to be refitted by
`errmodel` from an ensemble.

### Rotation

The data (`[rotation]` of the config; value, sd and the model error sd_model per quantity, the point's sd their
combination; σ_m 0):

| quantity | value | sd | sd_model | source |
|---|---|---|---|---|
| J̇2 [10⁻¹¹/yr] (unnormalized; J2 = -√5 C̄20) | -5.7 | 0.6 | 1.0 | observed -3.7 ± 0.1 (SLR 1976-1990, -3.4 to -4.0 by the window) less the melt of 1900-1990, +2.0 ± 0.4 (Nakada et al. 2015, Table 3; Mitrovica et al. 2015) |
| TPW rate [°/Myr] | 0.6 | 0.3 | 0.2 | observed 0.97 °/Myr (3.5 mas/yr) toward 79°W relative to the mean lithosphere (Argus and Gross 2004); the GIA part 33 ± 18 % (Adhikari et al. 2018) to 0.85 °/Myr (Nakada et al. 2015; Adhikari and Ivins 2016), the rest mantle convection; the 20th-century melt < 0.1 °/Myr |
| TPW direction [°E] | -77 | 15 | 15 | 75-79°W (Argus and Gross 2004; Nakada et al. 2015) |

The model: over the last interval of the run's VILMA2 records, J̇2 from the degree-2 zonal coefficient of the sea-surface
change `dz_ss` (J2 = -(5/4πa) ∫ N P2 dΩ on the Gauss grid; the uniform shift and VILMA2's rotational potential, order 1,
do not enter), and the polar motion from `rot_m_re`, `rot_m_im` (m1 toward Greenwich, m2 toward 90°E; at 20 ka every run
has its pole 2·10⁻³ rad toward 102-106°E, away from the Laurentide load, as it should) (`pm_compare.jl`
`rotation_rates`). The interval is 100 yr in a standalone run and 500 yr in a kryos run, and the run's own ice change
over it is in the rates (a GIA model's late-Holocene ice is part of its GIA; the 20th-century melt is in the data's
non-GIA part). sd_model: the loading before the run's start at 80 ka, the record interval (C7's J̇2 over the last two 500-yr
intervals: -1.5, -0.5). A kryos run's VILMA2 stream does not write the polar motion (`kryos_geo_vilma2.f90` `VARS_VILMA`
lacks `rot_m_re`, `rot_m_im`): its rotation is J̇2 alone (n = 1). The direction is that of the rate: where the rate is
near zero it is noise.

### ESL

The data: `sl_short` of the stack (FesmData `SeaLevel-Spratt2016`, 1 sigma `sl_short_sd`) relative to its own present
(+8.5 m at 0 ka), linear in time between its 1-kyr records at 25, 22.5, ..., 10 ka (22.5 ka is the mean of 22 and 23 ka).
The model: VILMA2's barystatic sea level `bsl` of the run's global instance (`geo_glob_vilma.nc`; `vilma_out.nc` of a
standalone run) relative to its last record (0 ka), linear between its records. `bsl` is the ocean mean of C·rsl,
∫C rsl dA / ∫C dA over the ocean of the time (VILMA2 `update_bsl`); by the mass balance of the sea-level equation it is the
melt water of the grounded-ice change plus the fill of newly flooded area. The ocean that fills the cavity of lost marine
ice does not raise it, so in effect the ice above flotation counts, with the change of the ocean area: the
ice-volume-equivalent sea level the stack stands for. The error model of the group has only a time term (σ_m 5 m, L_t 10
kyr). `[esl] model_var` of the scoring config selects the model's series (`tools/levante/plastic_ensemble.jl`).

Until 2026-11-18 the model's series was `bsl_ive`, −(ρ_i/ρ_w) of the whole grounded-column change over the ocean area, the
ice below flotation included (about 10 m more ice at the LGM than `bsl`), and the data were taken at their nearest record
(22, 17, 12 ka for 22.5, 17.5, 12.5 ka). `model_var = "bsl_ive"` gives that comparator back (with the data linear in
time). With `bsl` every simulation reads 10-11 m less ice at 20 ka and its ESL I rises by about 4-5 (where the run is
short of the stack); the interpolation alone moves I by up to ±0.8.

## The error model

Per data group g, the residuals d = model − observation have the covariance

R = diag(σ_i²) + σ_m,g² C,  C_ij = exp(−dist_ij / L_x,g) · exp(−|Δt_ij| / L_t,g)

with σ_i the data errors and a model error σ_m,g correlated over distance L_x,g and time L_t,g. The model error covers
what no simulation can fit (the resolution, the ice outside the reconstruction, the solver). Its parameters are fitted
from the residuals of the best members of an ensemble (ens-v2 for the RSL groups; wave A of the Earth search for GNSS:
Fennoscandia 0.7 mm/yr over 750 km, Hudson Bay 1.5 mm/yr over 500 km, Greenland 2.0 mm/yr over 1000 km, Antarctica 4.0
mm/yr over 250 km) and then held fixed for every simulation. The ice cores have starting values ([above](#ice-cores)), the
rotation group a model error per quantity.

## The scores

- **log L**, the Gaussian log-likelihood −½ d′R⁻¹d − ½ log|R| summed over the six RSL groups: higher is better.
  Differences of a few tens are large. This is the one scale of all the tables.
- **I**, the implausibility of a group, (χ² − n)/√(2n) with χ² = d′R⁻¹d: about 0 when the residuals are as large as the
  error model expects, ≤ 3 passes.
- **passed**, the number of judging groups within the bar: of 13 (six RSL groups, ESL, four GNSS groups, ice cores,
  rotation; of 9 until 2026-11-18), over the groups a series has a score in: a published reconstruction has no GNSS and
  no rotation, a run before 2026-11-12 no GNSS. In the Earth search ESL is reported, not judged.
- **RMSE** [m] and **RMS z** (the root mean square of the residuals normalized by σ_i, the score of the first
  comparisons) are reported for reference only. A few points with tiny σ_i dominate RMS z, Greenland above all.

## The bar

A group passes at I ≤ 3. Where both published reconstructions fail a group, the cut is raised to the better of their
two scores: Greenland 7.9 (PaleoMIST). Every other RSL group keeps 3, and so do the ice cores (GLAC-1D 2.6). No
publication has a GNSS or rotation score: those groups keep 3.

ESL's target is the Spratt and Lisiecki (2016) stack, and its bar is the plain cut 3. No published reconstruction sets
it: PaleoMIST no longer sets the ESL bar (until 2026-11-18 the Earth search and the Coulomb tables took PaleoMIST's own
score, 7.04, whose ESL is its whole ice volume, below flotation included). Where it is judged:

- **Our ice** (the Coulomb runs, the ensembles, the joint ice and Earth ensemble): judged at 3. ESL is an independent
  constraint on the ice, which the publications were not built against.
- **An Earth under imposed ice** (the Earth search): reported, not judged. Under imposed ice ESL is a property of the
  ice, not of the Earth (10.0-11.0 under PaleoMIST's ice on every Earth, 7.3-8.4 under GLAC-1D's), so it would rule
  Earths out for the ice's sake.

The code holds to this: the ESL group always keeps the plain cut (`plastic_ensemble.jl` `plain_cut_groups`), and the
Earth search's config has `[esl] judge = false` (`experiments/earth-search.sh`).

A simulation is *as good as the publications* when it passes every group, and *better* when it does so with a higher
log L. The Earth search uses the same bar to rule Earths out, ESL aside ([The Earth search](earth-search.md)).
