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
| far field | GAPSLIP 2.1, below 15 ka | 1787 | yes |
| far field LGM | GAPSLIP 2.1 index points and limiting samples, 15-30 ka, and the Red Sea proxy record ([below](#far-field-lgm-and-mis-3)) | 671 | yes |
| far field MIS 3 | the same at 30-60 ka, without Kikaijima's limiting samples, and the Karimata proxy | 176 | no (reported) |
| ESL | the Spratt and Lisiecki (2016) sea-level stack as given, 25-10 ka every 2.5 kyr, against the run's ice-mass sea level `bsl_ive`, with the spread of Spratt and Lambeck et al. (2014) in the model error ([below](#esl)) | 7 | yes |
| GNSS Fennoscandia | present-day uplift rates, BIFROST2015 (Kierulf et al. 2021) | 164 | yes |
| GNSS Hudson Bay | Schumacher et al. (2018), elastic-corrected, and MIDAS stations (Blewitt et al. 2016) | 61 | yes |
| GNSS North America forebulge | Schumacher et al. (2018), the Great Lakes, the St Lawrence and the US Atlantic coast ([below](#gnss-in-new-regions)) | 661 | yes |
| GNSS British Isles | Schumacher et al. (2018) | 112 | no (reported) |
| GNSS Greenland | GNET (Berg et al. 2024), less its elastic part ([below](#gnss-greenland-and-antarctica)) | 55 | yes |
| GNSS Antarctica | Martín-Español et al. (2016), elastic-corrected | 67 | yes |
| ice cores (Greenland) | Holocene surface elevation at Camp Century, NGRIP, GRIP and DYE-3 (Vinther et al. 2009), 11-3 ka every 2 kyr ([below](#ice-cores)) | 20 | no (reported) |
| rotation | GIA J̇2 and the true polar wander velocity along and across its observed direction, from the literature ([below](#rotation)) | 3 | yes |

The RSL points are those of age ≤ 25 ka, GLAC-1D's own RSL ending at 26 ka so that the published reconstructions can
be scored on the same points, except the far field: below 15 ka in its own group, 15-30 ka in the far field LGM and 30-60
ka in the far field MIS 3, each with its own window in every scoring window (a series without the records of a window has
no score there: `pub-GLAC1D` in the LGM and MIS 3 groups, our GLAC-1D runs, with their gap at 32.5 and 37.5 ka, in MIS 3). GNSS rates are compared with the run's rate of
bed change at its last record (elastic-corrected data; VILMA2's centre-of-mass frame against ITRF/IGS, a few tenths of
a mm/yr apart, not corrected). Until 2026-11-18 GNSS Greenland and Antarctica were reported, not judged; the far field
held every age (2248 points of age ≤ 25 ka), and the new groups of scores v4 (far field LGM/MIS 3, GNSS North America
forebulge and British Isles) did not exist. Scores v5 (2026-11-18): the far field LGM/MIS 3 split at 30 ka, MIS 3
reported only, Kikaijima's limiting samples out; ESL's spread from Spratt and Lambeck alone; the Greenland ice cores and
GNSS British Isles reported only.

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
of the sites by 30-380 m (RMS 240 m, C7), most of which cancels in the change relative to the present; kept after a
refit on the joint ensemble's wave 1 that the data do not allow ([the error model](#the-error-model)).

### Rotation

The data (`[rotation]` of the config; value, sd and the model error sd_model per quantity, the point's sd their
combination; σ_m 0). The true polar wander is a velocity, scored by its two components along the observed direction
(`tpw_axis_deg`, 79°W) and across it (90° east of it): where the model's rate is near zero its direction is noise, and a
direction scored on its own (until 2026-11-18, ±15° about 77°W) failed such runs for nothing.

| quantity | value | sd | sd_model | source |
|---|---|---|---|---|
| J̇2 [10⁻¹¹/yr] (unnormalized; J2 = -√5 C̄20) | -5.7 | 0.6 | 1.0 | observed -3.7 ± 0.1 (SLR 1976-1990, -3.4 to -4.0 by the window) less the melt of 1900-1990, +2.0 ± 0.4 (Nakada et al. 2015, Table 3; Mitrovica et al. 2015) |
| TPW along 79°W [°/Myr] | 0.6 | 0.3 | 0.2 | observed 0.97 °/Myr (3.5 mas/yr) toward 79°W relative to the mean lithosphere (Argus and Gross 2004); the GIA part 33 ± 18 % (Adhikari et al. 2018) to 0.85 °/Myr (Nakada et al. 2015; Adhikari and Ivins 2016), the rest mantle convection; the 20th-century melt < 0.1 °/Myr |
| TPW across (toward 11°E) [°/Myr] | 0.0 | 0.3 | 0.2 | the GIA models point along the observed motion, 75-79°W (Nakada et al. 2015); the non-GIA part is of unknown direction, so the across component has the along one's sd |

The model: over the last interval of the run's VILMA2 records, J̇2 from the degree-2 zonal coefficient of the sea-surface
change `dz_ss` (J2 = -(5/4πa) ∫ N P2 dΩ on the Gauss grid; the uniform shift and VILMA2's rotational potential, order 1,
do not enter), and the polar motion from `rot_m_re`, `rot_m_im` (m1 toward Greenwich, m2 toward 90°E; at 20 ka every run
has its pole 2·10⁻³ rad toward 102-106°E, away from the Laurentide load, as it should), its velocity projected on the
axis and across it (`pm_compare.jl` `rotation_rates`). The interval is 100 yr in a standalone run and 500 yr in a kryos
run, and the run's own ice change over it is in the rates (a GIA model's late-Holocene ice is part of its GIA; the
20th-century melt is in the data's non-GIA part). sd_model: the loading before the run's start at 80 ka, the record
interval (C7's J̇2 over the last two 500-yr intervals: -1.5, -0.5). Since 2026-11-18 a kryos run's VILMA2 stream writes
the polar motion when the rotation solver runs (`kryos_geo_vilma2.f90`, `rot_m_re`, `rot_m_im`); the runs of earlier
builds have J̇2 alone (n = 1). In the Earth search the group is judged: J̇2 and the polar wander are the Earth's (the
lower mantle above all) under either ice (`experiments/earth-search.sh`).

### ESL

The data: `sl_short` of the stack (FesmData `SeaLevel-Spratt2016`, 1 sigma `sl_short_sd`) as given, 0 at 5 ka by its
definition (-130 m at 24 ka, the anchor of Clark et al. 2009), linear in time between its 1-kyr records at 25, 22.5, ...,
10 ka (22.5 ka is the mean of 22 and 23 ka). Every SeaLevel curve is used as its authors give it (`pm_compare.jl`
`sealevel_at`): until 2026-11-18 the stack was re-zeroed at 0 ka, where it reads +8.5 m (the noise of the stack near the
present), which moved the whole glacial curve 8.5 m down.

The model: what the record measures. The stack is built on δ18O, which sees the mass of all grounded ice, below
flotation too: its comparator is VILMA2's `bsl_ive` of the run's global instance (`geo_glob_vilma.nc`; `vilma_out.nc` of a
standalone run), −(ρ_i/ρ_w) ∫ΔI_g dA over the ocean area of the time, relative to its last record (0 ka), linear between
its records (`[esl] model_var`, default `bsl_ive`). The same holds for the mass-based curves (Lambeck et al. 2014,
Rohling, PaleoMIST's own ESL). `bsl`, the barystatic sea level (the ocean mean of C·rsl: by the mass balance of the
sea-level equation the melt water plus the fill of newly flooded area, so in effect the ice above flotation), is what a
direct curve of the mean sea surface measures, and a site record (the Red Sea) is neither: site RSL, which the far
field LGM/MIS 3 group scores at its site. PR #55 (2026-11-18) took `bsl` for the stack's quantity: wrong, the stack
measures ice mass, not the sea surface; `bsl_ive` holds about 10 m more ice at the LGM than `bsl`.

The error model: a time term only, σ_m(t) = √(σ_m² + s(t)²) correlated over L_t, with σ_m 5 m and L_t 10 kyr and s(t)
the disagreement of the two records of the ice mass, Spratt and Lisiecki (2016) as given and Lambeck et al. (2014),
|Spratt − Lambeck|/√2, the variance averaged over ±1 kyr (`tools/levante/esl_spread.jl`, computed once; its output
`tools/levante/esl_spread.csv` is the input, `[esl] spread_file`): 2.8, 4.5, 9.3, 10.0, 11.6, 7.4, 10.4 m at 25, 22.5, ...,
10 ka. Not Grant et al. (2014: the Red Sea, site RSL with its GIA, scored at its site by the far field LGM) nor
Waelbroeck et al. (2002, a predecessor of the stack's own inputs). With them (scores v4) s(t) was 10-16 m and ESL did not
tell an 83-m LGM from a 110-m one (`kryos-work/esl-records`: JE-lo, the smoke test's low corner, 3.9 against C7's 0.6); now
JE-lo 15.7, C7 4.8, ED-pNE25 −1.0, PaleoMIST's ice on E117 2.6 ([results](results.md#scores-v5)).

### Far field LGM and MIS 3

The far-field regions of GAPSLIP 2.1 (the Caribbean, East, Southeast and South Asia, Australia, the Pacific islands)
split by age (`pm_compare.jl` `GROUPS`, `GROUP_AGES`): the far field below 15 ka (the deglaciation and the Holocene, dense
coral and mangrove index points: 1787 of age < 15 ka), the far field LGM, 15-30 ka, and the far field MIS 3, 30-60 ka.
The old far field group held every age, and its error model, fitted on the residuals of ens-v2's best members (dominated
by the Holocene points), was 23 m over 100 km: wide enough to hide that every run is 16-21 m short of the far-field index
points at 18-30 ka (`kryos-work/esl-records`). The two older groups hold:

- the primary index points (LGM 568, MIS 3 82);
- the limiting samples (LGM 97, MIS 3 49), one-sided (`pm_compare.jl` `limit_misfit`): d = 0 on the side the sample
  allows (marine: sea level above it; terrestrial: below), d = m − rsl beyond it, σ as an index point's. A sample that
  holds adds 0 to χ² and 1 to n, so every series is scored on the same points. Kikaijima's 46 marine-limiting corals of
  38-60 ka are left out (`LIMITING_EXCLUDE`): every series, PaleoMIST's own GIA and every ESL record included, has sea
  level 63-78 m below them, so they measure the site's uplift correction (GAPSLIP's primary rate; its alternatives of
  1.7-2.3 mm/yr are not primary data), not the ice; Kikaijima's index points stay;
- the two proxy records of GAPSLIP, FesmData `RSL-GAPSLIP_proxy` (FesmData #16, released in RSL v2.3.0;
  `[score] rsl_proxy` to read another file), as site RSL by age: the Red Sea at Bab-el-Mandeb (Grant et al. 2012, 2014; index
  points, ±12 m), thinned to the mean of each 1-kyr bin (`read_proxies`: a continuous record of 4-5 samples per kyr whose
  δ18O scatter is not independent; LGM 6 points, MIS 3 30), and the Karimata Strait (Weiss et al. 2022; terrestrial
  limiting, 60-32 ka, every 2 kyr: 15 in MIS 3).

Far-field samples older than 60 ka are in no group. **MIS 3 is reported, not judged** (`[ensemble] report_only`): its
data are few (82 index points, most at Huon and Cairns), the limiting samples and the Red Sea carry it, and the ice's
2.5-kyr records miss its millennial changes; every series misses it (C7 92, PaleoMIST 61) mostly by too little MIS 3 ice
(C7 26-42 m above the Red Sea). The far field LGM is judged.

The error model of both is a priori, not fitted on members: a fit to members that are all short of the LGM ice would take
their common deficit for model error, which is what hid it. σ_m 6 m: the GIA of a good Earth at a far-field site (site
RSL less ESL, −9 to +6 m at 30-20 ka over Earths and ice histories) and the tectonic uplift corrections (Huon, Barbados:
0.1-0.3 mm/yr over 20-60 kyr); L_x 2000 km, that term being of ocean-basin scale (ocean syphoning, continental levering);
L_t 5 kyr, the 2.5-kyr ice records. Under it the LGM's I is high for every series (C7 158.7, `pub-PMIST` 119.7,
PaleoMIST's ice on E117 121.8, JE-lo 457.8), and the bar is set accordingly ([below](#the-bar)). The index-point RMSE of
the far-field sites, the metric of the first comparisons, is kept for context (`sites.csv` of a score folder,
`plastic_ensemble.jl` `write_sites`; [results](results.md#scores-v5)).

### GNSS in new regions

Two GNSS groups from Schumacher et al. (2018) (FesmData `GPS-Schumacher2018`, global, elastic-corrected), each a box
(`[[gnss]]`): **North America forebulge**, 35-50°N, 95-65°W, south of the Hudson Bay box (the Great Lakes, the St
Lawrence, New England and the US Atlantic coast: the collapsing peripheral bulge, sensitive to the lower mantle and the
total Laurentide load), 661 stations; **British Isles**, 49-61°N, 11°W-2°E (a small ice sheet on a thin lithosphere),
112 stations. MIDAS (observed, not corrected) is left out there: Schumacher's stations already cover both regions, and
MIDAS's hydrology and groundwater signals (the US coastal plain) are not corrected. The error models are fitted as the
other GNSS groups', from the 6 best Earths of wave A of the Earth search (`plastic_ensemble.jl errmodel
refs=E1D-ES1-E`; `kryos-work/scores-v4/errfit`): North America forebulge 1.5 mm/yr over 750 km, British Isles 0.7
mm/yr over 175 km.

### Antarctic ice cores

Not added. Buizert et al. (2021) do not give their elevation changes as data: Data S1 of the supplement
(`abd2897_datas1.xlsx`, 564 KB) holds the firn-based temperature and accumulation of the seven cores (0-40 ka), the total
air content of WAIS Divide and Dome Fuji, δ15N, CH4, borehole temperatures and Δage, no elevation (checked 2026-11-18);
an elevation history would have to be derived from the total air content (two cores); the open NOAA archive of the study (doi:10.25921/hnt8-bv54) has
temperature and accumulation, no elevation. The main text has only the West-minus-East contrast (WAIS Divide 420 m,
280-590 m, higher than central East Antarctica). Skytrain Ice Rise (Grieman et al. 2024, CC BY 4.0) has a 430 ± 110 m
lowering between 9.5-8.5 and 6.6-5.0 ka (total air content), and no peer-reviewed Holocene elevation numbers were found
for Roosevelt Island or Siple Dome (flow-model estimates only, paywalled).

## The error model

Per data group g, the residuals d = model − observation have the covariance

R = diag(σ_i²) + s_i s_j C,  s_i = √(σ_m,g² + sm_i²),  C_ij = exp(−dist_ij / L_x,g) · exp(−|Δt_ij| / L_t,g)

with σ_i the data errors and a model error σ_m,g correlated over distance L_x,g and time L_t,g; sm_i is a point's own part
of the model error, the spread of the ESL records at its age ([ESL](#esl)), 0 in every other group. The model error
covers what no simulation can fit (the resolution, the ice outside the reconstruction, the solver). Its parameters are
fitted, then held fixed for every simulation:

- the RSL groups from the residuals of the best members of ens-v2, but the far field (below 15 ka since scores v4) from
  the imposed ice histories: the 6 best of PaleoMIST's and GLAC-1D's ice on the 200 Earths of wave A of the Earth search
  (`plastic_ensemble.jl errmodel refs=...`, the RSL groups from references since scores v4; `kryos-work/scores-v4/errfit`):
  16 m over 50 km and 20 kyr (all six PaleoMIST's; until scores v4 23 m over 100 km and 20 kyr, from ens-v2, every age).
  An error model fitted on our members takes their common misfit for model error; on the imposed ice it is what the
  best Earths cannot fit under a published ice;
- the far field LGM and MIS 3 a priori, 6 m over 2000 km and 5 kyr ([above](#far-field-lgm-and-mis-3));
- GNSS from wave A of the Earth search: Fennoscandia 0.7 mm/yr over 750 km, Hudson Bay 1.5 mm/yr over 500 km, Greenland
  2.0 mm/yr over 1000 km, Antarctica 4.0 mm/yr over 250 km, North America forebulge 1.5 mm/yr over 750 km, British Isles
  0.7 mm/yr over 175 km;
- ESL σ_m 5 m over 10 kyr with the records' spread; the rotation group a model error per quantity;
- the ice cores 75 m over 500 km and 5 kyr, the starting values ([above](#ice-cores)), kept: refit on the 6 best members
  of the joint ensemble's wave 1 (`errmodel refs=^tau-JE-`, `kryos-work/joint-ens/wave1/score-v4-corefit`) the maximum
  likelihood is on the edge of the grid (80 m, the largest; 10 km, the smallest) with χ²/n still 2.4: the best members'
  misfit is not of the form of the error model (the plastic surface on 32 km misses the Holocene thinning at single
  sites), and 20 points do not fit three parameters.

## The scores

- **log L**, the Gaussian log-likelihood −½ d′R⁻¹d − ½ log|R| summed over the RSL groups (seven since scores v4, the
  far field split in two): higher is better. Differences of a few tens are large. This is the one scale of all the tables.
- **I**, the implausibility of a group, (χ² − n)/√(2n) with χ² = d′R⁻¹d: about 0 when the residuals are as large as the
  error model expects, ≤ 3 passes.
- **passed**, the number of judging groups within the bar: of 14 (seven RSL groups, ESL, five GNSS groups, rotation; of 16
  in scores v4, of 13 until then, of 9 until 2026-11-18), over the groups a series has a score in: a published
  reconstruction has no GNSS and no rotation, `pub-GLAC1D` no far field LGM, a run before 2026-11-12 no GNSS. Reported,
  not judged: the far field MIS 3, the Greenland ice cores, GNSS British Isles (scores v5: no member of the joint
  ensemble passes the last two). In the Earth search ESL is reported, not judged.
- **RMSE** [m] and **RMS z** (the root mean square of the residuals normalized by σ_i, the score of the first
  comparisons) are reported for reference only. A few points with tiny σ_i dominate RMS z, Greenland above all.

## The bar

A group passes at I ≤ 3. Where both published reconstructions fail a group, the cut is raised to the better of their
two scores: Greenland 7.9 (PaleoMIST) and, since scores v4, the far field below 15 ka 3.3 (PaleoMIST; GLAC-1D 4.4).
The far field LGM has bar references of its own (`[ensemble] bar_groups`, `plastic_ensemble.jl` `bar_refs`,
`group_cut`): `pub-PMIST` and, for GLAC-1D, whose own RSL ends at 26 ka, its ice on our solver with Gowan's 1-D Earth
(I5, `ice-GLAC1D | E1D-PMIST`, 184.5); the cut is PaleoMIST's 119.7. Under its a priori error model every series misses
that group far beyond 3, the publication too: the cut asks for an LGM at the far-field sites no worse than PaleoMIST's
(C7 158.7). Every other RSL group keeps 3. No publication has a GNSS or rotation score: those groups keep 3.

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
