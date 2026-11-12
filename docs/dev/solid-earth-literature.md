# Solid-Earth structure for GIA: a literature review for kryos

This review supports one decision: which Earth structures, and which parameter ranges, to use for the next kryos Earth ensemble. That ensemble is an Earth basis built from the principal components of expert 3-D viscosity fields, with the lithosphere as an extra parameter and GNSS uplift rates as a second constraint. Our runs couple the plastic ice-sheet model to VILMA2 (3-D viscosity, degree 170) and score them against GAPSLIP 2.1 RSL index points in five groups (see [vilma2-paleomist.md](vilma2-paleomist.md) and [simulations.md](../simulations.md)).

Each source is marked by how much of it was read:

| tag | meaning |
|---|---|
| **FT** | full text read |
| **P** | part of the full text read |
| **A** | abstract only |
| **S** | secondary: numbers as cited by another paper, not checked in the original |

Every DOI in the reference list was checked against Crossref. Text marked *Inference* is our reading, not a statement of the cited paper.

## Summary for our decisions

1. **Pan et al. (2022) is a far-field study, not a study of near-field fit.** It quantifies how 40 realizations of 3-D Earth structure change the predicted LGM lowstand at Barbados, the Great Barrier Reef, the Bonaparte Gulf and the Sunda Shelf, relative to their 1-D spherical means. It does not fit near-field RSL, does not compute present-day uplift and does not rank structures against data.
2. **3-D structure moves far-field LGM RSL by a few metres, with a site-dependent sign.** At Barbados the 3-D lowstand is 2.5–11 m shallower than the 1-D one in all 40 runs (mean ~6.5 m). At Noggin Pass (Great Barrier Reef) it is 0.5–8 m deeper (mean ~2.9 m). At the Bonaparte Gulf it is about zero (−3.5 to +2.1 m).
3. **Our far-field misfit cannot be removed by Earth structure.** Our model is 10–20 m too high in the far field. Pan et al. find no consistent far-field signal, and at Barbados the 3-D effect has the wrong sign: it makes the predicted lowstand shallower still. The misfit points to too little LGM ice; PaleoMIST holds ~116 m of GMSL against 125–135 m in far-field GIA studies. Part of it is also a VILMA2 setup difference of about +4 to +8.5 m against Gowan's SELEN under the same ice.
4. **In Pan et al., the choice of seismic model and of the 1-D background profile matter most.** The lithosphere model at fixed global mean, and the temperature-to-viscosity scaling, matter least in most far-field regions. Halving the scaling halves the logarithmic range of lateral viscosity. *Inference:* the geometry of the anomalies (which tomography) and the radial profile are first-order parameters for our basis, and the amplitude is a second, separable one.
5. **Our current structure is Bagge's extreme member.** Bagge structure 6, v_1.0_s16, has the largest lateral contrast and the thickest cratonic lithosphere of the 18. Bagge's reference is v_0.4_s16, chosen from a rheological argument, not from a fit to data. Bagge et al. (2021) already show that 3-D structure gives less uplift than its 1-D mean at Ångermanland, and that the structures spread RSL there by ~150 m at 14 ka. Our sign flip around the Gulf of Bothnia fits this picture.
6. **3-D Fennoscandia fits worse with an ice history built on a 1-D Earth.** Blank et al. (2025) find that 3-D models only match the 1-D fit with ICE-6G, which was tuned with VM5a. With ice histories not tuned to a 1-D Earth they beat 1-D models, using a stiff craton of 1–5×10²¹ Pa s and a 122–189 km lid. PaleoMIST was built on Gowan's 1-D Earth, so the same co-adaptation most likely penalizes our 3-D runs.
7. **Regional 1-D values span less than our Earth ensemble does.** The literature gives a lithosphere of 65–160 km, an upper mantle of 3–7×10²⁰ Pa s and a lower mantle of ~10²¹–10²³ Pa s for Fennoscandia and North America. Greenland RSL prefers 120 km, 5×10²⁰ and a lower mantle of 1–10×10²¹ Pa s. The far field prefers a thin lithosphere (40–70 km), ~1.5×10²⁰ in the upper mantle and a bimodal lower mantle of ~2×10²¹ or ~7×10²² Pa s. The Barents Sea is softer than Fennoscandia, with half its upper-mantle viscosity and a 65–75 km lithosphere.
8. **The west Greenland and St. Lawrence misfits are ice-history problems.** Milne et al. (2018) find in a 3-D ensemble that lateral structure reduces the RSL fall at Disko Bugt, and conclude it cannot explain the misfits. ICE-6G and ICE-7G over-predict St. Lawrence RSL by up to ~60 m at ~9 ka, as our runs do, and thinning the local ice stream reduces this.
9. **No published GIA study uses a PCA basis of expert 3-D fields.** The established approaches are 1-D grid or Bayesian inversions with regional ice scaling, factorial ensembles of 3-D structures, and a scalar scaling of tomographic anomalies. Our basis is therefore new. It should keep the lid definition, the radial profile and the anomaly amplitude as explicit parameters, so that the principal components only carry the anomaly geometry.
10. **GNSS constrains the upper mantle where elastic and recent loads are small.** Fennoscandia (data errors ~0.2–0.6 mm/yr) and Hudson Bay are the useful regions. In Greenland and Antarctica, present-day uplift is dominated by elastic response, Little Ice Age loading and transient (low-viscosity) rheology. There it needs an elastic correction with ~1 mm/yr or more of uncertainty, and inflated model errors.

## 1. Pan et al. (2022)

Pan, Milne, Latychev, Goldberg, Austermann, Hoggard and Mitrovica (2022), *The influence of lateral Earth structure on inferences of global ice volume during the Last Glacial Maximum*, QSR 290, 107644. **FT** (accepted manuscript, 37 pp., including Tables 1–3; Figs 1 and 3 inspected).

### 1.1 Aim

Pan et al. ask whether 3-D Earth structure biases LGM ice volume when that volume is inferred from far-field RSL. Their motivation is the "missing ice" problem: far-field GIA studies need 125–135 m of GMSL-equivalent ice at the LGM, while the ice sheets reconstructed region by region sum to less. They list ICE-6G at ~127 m (~14 m in Antarctica, tuned to Barbados) and Lambeck et al. (2014) at ~134 m (~23 m in Antarctica). They also list Gowan et al. (2021), PaleoMIST, at 116 m, "reported to match observed LGM RSL lowstands". The study extends Austermann et al. (2013), who found that a high-viscosity slab beneath Barbados requires ~7 m more LGM ice than a 1-D Earth.

### 1.2 Earth structure

Each 3-D Earth is a 1-D background profile with lateral variations added to it. Pan et al. vary four ingredients (their Tables 1–2).

**Background profiles** (spherical means of the 3-D Earths):

| label | ice history | lithosphere | upper mantle (Pa s) | lower mantle (Pa s) |
|---|---|---|---|---|
| M1DA | ANU (Lambeck et al. 2014) | 71 km | 2×10²⁰ | 3×10²² |
| M1DB | ANU | 71 km | 3×10²⁰ | 2×10²¹ |
| 96VM2 | ICE-6G | 96 km | 4×10²⁰ | 2.2×10²¹ (top ~1100 km), 3.3×10²¹ below |

M1DA and M1DB sample Lambeck et al.'s two solution classes, with the upper-to-lower-mantle jump above two orders of magnitude and about one, respectively. Both lie within Lambeck's 1-σ ranges. For M1DA they chose a lower mantle of 3×10²², not Lambeck's optimal 7×10²², because "a growing number of studies" favour the smaller value. The third profile is a three-layer version of VM2, which Pan et al. call very similar to VM5a.

**Seismic models:** S40RTS (Ritsema et al. 2011), Savani (Auer et al. 2014), SEMUCB-WM1 (French & Romanowicz 2014) and SL2013sv (Schaeffer & Lebedev 2013). The same four were compared by Milne et al. (2018).

**Conversion** (following Latychev et al. 2005, with the parameters of Austermann et al. 2013):

δlnρ = (∂lnρ/∂lnv_s)(r) δlnv_s,   δT = −δlnρ / α(r),   η = η₀(r) exp(−ε δT),

where α is the depth-dependent thermal expansivity, ∂lnρ/∂lnv_s a depth-dependent velocity-to-density scaling and ε the temperature sensitivity of viscosity. They use ε = 0.04 or 0.02 °C⁻¹. They state that halving ε roughly halves the logarithmic range of lateral viscosity over part of the mantle. No composition term and no viscosity cap or floor are described.

**Size of the lateral variations** (Fig. 2, S40RTS and Savani, ε = 0.04):

- At 346 km, 99% of cells lie within roughly ±2 orders of magnitude of the background.
- At 1071 km, 99% of cells lie within −2.9 to +3.9 orders (S40RTS) and −3 to +3.4 orders (Savani).
- The authors regard these lower-mantle ranges as an upper bound, since mineral physics suggests less.
- High viscosity follows subduction zones (Malay Archipelago).
- The models can disagree in sign, for example beneath Québec at 346 km.

**Lithosphere**, an elastic layer of laterally varying thickness:

- AF (Afonso et al. 2019) is thermal, from a joint geophysical–geochemical inversion.
- YO (Yousefi et al. 2021) takes elastic thickness from gravity–topography coherence on the continents (Audet & Bürgmann 2011; Chen et al. 2017; Steffen et al. 2018 for Greenland), and a sea-floor-age geotherm in the oceans.
- Both are scaled to the global mean of the background profile (71 or 96 km).
- At a 71 km mean, the minimum is ~15 km (AF) and ~25 km (YO) at ridges, and the maximum ~200–350 km in cratons (Fig. 1). The authors note the cratonic maxima are partly an artefact of that scaling.

**Ensemble:** 2 ice histories × backgrounds × 2 lithospheres × 4 seismic models × 2 values of ε. M1DB is run with ε = 0.04 only. This gives 40 3-D runs (24 ANU, 16 ICE-6G) and three 1-D runs.

### 1.3 GIA model

The model is Seakon, the finite-volume code of Latychev et al. (2005), with the following set-up:

- **Physics:** a Maxwell viscoelastic Earth. It solves the generalized sea-level equation (Kendall et al. 2005) with time-varying shorelines and rotation (Mitrovica et al. 2005). Paleotopography takes three outer iterations.
- **Grid:** ~17 million nodes and 67 radial layers. Resolution is ~12 km at the surface and ~60 km at the base of the mantle.
- **Start:** runs start at 36 ka. A test from 80 ka changed far-field LGM RSL by only ~0.1 m.
- **Cost:** several days on ~100 cores per run.

### 1.4 Results

**Maps (Fig. 3):** the mean 3-D minus 1-D LGM RSL is similar for both ice histories at low latitudes. A positive difference means a shallower 3-D lowstand, so a 1-D analysis there underestimates ice volume. The difference is positive in the western North Atlantic and Caribbean and negative over most of the South Pacific and Southern Ocean. The colour scale saturates at ±5 m. The standard deviation exceeds 3 m at high northern latitudes and around Antarctica.

**By latitude (Fig. 4):** between ~20°N and ~35°S the peak mean difference is ~4–12 m. It reaches a few tens of metres poleward of ~20°N (the Laurentide peripheral bulge) and ~35°S (West Antarctica). Pan et al. attribute the northern signal to thick lithosphere and high shallow-mantle viscosity under North America and the western North Atlantic, which suppress the bulge. They attribute the southern signal to thin oceanic lithosphere and low West Antarctic viscosity. Large low-latitude differences cluster near subduction zones (Southeast Asia, Caribbean), on the northern Australian margin, in the Red Sea and at the Central Indian Ocean Triple Junction. At the peak locations, the ensemble standard deviation is ~1–4 m, i.e. 30–50% of the signal.

**Sites (Table 3, Fig. 5):** 3-D minus 1-D RSL at the LGM (26 ka for ICE-6G, 21 ka for ANU), in m. Group mean ± standard deviation:

| background (runs) | Barbados | Sunda Shelf | Bonaparte Gulf | Noggin Pass | Hydrographer's Passage |
|---|---|---|---|---|---|
| M1DA (16) | 4.98 ± 1.66 | −1.11 ± 1.61 | 0.29 ± 1.00 | −1.95 ± 0.95 | −0.85 ± 1.09 |
| M1DB (8) | 8.94 ± 1.34 | −3.35 ± 2.65 | −2.03 ± 1.16 | −5.27 ± 1.80 | −2.92 ± 1.74 |
| 96VM2 (16) | 6.89 ± 2.41 | −1.95 ± 2.02 | 0.32 ± 0.92 | −2.53 ± 1.31 | −1.42 ± 1.22 |
| all (40) | 6.54 ± 2.42 | −1.90 ± 2.12 | −0.16 ± 1.36 | −2.85 ± 1.78 | −1.49 ± 1.47 |
| range | 2.6 to 11.0 | −6.7 to 1.8 | −3.5 to 2.1 | −8.2 to −0.5 | −5.3 to 1.1 |

- **Barbados:** every run gives a shallower 3-D lowstand. Pan et al. attribute this to the subducting slab, which reduces the subsidence under the ocean load, and/or to altered peripheral-bulge dynamics, consistent with Austermann et al. (2013).
- **Noggin Pass:** every run gives a deeper lowstand. The local lithosphere is thin in both models (~38 km AF, ~56 km YO), so continental levering amplifies the offshore subsidence.
- **Sunda Shelf:** most runs give a deeper lowstand. The lithosphere there is 40–45 km, and SEMUCB-WM1 is ~1 order of magnitude softer than S40RTS at 346 km.
- **Bonaparte Gulf:** either sign, small on average.
- **All sites together:** the mean effect over the five sites and all runs is ~0 (mean 0.03 m, standard deviation ~3.9 m).
- **Lower-mantle contrast:** the M1DB group, with its smaller upper-to-lower-mantle contrast, gives the largest effects at every site.
- **Scaling:** halving ε halves the Barbados effect (2.6–4.4 m against 5.9–7.1 m in the M1DA group).

**Sensitivity (Fig. 6):** Pan et al. change one ingredient of a reference run (S40RTS, AF, ε = 0.04, M1DA) at a time.

- Changing the lithosphere model, at the same global mean, has the smallest effect in most regions. The exception is near mid-ocean ridges, where YO gives lower RSL.
- Changing ε also has a small effect in most regions. It mainly rescales the pattern, which they suggest reflects the long wavelengths of S40RTS.
- The seismic model and the background profile have the largest effects.

All four ingredients matter at some far-field locations (Fig. S4). They conclude that 3-D structure gives "no consistent, high magnitude signal" across sites and models, so it does not resolve the missing-ice problem. They still recommend 3-D Earth models and multiple sites for far-field ice-volume estimates.

### 1.5 What Pan et al. do not show

Pan et al. show no near-field RSL, no present-day uplift and no misfit to data. Their maps mask the near field, and their discussion stops at 45°S–45°N. They therefore say nothing direct about Fennoscandia, Greenland or North American uplift.

The same group and code (Seakon, the same four tomographies and the AF/YO lithospheres) cover the near field in other papers:

- **Greenland (Milne et al. 2018, P):** lateral structure changes Holocene RSL by several tens of metres. It reduces the RSL fall at Disko Bugt in all runs, and it is "not a viable explanation" for the east–west misfit pattern. It changes uplift by up to >2 mm/yr, the ensemble spread is as large as the mean effect, and it does not improve the fit at 54 elastic-corrected GNET sites.
- **Pacific coast of North America (Yousefi et al. 2021, A):** lateral structure changes RSL significantly but does not improve the data–model fit.
- **Eastern North America (Parang et al. 2024, A):** 3-D models support a split of the region into two Earth domains but do not improve on the regional 1-D fits.
- **North America, MIS 5a (Thompson et al. 2023, P, via a separate agent):** with the Richards/Hoggard/Austermann 3-D model, the χ² misfit falls from 12.6 (1-D) to 2.0 (3-D) on the Atlantic–Caribbean coast and from 5.1 to 3.5 on the Pacific coast.
- **Greenland, transient rheology (Pan et al. 2024, P):** a thin low-viscosity layer below the lithosphere (45 km at 10¹⁹ or 85 km at 4×10¹⁹ Pa s) barely changes Holocene RSL but raises present-day uplift by an order of magnitude or more.

### 1.6 The Pan field in FesmData

FesmData holds one Pan et al. lithosphere map and one 3-D viscosity field (`Pan2022_litho`, products `Pan2022-litho` and `Pan2022-visc`, not released). These came by personal communication via FastIsostasy and are archived in isostasy_data v0.2.0. The field has log10 η on 64 radii from 24 to 2841 km depth, at 0.5°.

Pan et al. publish no data, and neither the archive nor Swierczek-Jereczek et al. (2024) states which of the 40 realizations it is. Swierczek-Jereczek et al. state only that their Antarctic test uses "the lithospheric thickness and the mantle viscosity fields from Pan et al.". Their FastIsostasy comparison against Seakon (Test 4, a 122 kyr Antarctic run) gives a maximum relative error of ~0.14 for their laterally variable model and ~0.22 for a 1-D Earth.

*Inference:* before using the field, identify its background from its spherical means:

- a 71 km lithosphere with ~2×10²⁰ / 3×10²² Pa s means M1DA;
- 71 km with ~3×10²⁰ / 2×10²¹ means M1DB;
- 96 km with ~4×10²⁰ / ~2–3×10²¹ means 96VM2.

The ε value follows from the log range at ~350 km: about ±2 orders for 0.04. The tomography is not identifiable from statistics alone and needs the authors.

Seakon treats the lithosphere as an elastic layer of the mapped thickness. The natural VILMA2 lid for this field is therefore `z_litho` itself, not a viscosity threshold applied to `log10_visc`.

### 1.7 Implications for our results

- **Far field.** Our runs are 10–20 m too high at Barbados, the Great Barrier Reef and the Bonaparte Gulf. On Pan et al.'s numbers, 3-D structure would add 2.5–11 m at Barbados (worse), remove 0.5–8 m at Noggin Pass and ~1.5 m at Hydrographer's Passage (better), and leave the Bonaparte Gulf unchanged. A uniform +10–20 m offset is therefore not an Earth-structure signal.
  - *Inference:* it is most likely too little LGM ice, since PaleoMIST's 116 m sits below the 125–135 m of far-field GIA studies.
  - *Inference:* part of it is also the VILMA2 ocean bias of +8.5 m against Gowan's SELEN at 22.5 ka under the same ice ([vilma2-paleomist.md](vilma2-paleomist.md), finding 2), which our own setup must resolve.
- **Lower mantle.** In Pan's runs, the smaller upper-to-lower-mantle contrast (M1DB) gives the largest 3-D effects. Lambeck et al. (2014) find that their low-viscosity solution needs ~7 m more ice.
  - *Inference:* lowering our lower mantle from Gowan's 4×10²² towards ~2×10²¹ Pa s would make the far field shallower for the same ice, which is the wrong direction for our misfit.
- **North American bulge.** Pan et al. find that 3-D structure suppresses the Laurentide peripheral bulge in the western North Atlantic by up to a few tens of metres at the LGM. *Inference:* sites on the US Atlantic margin and the Caribbean are where a stiff North American lithosphere should change our fit.
- **First-order ingredients.** Pan et al. identify the seismic model and the background profile as first-order, and the lithosphere map (at fixed mean) and ε as second-order. *Inference:* a basis built from Bagge's 8 structures alone samples ε-like amplitude and radial-profile variants of a single tomography (SL2013SV over Grand). It therefore misses the ingredient Pan et al. find most important. Adding Lloyd (GLAD-M25) and Pan (one of four tomographies) gives the basis different anomaly geometries.

## 2. Bagge et al. (2021) and the GFZ structures

Bagge, Klemann, Steinberger, Latinović and Thomas (2021), G-Cubed 22, e2021GC009853. **FT** (main text via a separate agent; supporting information not read).

### 2.1 How the 18 structures are built

**Tomography and conversion:**

- The structures follow Steinberger (2016). They use SL2013SV above 200 km and the 2010 update of Grand (2002) below, with a smooth transition at 200 km.
- Velocity becomes temperature as in model M2 of Steinberger & Calderwood (2006). In the upper mantle, −(∂v_s/∂T)/v_s ≈ 10–15×10⁻⁵ K⁻¹, i.e. a −1 to −1.5% anomaly corresponds to +100 K.
- Viscosity follows η ∝ exp(rH/RT), with H ≈ 500–700 kJ/mol in the upper mantle (525 kJ/mol for dry olivine dislocation creep).

**The reduction factor r** takes the values 0.2857, 0.3, 0.4, 0.5, 0.75 and 1.0.

- It mimics non-Newtonian (dislocation) creep with a Newtonian law of reduced activation enthalpy (Christensen 1983). r = 1 corresponds to diffusion creep.
- A larger r gives a stronger temperature dependence, i.e. larger lateral contrasts and a thicker lithosphere.
- r = 0.4 is the reference. It was chosen because Christensen (1983) showed that n = 3 power-law flow is imitated by a reduction of 0.3–0.5. The choice does not come from a fit to data.

**Radial viscosity profiles**, constrained by the geoid, the radial heat flux and a Haskell average of 10²¹ Pa s:

- **s16** has a viscosity increase at 410 km, which GIA profiles such as VM2 and VM5a do not have (class I).
- **sc06** has a softer transition zone than the layers above and below (class II).
- **sc06b** has no step at 410 km, the closest to the homogeneous upper mantle usual in GIA (class III).
- The 18 structures are named v_[r]_[profile].

**Lithosphere:** the depth below which η < 10^23.5 Pa s, i.e. ~1087 K at 100 km and 1153 K at 200 km for v_0.4_s16. The authors stress that this is not the effective elastic thickness of GIA studies.

**Numerical limits:**

- Viscosity is floored at 10¹⁹ Pa s for cost, so West Antarctica is not fully resolved.
- Lateral variations are kept only above 870 km; deeper layers hold their depth mean.
- The fields are resolved to degree 63 (~320 km).

**Ranges across the ensemble:**

- The lithosphere is <30 km with an upper mantle of ~10²⁰ Pa s in Cascadia and West Antarctica.
- It is >200 km with >10²³ Pa s in the "asthenosphere" layer (lithosphere base to 225 km) in the cratons of Fennoscandia, North America and East Antarctica.

FesmData holds 8 of the 18 structures. They are the six s16 structures plus v_0.4_sc06 and v_0.4_sc06b, the same 8 that Bagge et al. ran with the NAICE ice history. Our structure 6 is v_1.0_s16, the largest r. *Inference:* it has the stiffest cratons, the thickest lithosphere and the s16 step at 410 km, and so the stiffest Fennoscandian mantle of the set.

### 2.2 What they compared and found

**Set-up:**

- VILMA: incompressible, PREM elasticity, degree 170, explicit 2-yr time step and four spin-up loops over the glacial history.
- Ice: ICE-5G and ICE-6G with all 18 structures, NAICE with 8. Together with 1-D global-mean and regional-mean variants, this makes ~300 runs.
- RSL is compared visually at six sites: Churchill, Ångermanland, Ross Sea, San Jorge Gulf, Central Oregon and Pioneer Bay. There is no misfit statistic, no GNSS and no ranking of structures.

**Spread across the structures:**

- At 14 ka the RSL spread is ~150 m at Churchill (ICE-5G), ~100 m (ICE-6G) and ~50 m (NAICE). It is ~150 m at Ångermanland and ~5 m in the far field.
- At ice-centre sites a higher r gives lower RSL.
- sc06, with its soft transition zone, gives the highest RSL at Churchill and Ångermanland, and s16 the lowest.

**3-D against 1-D:** under Europe the 3-D models give less uplift than their 1-D global means, because the lithosphere is thick and the mantle stiff. The Ross Sea shows the opposite.

**Regional 1-D means:** they reproduce the 3-D result under North America (e.g. within ~15 m at 11 ka for m_1.0_s16), but not in Oregon, Patagonia or Antarctica.

**Data fit:** 3-D structures fit Oregon and Patagonia better than VM5a (EGU21 poster, via a separate agent).

### 2.3 Related GFZ work

- **Albrecht, Bagge & Klemann (2024, P)** couple PISM to VILMA with v_0.4_s16 as reference.
  - They describe that structure as optimized to minimize the misfit to a global RSL compilation. *Inference:* Bagge et al. (2021) do not show such an optimization, so this may refer to unpublished work.
  - In Antarctica, the structure has a lithosphere <50 km in the west and >130 km in the east, and 10¹⁹–10²¹ Pa s at 280 km.
- **Schachtschneider et al. (2022, P)** assimilate synthetic RSL into a 1-D VILMA with a particle filter. They recover the upper mantle to ~1% and the lower mantle to ~5%. Laurentide data constrain the lower mantle better than Fennoscandian data, and ice-model errors strongly bias the upper-mantle estimate.
- **Tanaka et al. (2011, A via search summary):** compressibility barely changes vertical rates but roughly doubles horizontal rates around Fennoscandia and the Laurentide. The standard VILMA is incompressible.

No GFZ paper was found that scores the 18 structures against RSL or GNSS. No recommended structure based on data exists in the literature we read.

## 3. Seismically derived viscosity and its uncertainty

**Richards et al. (2020), JGR (P, pp. 1–2).** They calibrate the Yamauchi & Takei (2016) anelasticity parameterization to convert shear-wave speed into temperature, density and diffusion-creep viscosity. The calibration fits four constraints: an oceanic plate geotherm, an adiabat at 225–400 km, the QL6 attenuation profile and a mean viscosity of ~3×10²⁰ Pa s at 225–400 km, as described by Lloyd et al. (2024). They find locally a ~150 km thick low-viscosity asthenospheric channel. Lloyd, Austermann and Hazzard all build on this calibration.

**Austermann et al. (2021), GJI (P, via fetch summaries; numbers to be checked before quoting).**

- **Construction:** SL2013sv (with SL2013NA for North America) above 400 km and SEMUCB-WM1 below, blended over 300–400 km. The upper mantle uses the Richards et al. calibration. The lower mantle uses a Perple_X conversion for pyrolite with the Q5 anelasticity correction, and viscosity follows the activation-enthalpy profile of Steinberger & Calderwood (2006).
- **Background:** 5×10²⁰ / 5×10²¹ Pa s, with the lithosphere base at the 1175 °C isotherm (~100 km mean, ~300 km under cratons, <25 km at ridges).
- **Uncertainty:** the inverse calibration constrains log10 η to within ~0.2 (one standard deviation) over 1000 parameter sets. A forward mapping without calibration spreads it about five times more. Including anelasticity raises viscosity by 1–1.5 orders of magnitude in the slowest regions, compared with an anharmonic conversion.
- **Results:** for Last Interglacial sea level, 3-D minus 1-D is metres to tens of metres in the near field and up to a few metres in the far field. A 1-D profile can match the 3-D result at a site, but it is not the local average of the 3-D viscosity beneath that site.

**Lloyd et al. (2024), GJI (P, pp. 1139–1149).**

- **Field:** viscosity from the Voigt-average v_s of GLAD-M25, built following Austermann et al. (2021). The calibrated parameters are μ₀ = 80.82 GPa, log10 η_r = 23.301 and E_a = 545 kJ/mol.
- **Lithosphere:** the 1175 °C isotherm (~100 km mean), assigned ~1.5×10²⁶ Pa s, with no explicit elastic lid.
- **Means:** the upper mantle is shifted to 5×10²⁰ and the lower mantle to 5×10²¹ Pa s.
- **Simulation:** the minimum viscosity is truncated at 2×10¹⁹ Pa s, which allows a 1-yr time step instead of 0.05 yr. Degree 64, ICE-6G.
- **3-D against 1-D:** peak RSL fall over 26 kyr is ~800 m in 1-D, against ~700 m (Laurentide) and ~500 m (Fennoscandia) in 3-D. The far-field ocean is up to 10 m higher in 3-D.
- **Kernels:** sensitivity moves towards weak regions and vanishes in stiff (effectively elastic) ones. Far-field RSL is sensitive to the deep mantle, so 1-D inferences from far-field data can be biased.
- *Inference:* this field puts the Fennoscandian shield in the same stiff class as Bagge's high-r structures.

**Antarctica.**

- **Hazzard et al. (2023, JGR; P).** A Bayesian version of the Richards calibration applied to ANT-20 (~100 km resolution) gives the following:
  - Steady-state viscosity reaches a lower threshold of ~10¹⁹ Pa s at 150–350 km under West Antarctica, and 10^(19.5±0.3) Pa s at 150 km towards the Amundsen Sea.
  - The lithosphere–asthenosphere boundary is 63 ± 13 km in West Antarctica and 233 ± 41 km in East Antarctica.
  - Uncertainty within West Antarctica stays within one order of magnitude at 0–400 km.
  - GNSS-derived viscosities of 10^18.4–10^19.4 (Amundsen Sea; Barletta et al. 2018) and lower on the Peninsula agree only once transient (load-duration-dependent) viscosity is used.
- **Ivins et al. (2021, Geol. Soc. Memoir 56; P).** Three mappings of ANT-20 give consistent results:
  - At 150 km: East Antarctica >10²¹ Pa s (10²²–10²³ in the cratonic core) and West Antarctica ~2×10¹⁹–10²⁰.
  - At 300 km: ~4×10²⁰ in a swath from the Wilkes coast to the Weddell Sea, and 10¹⁸–10¹⁹ under Marie Byrd Land and the coast.
  - The lateral range is about four orders of magnitude at 150 km and 1–1.5 at 550 km.
  - The northern Peninsula stays about an order of magnitude stiffer than geodetic estimates.
- **Wiens et al. (2021, Memoir 56; A).**
  - East Antarctic lithosphere reaches ~200 km, thickest 500–1000 km inland of the Transantarctic Mountains.
  - The mountains form a sharp boundary to the warm West Antarctic mantle, and slow velocities under Marie Byrd Land extend through the transition zone.

*Inference:* the Bagge floor (10¹⁹) and the Lloyd cap (2×10¹⁹) lie at about the steady-state West Antarctic value. They are adequate for glacial-cycle RSL. They are not adequate for decadal GNSS in the Amundsen Sea or the Peninsula, which reflects transient rheology.

## 4. Regional constraints

### 4.1 Fennoscandia and the Barents Sea

- **Lambeck et al. (1998, A):** RSL from northern Europe gives an effective lithosphere of 65–85 km and an upper mantle of 3–4×10²⁰ Pa s. The lower mantle is at least ten times larger. The response is "spatially fairly uniform" from Scandinavia to the British Isles.
- **Steffen & Kaufmann (2005, FT via OUP page)** use Neighbourhood Algorithm inversions:
  - Scandinavian RSL: 120 km (160 central, 100 peripheral), 4×10²⁰, and a poorly constrained 10²³ Pa s in the lower mantle.
  - BIFROST GPS: 120 km, 7×10²⁰, 10²².
  - Barents Sea: 65–75 km, 4×10²⁰, 2×10²², with a low-viscosity zone of (1–4)×10¹⁹ Pa s at 130–200 km. No such zone appears under Scandinavia.
- **Milne et al. (2004, A):** BIFROST gives 95% ranges of 5×10²⁰–10²¹ (upper mantle) and 5×10²¹–5×10²² Pa s (lower mantle), with an elastic lithosphere of 90–170 km (Milne et al. 2001). Resolution is ~200 km near the base of the upper mantle and ~700 km at the top of the lower mantle.
- **Zhao et al. (2012, S/snippet):** from GPS, 93–110 km, 3.4–5.0×10²⁰ and 7–13×10²¹ Pa s.
- **Kierulf et al. (2014) and Vestøl et al. (2019, S via Blank et al.):** 140 km / 7×10²⁰ and 160 km / 7×10²⁰ Pa s, respectively. Kierulf et al. found their 1-D model fitted vertical rates better than their 3-D model, except in northern Norway.
- **Rovira-Navarro et al. (2020, FT landing page):** GRACE requires ≥3×10²⁰ Pa s under the Barents Sea and ≥5×10²⁰ under Fennoscandia, a ratio of about two that matches seismic models. The lower mantle was fixed at 2.6×10²¹, and no lithosphere in 40–180 km was preferred.
- **Auriac et al. (2016, A):** Barents RSL prefers an upper mantle of 0.2–2×10²¹ Pa s and is insensitive to the lithosphere (71–120 km) and the lower mantle.
- **Blank, Bradley & van der Wal (2025, GJI; P)** run 42 models (7 1-D, 35 3-D from WINTERC-G or seismic temperatures with olivine flow laws), each with three ice histories, against 4090 RSL points and GNSS:
  - With ICE-6G (tuned with VM5a), 3-D models match but do not improve on 1-D.
  - With BRITICE-CHRONO and Patton et al. (2017) ice, 3-D beats 1-D, consistently in central and eastern Fennoscandia.
  - The ten best 3-D models have medians of 4.9×10²¹ (70–180 km) and 1.3×10²¹ Pa s (180–400 km) under Fennoscandia, with 122–189 km elastic thickness. Under the British Isles the medians are 0.56×10²¹ and 0.5×10²¹, with 70–115 km.
  - Dry rheology fits best.
  - A 1-D VM5a with a 3×10²⁰ upper mantle fits RSL best but predicts 6.5 mm/yr peak uplift against ~11 mm/yr observed.
  - Almost all 3-D models fail on the SW Norwegian coast and in the Danish straits, which need faster relaxation.
- **Steffen, Kaufmann & Wu (2006) and Steffen, Wu & Kaufmann (2007), S/snippet:** tomography-based 3-D viscosity changes vertical velocities by up to 7 mm/yr relative to 1-D. Lateral variation in the lower mantle barely affects GNSS.

*Inference:* the 1-D literature agrees on 3.4–7×10²⁰ Pa s and 75–160 km for Fennoscandia, which brackets Gowan's 4×10²⁰ and 120 km. Bagge's high-r cratons are an order of magnitude stiffer still. Blank et al. show such stiffness can fit, but only with ice not built on a 1-D Earth.

Our Fennoscandia/Barents group mixes two Earths, a stiff shield and a softer, thinner Barents shelf. Its score is therefore pulled in two directions by any single regional shift.

### 4.2 Greenland

- **Lecavalier et al. (2014, P)** test 243 1-D Earths (71/96/120 km, upper mantle 0.1–5×10²¹, lower mantle 1–50×10²¹ Pa s):
  - The optimal Huy3 Earth is 120 km, 5×10²⁰ and 2×10²¹ Pa s. All members of the 95% set have 120 km and 5×10²⁰, with a lower mantle of 1–10×10²¹.
  - Fits are more sensitive to the upper than to the lower mantle.
  - Fitting the east-coast data alone gave 120 km, 3×10²⁰ and 5×10²² for Huy2 (Simpson et al. 2009).
  - Persistent west-Greenland misfits remain (Kangerluarsuk, Pakitsoq, mid-Holocene RSL under-predicted). A shelf-edge LGM extent improved Disko Bugt.
  - The North American ice complex changes west and NW Greenland RSL at 16 ka by tens of metres (Fig. 9).
- **Milne et al. (2018, P)**, eight 3-D realizations:
  - On average, viscosity is higher under Greenland and lower to the east.
  - Holocene RSL changes by several tens of metres. The RSL fall at Disko Bugt is reduced in all runs.
  - Lateral structure is "not a viable explanation for the east–west issue".
- **Khan et al. (2016, FT):** GNSS and RSL with seven regional profiles give 80 km, 5×10²⁰ and 10²² Pa s (Green1). In southeast Greenland (Kangerlussuaq) they find 40 km and a 200 km asthenosphere at 10¹⁹ Pa s, attributed to the Iceland hotspot track. That result rests on two GPS sites.
- **Mordret (2018, A):** pronounced low velocities in the upper mantle beneath the east coast, the Iceland hotspot track.
- **Adhikari et al. (2021, FT) and Lewright et al. (2026, P):**
  - Present-day uplift residuals need an upper mantle of 6–11×10¹⁹ Pa s (Adhikari) or 6.3–7.0×10¹⁹ Pa s (Lewright, with a 200 km lithosphere preferred) once Little Ice Age loading is included, i.e. 4–8 times below glacial values.
  - Lewright et al. use the GNSS uplift of Berg et al. (2024, not read), 1.5–13.8 mm/yr.
- **Pan et al. (2024, P):** a thin low-viscosity layer can reconcile Holocene RSL (75 km, 5×10²⁰, 2×10²¹) with fast present-day uplift.

*Inference:* Greenland's glacial-cycle data prefer an upper mantle near 5×10²⁰ and a lower mantle well below Gowan's 4×10²² Pa s. This is consistent with our ensemble preferring stiffer Greenland mantle than Bagge structure 6 provides under the east coast. The west Greenland misfit is not an Earth-structure problem (Milne et al. 2018; Lecavalier et al. 2014).

### 4.3 North America

- **Lambeck, Purcell & Zhao (2017, A via EGU abstract):** RSL, tilted proglacial-lake shorelines and GPS give 102 (85–120) km, 5.1 (3.5–7.5)×10²⁰ and 1.3 (0.8–2.8)×10²² Pa s (95% limits). The upper mantle differs significantly from the oceanic value.
- **VM5a** (Peltier et al. 2015, A; Argus et al. 2014, P via a separate agent; Pan et al. 2022):
  - 60 km elastic over 40 km at 10²² Pa s, with 5×10²⁰ to 660 km, rising to ~3×10²¹ Pa s in the deep mantle.
  - ICE-6G was tuned with this profile.
  - The VM6 and VM7 values were not found.
- **Mitrovica & Forte (2004, A):** the joint convection–GIA inversion gives ~4×10²⁰ in the upper mantle, rising by about three orders of magnitude to >10²³ Pa s near 2000 km.
- **Parang et al. (2024, A):** 1013 RSL points at 38 eastern sites, 34 ice histories and 440 1-D Earths. Fits improve markedly when Hudson Bay/Ungava/Labrador and the Atlantic coast (including the St. Lawrence) have separate Earths. 3-D models do not improve on these.
- **Li et al. (2022, EGU abstract and slides):** ICE-6G (VM5a) and ICE-7G (VM7) over-predict St. Lawrence RSL by up to ~60 m at ~9 ka. Local thinning of an ice stream in the St. Lawrence channel reduces the misfit by ~40%.
- **Li et al. (2020, A):** across 250 1-D and 3-D models, North American present-day uplift is uncertain by up to ~2.4 mm/yr in SW Hudson Bay, mainly from 3-D lower-mantle structure.

The GLAC-1D Earth (Tarasov et al. 2012) was not found in what we read. Tarasov et al. list Earth rheology as a major unquantified uncertainty of their calibration (A).

### 4.4 Far field

- **Lambeck et al. (2014, FT):** ~1000 far-field points give a lithosphere of 40–70 km (poorly resolved) and an upper mantle of 1.5 (1–2)×10²⁰ Pa s (well constrained).
  - The lower mantle is bimodal: 2 (0.7–4)×10²¹ or ~7 (1–20)×10²² Pa s. The higher solution is adopted.
  - The low-viscosity solution needs ~7 m more LGM ice.
  - Results are nearly identical with or without Barbados.
  - Separate continental-margin and island solutions "do not require lateral variation", but islands resolve depth poorly.
- **Caron et al. (2017, P):** a Bayesian 1-D inversion of 5720 RSL records gives log10 η_UM = 20.76 ± 0.24 and a bimodal lower mantle (10^21.57 and 10^22.38), with a 104 ± 29 km lid.
  - Moving the lower mantle between the two modes trades against ~20% more East Laurentide ice and ~30% less Antarctic ice (ANU).
  - The LGM lowstand differs by up to ~15 m between the modes.
- **Austermann et al. (2013, A):** the Caribbean slab makes the post-LGM RSL change at Barbados ~10 m smaller than in 1-D, which implies ~7 m more ice (Pan et al. 2022).
- **Yokoyama et al. (2018, A):** the Great Barrier Reef lowstand reaches −118 m RSL, with a GIA-derived GMSL of ~−125 to −130 m.

### 4.5 Regional Earth parameters preferred in the literature

Values are 1-D unless marked 3-D. LM is lower mantle and UM upper mantle.

| region | lithosphere | UM (Pa s) | LM (Pa s) | source (read) | data | strength |
|---|---|---|---|---|---|---|
| N Europe | 65–85 km | 3–4×10²⁰ | ≥10× UM | Lambeck et al. 1998 (A) | RSL | UM good, LM weak |
| Scandinavia | 120 (100–160) km | 4×10²⁰ | 10²³ | Steffen & Kaufmann 2005 (FT) | RSL | LM poor |
| Fennoscandia | 120 km | 7×10²⁰ | 10²² | Steffen & Kaufmann 2005 (FT) | BIFROST GPS | moderate |
| Fennoscandia | 90–170 km | 5×10²⁰–10²¹ | 5×10²¹–5×10²² | Milne et al. 2001, 2004 (A) | BIFROST GPS | 95% ranges |
| Fennoscandia | 93–110 km | 3.4–5.0×10²⁰ | 7–13×10²¹ | Zhao et al. 2012 (S) | GPS | 3-layer only |
| Fennoscandia | 140 / 160 km | 7×10²⁰ | not read | Kierulf et al. 2014 / Vestøl et al. 2019 (S) | GNSS, levelling | model choice |
| Fennoscandia, 3-D | 122–189 km | 4.9×10²¹ (70–180 km), 1.3×10²¹ (180–400 km) | 3×10²¹ (fixed) | Blank et al. 2025 (P) | RSL, GNSS | best 10 of 35; depends on ice |
| British Isles, 3-D | 70–115 km | 0.56×10²¹, 0.5×10²¹ | 3×10²¹ (fixed) | Blank et al. 2025 (P) | RSL, GNSS | as above |
| Barents Sea | 65–75 km | 4×10²⁰, LVZ (1–4)×10¹⁹ at 130–200 km | 2×10²² | Steffen & Kaufmann 2005 (FT) | RSL | LVZ depth uncertain |
| Barents Sea | insensitive (71–120 km) | 0.2–2×10²¹ | insensitive | Auriac et al. 2016 (A) | RSL | UM only |
| Barents / Fennoscandia | none preferred | ≥3×10²⁰ / ≥5×10²⁰ | 2.6×10²¹ (fixed) | Rovira-Navarro et al. 2020 (FT) | GRACE | lower bounds |
| Greenland | 120 km (top of range tested) | 5×10²⁰ | 2×10²¹ (1–10×10²¹) | Lecavalier et al. 2014 (P) | RSL | UM good, LM weak |
| E Greenland | 120 km | 3×10²⁰ | 5×10²² | Simpson et al. 2009 via Lecavalier (S) | RSL, east coast | superseded |
| Greenland | 80 km | 5×10²⁰ | 10²² | Khan et al. 2016 (FT) | GNSS, RSL | moderate |
| SE Greenland | 40 km | 10¹⁹ (200 km layer) | 2×10²² | Khan et al. 2016 (FT) | GNSS (2 sites) | weak |
| Greenland, decadal | 200 km (60–240) | 6–11×10¹⁹ | not given | Adhikari et al. 2021; Lewright et al. 2026 (FT, P) | GNSS + LIA load | transient, not glacial |
| North America | 102 (85–120) km | 5.1 (3.5–7.5)×10²⁰ | 1.3 (0.8–2.8)×10²² | Lambeck et al. 2017 (A) | RSL, lake tilts, GPS | 95% limits |
| North America (VM5a) | 60 km + 40 km at 10²² | 5×10²⁰ | ~1.6–3×10²¹ | Peltier et al. 2015; Argus et al. 2014 (A, P) | RSL, GPS, rotation | tuned jointly with ICE-6G |
| global | 104 ± 29 km | 10^(20.76±0.24) | 10^21.57 or 10^22.38 | Caron et al. 2017 (P) | 5720 RSL | bimodal LM |
| far field | 40–70 km | 1.5 (1–2)×10²⁰ | 7 (1–20)×10²² or 2 (0.7–4)×10²¹ | Lambeck et al. 2014 (FT) | ~1000 far-field RSL | lithosphere poor, UM good, LM bimodal |
| W Antarctica, 3-D | 63 ± 13 km (LAB) | ~10¹⁹ at 150–350 km | not addressed | Hazzard et al. 2023 (P) | seismic (ANT-20) | ≤1 order |
| E Antarctica, 3-D | 233 ± 41 km (LAB) | >10²¹, up to 10²³ in the craton | not addressed | Hazzard et al. 2023; Ivins et al. 2021 (P) | seismic | |
| Amundsen Sea | not read | ~4×10¹⁸ | not read | Barletta et al. 2018 (A) | GNSS (decadal) | transient |
| Gowan (PaleoMIST) | 120 km | 4×10²⁰ | 4×10²² | Gowan et al. 2021 | (our reference) | |

## 5. Methods for exploring Earth uncertainty

### 5.1 Established approaches

- **1-D grid searches and Bayesian inversions with regional ice scaling** are the established way to estimate Earth parameters. Examples are Lambeck et al. (1998, 2014, 2017), Caron et al. (2017, 2018) and Lau et al. (2016).
  - Caron et al. (2017) solve regional ice scale factors by least squares for each sampled Earth. This makes the ice–Earth trade-off explicit: ±20–60% regional ice for a shift of the lower-mantle mode.
  - Caron et al. (2018, A) use 128,000 models with GPS and RSL. They find that polar GIA uncertainty exceeds that of 14-yr GRACE trends.
- **Scalar scaling of tomographic anomalies** is the established way to build 3-D Earths: log Δη ∝ β δv_s/v_s (Li et al. 2018), ε (Pan et al. 2022; Milne et al. 2018) or r (Bagge et al. 2021). It is explored by grid search or small factorial ensembles.
  - Li et al. (2018, P) find that the best RSL fit with ICE-6G has β_UM ≈ 0 and β_LM ≈ 0.6. RSL alone poorly constrains β_LM, no model fits Fennoscandian uplift and gravity rates together, and β trades off with the background viscosity.
- **History matching with emulators and explicit structural-error terms** is now in use. Lecavalier & Tarasov (2025, FT) run 9293 Antarctic simulations with a 1-D Earth (lithosphere, UM, LM). They add a structural error of +1/−0.5 mm/yr on uplift for neglecting 3-D structure (from Fennoscandian 1-D/3-D comparisons) and 5% on RSL.
- **Adjoint sensitivity kernels** for 3-D viscosity (Crawford et al. 2018; Lloyd et al. 2024) are maturing. The announced inversion (Part II) is not yet published.

### 5.2 Experimental or not found

- **PCA or EOF bases over expert 3-D viscosity fields:** none found. EOFs appear in GRACE-based GIA separation only for hydrology (Sun & Riva 2020, via a separate agent).
- **Transdimensional GIA inversion:** none found. The nearest is a trans-D inversion of the geoid for radial viscosity (Rudolph et al. 2015, A via a separate agent).
- **Joint Bayesian inversion of 3-D Earth and ice:** none found.

### 5.3 Pitfalls

Each pitfall is listed with its source.

1. **Ice–Earth trade-off.** The lower mantle is non-unique given ice uncertainty (Caron et al. 2017; Whitehouse 2018). Fits of 3-D Earths depend on whether the ice history was tuned on a 1-D Earth (Blank et al. 2025).
2. **Effective 1-D viscosity depends on the observable.** An ice-weighted average works better than a global one, and 1-D reproduces 3-D to no better than ~10% (Paulson et al. 2005). The bias of 1-D inversions concentrates in the upper mantle and transition zone (Lau et al. 2018, A). A local average under a site is not the right 1-D equivalent (Austermann et al. 2021).
3. **The amplitude of lateral variations is poorly known.** Halving ε halves the log range (Pan et al. 2022). The r range of Bagge spans a factor of 3.5. Calibrated anelastic conversions narrow it to ~0.2 log units in the upper mantle (Austermann et al. 2021; Hazzard et al. 2023), but these calibrations rest on assumed reference viscosities. Lower-mantle ranges of 6–7 orders are an upper bound (Pan et al. 2022).
4. **The lithosphere definition differs between fields.** Bagge uses a 10^23.5 Pa s contour, Lloyd an isotherm with no elastic lid, and Pan an elastic layer of mapped thickness scaled to a global mean. Nield et al. (2018, P) cite elastic thresholds of 10²²–10²⁵ Pa s and show that a viscously layered lithosphere acts as a thinner elastic one. They also find 1-D/3-D uplift differences of ±3 mm/yr in West Antarctica.
5. **Floors and caps.** Bagge floors at 10¹⁹ and Lloyd caps at 2×10¹⁹ Pa s for time-step reasons. Both clip the softest West Antarctic mantle (~10¹⁸–10¹⁹; Barletta et al. 2018; Hazzard et al. 2023).
6. **Lateral resolution.** The fields are resolved to degree 63 (Bagge) and degree 64 (Lloyd). Seakon runs reach ~12 km at the surface. Latychev et al. (2005, A) note that uplift above 1 cm/yr needs finer grids. *Inference:* VILMA2 at degree 170 is finer than the input structures. The limit is the structures, not the solver, except where lithosphere gradients are sharp (West Antarctica; Nield et al. 2018).
7. **Transition-zone viscosity.** It carries the largest 1-D bias (Lau et al. 2018). In Bagge's runs, s16 against sc06 changes ice-centre RSL as much as r does.
8. **Compressibility.** It roughly doubles horizontal rates and leaves vertical rates nearly unchanged (Tanaka et al. 2011). Our incompressible VILMA2 should therefore be scored on vertical GNSS only.
9. **Toroidal flow and plate boundaries.** Weak plate boundaries produce toroidal surface motion at GPS accuracy (Klemann et al. 2008, A via a separate agent), and 3-D viscosity perturbs horizontal velocities by ~1 mm/yr (Latychev et al. 2005). This again matters for horizontal rates only.
10. **3-D structure does not guarantee a better fit.** No improvement appears in Greenland GNSS (Milne et al. 2018), the US Pacific coast (Yousefi et al. 2021), eastern North America (Parang et al. 2024) or Fennoscandia with ICE-6G (Blank et al. 2025).

## 6. GNSS as a GIA constraint

**Depth sensitivity.**

- Present-day uplift senses the viscosity beneath the site, beneath the former load and in the corridor between them, and is insensitive elsewhere (Paulson et al. 2005).
- In Fennoscandia, BIFROST resolves ~200 km near the base of the upper mantle and ~700 km at the top of the lower mantle, with little constraint on the deeper half of the mantle (Milne et al. 2004).
- Hudson Bay RSL constrains mainly the top half of the lower mantle (Mitrovica & Peltier, as summarized by a separate agent; not read).
- A thin low-viscosity layer below the lithosphere barely changes Holocene RSL but strongly raises present-day uplift (Pan et al. 2024).

*Inference:* GNSS and RSL are complementary for the upper mantle. However, GNSS also responds to shallow, short-timescale structure that a glacial-cycle Maxwell model does not represent.

**Elastic and recent-load corrections.**

- **Greenland:** uplift is dominated by the elastic response to present-day ice loss (Bevis et al. 2012, A). Adhikari et al. (2021) find an elastic-correction uncertainty (1σ) of up to ~1.2 mm/yr, dominated by the load, not by the Earth model. After that correction, standard GIA leaves a residual of 3.4 ± 1.9 mm/yr, which Little Ice Age loading with a transient viscosity of 6–11×10¹⁹ Pa s explains.
- **Antarctica:** the Maxwell time at 10¹⁹ Pa s is ~5 yr. At ≤10¹⁸ Pa s the viscous response to recent change exceeds the elastic one by about ten times (Whitehouse et al. 2019, P). The GIANT-REGAIN data set (Buchta et al. 2025: 286 sites, IGb14) provides no elastic correction. Lecavalier & Tarasov (2025) call the uncertainty of their elastic correction "ill-defined".

**Error assumptions.** We found no source for a generic 0.5–1 mm/yr floor. The documented components are:

- the reference frame (CM to CE and origin rate), ~0.2 mm/yr (Schumacher et al. 2018; Buchta et al. 2025);
- NKG2016LU, ~0.2 mm/yr where data are dense and ~0.6 mm/yr including frame stability (Lantmäteriet memo citing Vestøl et al. 2019, via a separate agent);
- the 1-D/3-D structural term, +1/−0.5 mm/yr (Lecavalier & Tarasov 2025);
- the elastic correction in Greenland, ≤1.2 mm/yr (Adhikari et al. 2021).

Schumacher et al. (2018) report ICE-6G misfits (RMS) of 1.67 mm/yr globally, 3.70 in Greenland and 3.84 in Antarctica, with data errors that exclude the elastic and frame terms. Peak observed uplift in Fennoscandia is ~10.3 mm/yr (Kierulf et al. 2021; Vestøl et al. 2019).

## 7. Implications for kryos

### 7.1 Fields to include in the Earth basis

- **Bagge, all 8 FesmData structures.** They span r (0.2857–1.0) and the three radial-profile classes. *Inference:* within this set, the principal components will mostly be "amplitude" (r) and "transition-zone profile" (s16/sc06/sc06b) of one tomography.
- **Lloyd et al. (2024), GLAD-M25.** A different tomography and a calibrated anelastic conversion, with a lower mantle that varies laterally.
- **Pan et al. (2022), as held in FesmData.** A third, unidentified tomography with an exponential conversion. It also has an independent elastic-lithosphere map (AF or YO), the only one of the three with a lithosphere not derived from its own viscosity field. First identify its background and ε from its spherical means (Sect. 1.6), and ask the authors which realization it is.
- **Not to include as a field:** decadal GNSS-based viscosities (Amundsen Sea, Greenland Little Ice Age). They describe transient rheology, not the glacial-cycle Maxwell response.

### 7.2 Construction of the basis

*Inference*, following the pitfalls in Sect. 5.3:

1. **Bring the lithosphere to one definition first.** Remove each field's own lithosphere (Bagge's 10^23.5 contour, Lloyd's 1.5×10²⁶ Pa s layer, Pan's map). Treat the lid as a separate parameter, either a viscosity threshold or a thickness map scaled to a global mean as Pan et al. do. Otherwise the first principal component is the lithosphere definition.
2. **Split each field into a radial mean profile and a lateral anomaly in log10 η.** Build the principal components from the anomalies only. Keep the radial profile as explicit global shifts (UM, LM, optionally the transition zone), because Pan et al. find the background profile first-order and it has direct literature ranges.
3. **Separate amplitude from geometry.** Normalize each anomaly field and keep an amplitude factor as a parameter, the analogue of ε, r and β. Pan et al. and Bagge et al. show that amplitude mostly rescales the pattern.
4. **Weight the three sources equally.** Eight Bagge fields against one each from Lloyd and Pan would otherwise make the basis mostly Bagge.
5. **Treat lateral variation in the lower mantle separately.** Bagge has none below 870 km. Pan's 6–7 orders at 1071 km are, by the authors' account, an upper bound. Lateral lower-mantle structure barely affects Fennoscandian GNSS (Steffen et al. 2007), but it does affect far-field RSL (Lloyd et al. 2024) and Hudson Bay uplift (Li et al. 2020). A separate amplitude factor for the lower mantle, possibly zero, keeps this explicit.
6. **Keep the 10¹⁹ Pa s floor.** It sits at the steady-state West Antarctic value (Hazzard et al. 2023) and limits the time step.

### 7.3 Parameter ranges

*Inference* from Sect. 4.5. These are ranges for the spherical means and shifts, not hard limits.

| parameter | suggested range | basis |
|---|---|---|
| lid, global mean thickness | 60–120 km | 1-D regional preferences: 40–70 km (far field), 65–160 km (Fennoscandia), 85–120 km (North America), 120 km (Greenland); Pan's means of 71 and 96 km |
| lid, as viscosity threshold | 10²²–10²⁵ Pa s | elastic thresholds cited by Nield et al. (2018); Bagge's 10^23.5 Pa s inside the range |
| UM mean | log10 η = 20.2–20.9 | 1.5×10²⁰ (far field) to 7×10²⁰ (Fennoscandia, GNSS); Caron et al. 2017: 20.76 ± 0.24 |
| LM mean | log10 η = 21.3–22.8 | bimodal: ~2×10²¹ (Greenland, VM5a, Lambeck low) to ~7×10²² (Lambeck high); Lambeck et al. 2017: 1.3×10²² Pa s; Gowan's 4×10²² is at the upper end |
| transition-zone shift | ±0.5 in log10 η, or Bagge's three profile classes | Lau et al. 2018; Bagge s16/sc06/sc06b |
| anomaly amplitude | 0.5–1.5 times the reference amplitude (at least a factor of 2–3 in total) | ε 0.02–0.04 (Pan); r 0.2857–1.0 (Bagge) |
| regional UM shifts | ±0.5 in log10 η | Fennoscandia 1-D 3.4–7×10²⁰ against 3-D 1.3–4.9×10²¹; Barents about half of Fennoscandia (Rovira-Navarro et al. 2020) |

### 7.4 Regions where a better Earth should help

- **Fennoscandia.**
  - *Inference:* moving away from structure 6 (lower r, or sc06b) should bring the Gulf of Bothnia back towards the 1-D result, since Bagge et al. show the 3-D craton reduces uplift there.
  - As long as the ice is PaleoMIST, built on Gowan's 1-D Earth, the 1-D Earth keeps an advantage (Blank et al. 2025). A fair test of 3-D structure needs ice adjusted on the same Earth, which the coupled ensemble provides.
- **Barents Sea.** A thinner lid and a softer upper mantle than Fennoscandia. Splitting the Fennoscandia/Barents score would show whether the regional shift is pulled by the Barents data.
- **Hudson Bay.** Our science Earth already fits better than Gowan's ([vilma2-paleomist.md](vilma2-paleomist.md), finding 9).
- **US Atlantic margin and Caribbean.** Pan et al. find that stiff North American lithosphere and mantle suppress the peripheral bulge by up to a few tens of metres at the LGM. Thompson et al. (2023) find 3-D clearly better than 1-D there for MIS 5a.
- **East versus west Greenland.** Seismic and GNSS studies support a softer east (Iceland hotspot track) and a stiffer west and centre. *Inference:* the effect on our scores is likely modest, since most Greenland misfit is in the west.
- **Antarctica.** West–east contrasts are large, but GAPSLIP has few points (47). The main effect is on the coupled ice dynamics (Albrecht et al. 2024), not on the score.

### 7.5 What the Earth cannot fix

- **West Greenland** (Sisimiut, Disko Bay, model 40–100 m too low). 3-D structure reduces the RSL fall at Disko Bugt (Milne et al. 2018), the wrong direction. The literature points to Greenland ice extent and thickness and to the North American ice load (Lecavalier et al. 2014).
- **St. Lawrence** (+16 to +56 m). ICE-6G and ICE-7G show the same over-prediction of up to ~60 m. It is reduced by thinning the St. Lawrence ice stream (Li et al. 2022), and 3-D Earths do not improve eastern North America (Parang et al. 2024).
- **Far field** (+10 to +20 m). The 3-D effect is a few metres of either sign (Pan et al. 2022) and has the wrong sign at Barbados. The lower-mantle trade-off also goes the wrong way. The remedy is more LGM ice, plus the VILMA2 ocean bias against SELEN (vilma2-paleomist finding 2).

### 7.6 GNSS

*Inference*, from Sect. 6:

1. **Score vertical rates only.** VILMA2 is incompressible, which affects horizontal rates (Tanaka et al. 2011).
2. **Start with Fennoscandia and Hudson Bay.** There the data are dense and the elastic signal small. Use a model–data error of about 1 mm/yr: data, frame and the +1/−0.5 mm/yr structural term (Lecavalier & Tarasov 2025).
3. **Treat Greenland and Antarctica separately.**
   - Use them only after an elastic correction for present-day ice loss, with that correction's uncertainty in the error (~1 mm/yr or more).
   - Without a Little Ice Age load and transient rheology, expect a systematic positive residual of a few mm/yr in Greenland (Adhikari et al. 2021) and in the Amundsen Sea.
   - Either exclude these regions from the Earth score or give them errors of ~2–3 mm/yr, comparable to the ICE-6G RMS misfits there (Schumacher et al. 2018).
4. **Use the 1-D literature as a check.** NKG2016LU's Earth (160 km, 7×10²⁰ Pa s) is an independent check on what our Fennoscandian shift should approach.

## References

Read level in brackets.

- Adhikari, S., et al. (2021). Decadal to centennial timescale mantle viscosity inferred from modern crustal uplift rates in Greenland. GRL 48, e2021GL094040. doi:10.1029/2021GL094040 [FT]
- Albrecht, T., Bagge, M., Klemann, V. (2024). Feedback mechanisms controlling Antarctic glacial-cycle dynamics simulated with a coupled ice sheet–solid Earth model. The Cryosphere 18, 4233–4255. doi:10.5194/tc-18-4233-2024 [P]
- Argus, D. F., Peltier, W. R., Drummond, R., Moore, A. W. (2014). The Antarctica component of postglacial rebound model ICE-6G_C (VM5a) based on GPS positioning, exposure age dating of ice thicknesses, and relative sea level histories. GJI 198, 537–563. doi:10.1093/gji/ggu140 [P]
- Auriac, A., et al. (2016). Glacial isostatic adjustment associated with the Barents Sea ice sheet: a modelling inter-comparison. QSR 147, 122–135. doi:10.1016/j.quascirev.2016.02.011 [A]
- Austermann, J., Mitrovica, J. X., Latychev, K., Milne, G. A. (2013). Barbados-based estimate of ice volume at Last Glacial Maximum affected by subducted plate. Nat. Geosci. 6, 553–557. doi:10.1038/ngeo1859 [A]
- Austermann, J., Hoggard, M. J., Latychev, K., Richards, F. D., Mitrovica, J. X. (2021). The effect of lateral variations in Earth structure on Last Interglacial sea level. GJI 227, 1938–1960. doi:10.1093/gji/ggab289 [P]
- Bagge, M., Klemann, V., Steinberger, B., Latinović, M., Thomas, M. (2021). Glacial-isostatic adjustment models using geodynamically constrained 3D Earth structures. G-Cubed 22, e2021GC009853. doi:10.1029/2021GC009853 [FT]
- Barletta, V. R., et al. (2018). Observed rapid bedrock uplift in Amundsen Sea Embayment promotes ice-sheet stability. Science 360, 1335–1339. doi:10.1126/science.aao1447 [A]
- Bevis, M., et al. (2012). Bedrock displacements in Greenland manifest ice mass variations, climate cycles and climate change. PNAS 109, 11944–11948. doi:10.1073/pnas.1204664109 [A]
- Blank, B., Bradley, S. L., van der Wal, W. (2025). Uplift and sea level constraints on 3-D upper mantle viscosity in Northern Europe. GJI 244, ggaf405. doi:10.1093/gji/ggaf405 [P]
- Buchta, E., et al. (2025). Advancing geodynamic research in Antarctica: reprocessing GNSS data to infer consistent coordinate time series (GIANT-REGAIN). ESSD 17, 1761–1780. doi:10.5194/essd-17-1761-2025 [FT]
- Caron, L., Métivier, L., Greff-Lefftz, M., Fleitout, L., Rouby, H. (2017). Inverting glacial isostatic adjustment signal using Bayesian framework and two linearly relaxing rheologies. GJI 209, 1126–1147. doi:10.1093/gji/ggx083 [P]
- Caron, L., Ivins, E. R., Larour, E., Adhikari, S., Nilsson, J., Blewitt, G. (2018). GIA model statistics for GRACE hydrology, cryosphere, and ocean science. GRL 45, 2203–2212. doi:10.1002/2017GL076644 [A]
- Crawford, O., et al. (2018). Quantifying the sensitivity of post-glacial sea level change to laterally varying viscosity. GJI 214, 1324–1363. doi:10.1093/gji/ggy184 [A]
- Gowan, E. J., et al. (2021). A new global ice sheet reconstruction for the past 80 000 years. Nat. Commun. 12, 1199. doi:10.1038/s41467-021-21469-w
- Hazzard, J. A. N., Richards, F. D., Goes, S. D. B., Roberts, G. G. (2023). Probabilistic assessment of Antarctic thermomechanical structure: impacts on ice sheet stability. JGR Solid Earth 128, e2023JB026653. doi:10.1029/2023JB026653 [P]
- Ivins, E. R., van der Wal, W., Wiens, D. A., Lloyd, A. J., Caron, L. (2021). Antarctic upper mantle rheology. Geol. Soc. London Memoirs 56, 267–294. doi:10.1144/M56-2020-19 [P]
- Khan, S. A., et al. (2016). Geodetic measurements reveal similarities between post–Last Glacial Maximum and present-day mass loss from the Greenland ice sheet. Sci. Adv. 2, e1600931. doi:10.1126/sciadv.1600931 [FT]
- Kierulf, H. P., et al. (2014). A GPS velocity field for Fennoscandia and a consistent comparison to glacial isostatic adjustment models. JGR Solid Earth 119, 6613–6629. doi:10.1002/2013JB010889 [S]
- Kierulf, H. P., et al. (2021). A GNSS velocity field for geophysical applications in Fennoscandia. J. Geodyn. 146, 101845. doi:10.1016/j.jog.2021.101845 [snippet]
- Lambeck, K., Smither, C., Johnston, P. (1998). Sea-level change, glacial rebound and mantle viscosity for northern Europe. GJI 134, 102–144. doi:10.1046/j.1365-246x.1998.00541.x [A]
- Lambeck, K., Rouby, H., Purcell, A., Sun, Y., Sambridge, M. (2014). Sea level and global ice volumes from the Last Glacial Maximum to the Holocene. PNAS 111, 15296–15303. doi:10.1073/pnas.1411762111 [FT]
- Lambeck, K., Purcell, A., Zhao, S. (2017). The North American Late Wisconsin ice sheet and mantle viscosity from glacial rebound analyses. QSR 158, 172–210. doi:10.1016/j.quascirev.2016.11.033 [A, EGU2017-5721]
- Latychev, K., et al. (2005). Glacial isostatic adjustment on 3-D Earth models: a finite-volume formulation. GJI 161, 421–444. doi:10.1111/j.1365-246X.2005.02536.x [A]
- Lau, H. C. P., et al. (2016). Inferences of mantle viscosity based on ice age data sets: radial structure. JGR Solid Earth 121, 6991–7012. doi:10.1002/2016JB013043 [A]
- Lecavalier, B. S., et al. (2014). A model of Greenland ice sheet deglaciation constrained by observations of relative sea level and ice extent. QSR 102, 54–84. doi:10.1016/j.quascirev.2014.07.018 [P]
- Lecavalier, B. S., Tarasov, L. (2025). A history-matching analysis of the Antarctic Ice Sheet since the last interglacial – Part 2: Glacial isostatic adjustment. The Cryosphere 19, 6673–6689. doi:10.5194/tc-19-6673-2025 [FT]
- Lewright, L., et al. (2026). Projections of 21st-century sea-level fall along coastal Greenland. Nat. Commun. 17, 353. doi:10.1038/s41467-025-68182-6 [P]
- Li, T., Wu, P., Steffen, H., Wang, H. (2018). In search of laterally heterogeneous viscosity models of glacial isostatic adjustment with the ICE-6G_C global ice history model. GJI 214, 1191–1205. doi:10.1093/gji/ggy181 [P]
- Li, T., et al. (2020). Uncertainties of glacial isostatic adjustment model predictions in North America associated with 3D structure. GRL 47, e2020GL087944. doi:10.1029/2020GL087944 [A]
- Li, T., Peltier, W. R., et al. (2022). Resolving the influence of ice stream instability on postglacial RSL histories: the St Lawrence river channel ice stream. EGU General Assembly, EGU22-1343 [A and slides]
- Lidberg, M., Johansson, J. M., Scherneck, H.-G., Milne, G. A. (2010). Recent results based on continuous GPS observations of the GIA process in Fennoscandia from BIFROST. J. Geodyn. 50, 8–18. doi:10.1016/j.jog.2009.11.010 [A]
- Lloyd, A. J., Crawford, O., Al-Attar, D., Austermann, J., Hoggard, M. J., Richards, F. D., Syvret, F. (2024). GIA imaging of 3-D mantle viscosity based on palaeo sea level observations – Part I: Sensitivity kernels for an Earth with laterally varying viscosity. GJI 236, 1139–1171. doi:10.1093/gji/ggad455 [P]
- Milne, G. A., et al. (2004). Continuous GPS measurements of postglacial adjustment in Fennoscandia: 2. Modeling results. JGR 109, B02412. doi:10.1029/2003JB002619 [A]
- Milne, G. A., Latychev, K., Schaeffer, A., Crowley, J. W., Lecavalier, B. S., Audette, A. (2018). The influence of lateral Earth structure on glacial isostatic adjustment in Greenland. GJI 214, 1252–1266. doi:10.1093/gji/ggy189 [P]
- Mitrovica, J. X., Forte, A. M. (2004). A new inference of mantle viscosity based upon joint inversion of convection and glacial isostatic adjustment data. EPSL 225, 177–189. doi:10.1016/j.epsl.2004.06.005 [A]
- Mordret, A. (2018). Uncovering the Iceland hot spot track beneath Greenland. JGR Solid Earth 123, 4922–4941. doi:10.1029/2017JB015104 [A]
- Nield, G. A., et al. (2018). The impact of lateral variations in lithospheric thickness on glacial isostatic adjustment in West Antarctica. GJI 214, 811–824. doi:10.1093/gji/ggy158 [P]
- Pan, L., Milne, G. A., Latychev, K., Goldberg, S. L., Austermann, J., Hoggard, M. J., Mitrovica, J. X. (2022). The influence of lateral Earth structure on inferences of global ice volume during the Last Glacial Maximum. QSR 290, 107644. doi:10.1016/j.quascirev.2022.107644 [FT]
- Pan, L., Mitrovica, J. X., Milne, G. A., Hoggard, M. J., Woodroffe, S. A. (2024). Timescales of glacial isostatic adjustment in Greenland: is transient rheology required? GJI 237, 989–995. doi:10.1093/gji/ggae095 [P]
- Parang, S., et al. (2024). Constraining models of glacial isostatic adjustment in eastern North America. QSR 334, 108708. doi:10.1016/j.quascirev.2024.108708 [A]
- Paulson, A., Zhong, S., Wahr, J. (2005). Modelling post-glacial rebound with lateral viscosity variations. GJI 163, 357–371. doi:10.1111/j.1365-246X.2005.02645.x [FT]
- Paulson, A., Zhong, S., Wahr, J. (2007). Inference of mantle viscosity from GRACE and relative sea level data. GJI 171, 497–508. doi:10.1111/j.1365-246X.2007.03556.x [A]
- Peltier, W. R., Argus, D. F., Drummond, R. (2015). Space geodesy constrains ice age terminal deglaciation: the global ICE-6G_C (VM5a) model. JGR Solid Earth 120, 450–487. doi:10.1002/2014JB011176 [A]
- Richards, F. D., Hoggard, M. J., White, N., Ghelichkhan, S. (2020). Quantifying the relationship between short-wavelength dynamic topography and thermomechanical structure of the upper mantle using calibrated parameterization of anelasticity. JGR Solid Earth 125, e2019JB019062. doi:10.1029/2019JB019062 [P]
- Rovira-Navarro, M., et al. (2020). GRACE constraints on Earth rheology of the Barents Sea and Fennoscandia. Solid Earth 11, 379–395. doi:10.5194/se-11-379-2020 [FT]
- Schachtschneider, R., et al. (2022). An approach for constraining mantle viscosities through assimilation of palaeo sea level data into a glacial isostatic adjustment model. NPG 29, 53–75. doi:10.5194/npg-29-53-2022 [P]
- Schumacher, M., et al. (2018). A new global GPS data set for testing and improving modelled GIA uplift rates. GJI 214, 2164–2176. doi:10.1093/gji/ggy235 [FT]
- Steffen, H., Kaufmann, G. (2005). Glacial isostatic adjustment of Scandinavia and northwestern Europe and the radial viscosity structure of the Earth's mantle. GJI 163, 801–812. doi:10.1111/j.1365-246X.2005.02740.x [FT]
- Steffen, H., Wu, P. (2011). Glacial isostatic adjustment in Fennoscandia – a review of data and modeling. J. Geodyn. 52, 169–204. doi:10.1016/j.jog.2011.03.002 [A]
- Swierczek-Jereczek, J., et al. (2024). FastIsostasy v1.0 – a regional, accelerated 2D glacial isostatic adjustment (GIA) model accounting for the lateral variability of the solid Earth. GMD 17, 5263–5290. doi:10.5194/gmd-17-5263-2024 [P]
- Tanaka, Y., Klemann, V., Martinec, Z., Riva, R. E. M. (2011). Spectral-finite element approach to viscoelastic relaxation in a spherical compressible Earth: application to GIA modelling. GJI 184, 220–234. doi:10.1111/j.1365-246X.2010.04854.x [A, search summary]
- Tarasov, L., Dyke, A. S., Neal, R. M., Peltier, W. R. (2012). A data-calibrated distribution of deglacial chronologies for the North American ice complex from glaciological modeling. EPSL 315–316, 30–40. doi:10.1016/j.epsl.2011.09.010 [A]
- Thompson, S. B., Creveling, J. R., Latychev, K., Mitrovica, J. X. (2023). Three-dimensional glacial isostatic adjustment modeling reconciles conflicting geographic trends in North American marine isotope stage 5a relative sea level observations. Geology 51, 808–812. doi:10.1130/G51257.1 [FT]
- van der Wal, W., et al. (2013). Glacial isostatic adjustment model with composite 3-D Earth rheology for Fennoscandia. GJI 194, 61–77. doi:10.1093/gji/ggt099 [A]
- Vestøl, O., Ågren, J., Steffen, H., Kierulf, H., Tarasov, L. (2019). NKG2016LU: a new land uplift model for Fennoscandia and the Baltic Region. J. Geod. 93, 1759–1779. doi:10.1007/s00190-019-01280-8 [S]
- Whitehouse, P. L. (2018). Glacial isostatic adjustment modelling: historical perspectives, recent advances, and future directions. Earth Surf. Dynam. 6, 401–429. doi:10.5194/esurf-6-401-2018 [FT]
- Whitehouse, P. L., Gomez, N., King, M. A., Wiens, D. A. (2019). Solid Earth change and the evolution of the Antarctic Ice Sheet. Nat. Commun. 10, 503. doi:10.1038/s41467-018-08068-y [P]
- Wiens, D. A., Shen, W., Lloyd, A. J. (2021). The seismic structure of the Antarctic upper mantle. Geol. Soc. London Memoirs 56, 195–212. doi:10.1144/M56-2020-18 [A]
- Yokoyama, Y., et al. (2018). Rapid glaciation and a two-step sea level plunge into the Last Glacial Maximum. Nature 559, 603–607. doi:10.1038/s41586-018-0335-4 [A]
- Yousefi, M., Milne, G. A., Latychev, K. (2021). Glacial isostatic adjustment of the Pacific Coast of North America: the influence of lateral Earth structure. GJI 226, 91–113. doi:10.1093/gji/ggab053 [A]
- Zhao, S., Lambeck, K., Lidberg, M. (2012). Lithosphere thickness and mantle viscosity inverted from GPS-derived deformation rates in Fennoscandia. GJI 190, 278–292. doi:10.1111/j.1365-246X.2012.05454.x [S]

**Not read, mentioned only through other papers:** Berg et al. (2024), Greenland GNSS uplift, via Lewright et al. (2026). Kuchar et al. (2019) and Lambeck et al. (2010) were not retrievable. Klemann, Martinec & Ivins (2008), Rudolph et al. (2015), Lau et al. (2018), Steffen et al. (2006, 2007), Mitrovica & Peltier (1991–1993), Sun & Riva (2020) and Weerdesteijn & Conrad (2024) are reported from search summaries by a separate agent and not checked here.
