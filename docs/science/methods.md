---
title: "Methods"
---

The model, the kinds of simulation and how they are named, and the ladder that makes each change of score
attributable. The data and the scores are in [Data and scores](data-scores.md).

## The model

**Solid Earth and sea level: VILMA2.** A global, viscoelastic (Maxwell) Earth with self-gravitating sea level,
rotational feedback and the sea-level equation solved on a moving ocean function, at spherical-harmonic degree 170. The
viscosity is either 1-D (radial layers) or 3-D (a field of log10 η on a longitude-latitude-radius grid); above a
threshold it is elastic, which sets the lithosphere ([VILMA2 backend](https://fesmc.github.io/kryos/backend-vilma2.html)).

**Ice: PlasticIceSheet.** A perfectly-plastic ice sheet: the ice surface follows from the basal shear stress τ_b and
the margin, on the bed and the sea level of the time. The margin follows a chronology (PaleoMIST's), τ_b follows Gowan's
domains times factors that we calibrate or sample. Northern Hemisphere on 32 km, Antarctica on 16 km, Patagonia on
8 km ([PlasticIceSheet backend](https://fesmc.github.io/kryos/backend-plastic.html)).

**Coupling.** In kryos the regional ice domains load one global VILMA2 Earth, which returns the bed and the sea level to
each domain. Outside the domains the Earth sees an ice history read from file
([Global hub and regional domains](https://fesmc.github.io/kryos/global-hub.html)).

## Two kinds of simulation

- **Imposed ice.** A published ice history forced onto VILMA2. The ice does not respond to the Earth. This isolates the
  solver and the Earth, and costs minutes (1-D Earth) to hours (3-D Earth).
- **Online ice.** PlasticIceSheet computes the ice every step from τ_b and the margin, on the bed that VILMA2 deforms.
  This is our reconstruction method. A coupled run on a 3-D Earth costs about 240 core-h.

## The present-day correction r

A glacial-cycle run must end with today's bed. The run starts from a reference topography, today's plus a correction r,
and r is iterated until the run's 0 ka bed matches today's to 1.5 m, in 3-5 passes. So r belongs to one ice history and
one Earth: a score then measures the ice and the Earth, not an offset of the topography. A run may borrow the r of a
close relative (for example a 3-D Earth the r of its radial profile) and run once
([the reference topography](https://fesmc.github.io/kryos/global-hub.html#the-reference-topography-refcorr)).

## Earths

| Earth | what it is |
|---|---|
| Gowan's 1-D | PaleoMIST's Earth: 120 km elastic lithosphere, upper mantle 4×10²⁰ Pa s, lower mantle 4×10²² Pa s |
| Bagge's structures | 8 of the 3-D viscosity structures of Bagge et al. (2021) (below), and their radial means as 1-D Earths |
| Lloyd, Pan | the 3-D fields of Lloyd et al. (2024) and Pan et al. (2022); used through the Earth basis |
| sampled Earths | a radial profile plus the experts' lateral pattern, scaled, plus regional upper-mantle shifts ([The Earth search](earth-search.md#the-parameters)) |

The Bagge structures used, by their lateral reduction factor and radial profile:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| Bagge name | v_0.2857_s16 | v_0.3_s16 | v_0.4_s16 | v_0.5_s16 | v_0.75_s16 | v_1.0_s16 | v_0.4_sc06 | v_0.4_sc06b |

Bagge's reference is structure 3. Structure 6 has the largest lateral contrast; ens-v1 and ens-v2 ran on it.

## Naming

Every component is named `<thing>-<source>-<variant>`, the variant left out when there is only one. A simulation is its
ice and its Earth, written `<ice> | <earth>`.

| thing | what it is |
|---|---|
| `pub` | a published reconstruction as a whole: its ice, its Earth, its own GIA solver and the RSL it published |
| `ice` | an ice history imposed on VILMA2 |
| `tau` | τ_b for PlasticIceSheet online |
| `E1D`, `E3D` | a 1-D or a 3-D Earth |

| source | what it is |
|---|---|
| `PMIST` | PaleoMIST 1.0 (Gowan et al. 2021), the minimal MIS 3 scenario: its ice, its τ_b through time, its 1-D Earth |
| `GLAC1D` | GLAC-1D (Tarasov et al. 2012, Briggs et al. 2014) from 80 ka, North America nn9894 |
| `B21` | Bagge et al. (2021) |
| `ENS1`, `ENS2` | our weighted τ_b ensembles ens-v1 and ens-v2 |
| `ES1`, `ES2` | the Earth search: wave A (1-D Earths), wave B and B2 (3-D Earths; `E1D-ES2-E<nnn>` is the radial profile of `E3D-ES2-E<nnn>`) |

| variant | of | meaning |
|---|---|---|
| `9927` | `ice-GLAC1D` | North America nn9927 instead of nn9894 |
| `PD` | `tau-PMIST` | times the present-day factors (fitted to today's Greenland and Antarctica) |
| `E<nnn>` | `tau-ENS2`, `E3D-B21-6`, `E1D-ES1`, `E3D-ES2` | member or Earth nnn |
| `1` ... `8` | `E3D-B21`, `E1D-B21` | Bagge structure k; `E1D-B21-k` is its radial mean |

**Defaults.** The GIA solver is VILMA2 except for `pub`: toroidal flow on, ice density 910 kg m⁻³, degree 170, and the r
of the same ice and Earth. `@NH16` marks the Northern Hemisphere on 16 km.

Examples: `ice-PMIST | E1D-PMIST` is PaleoMIST's ice and Earth on our solver; `tau-PMIST-PD | E3D-B21-6` is our ice from
PaleoMIST's τ_b with the present-day factors, on Bagge's structure 6.

## The ladder

Each simulation changes one thing against its parent, so the change of score is that thing's:

| step | changes | isolates |
|---|---|---|
| **P** | the published reconstructions | the bar |
| **I** | the published ice on VILMA2, then the Earth varied | the solver, then the Earth |
| **T** | our ice online from PaleoMIST's τ_b | the ice model |
| **D** | with the present-day factors | today's ice |
| **E** | ensembles of τ_b (and the Earth) | the ice history |
| **S** | the Earth search: imposed ice on sampled Earths | the Earth, independently of our ice |

The scores of every step: [Results](results.md).
