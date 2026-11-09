# Margin targets from original sources (future option)

Status (2026-11-09): not pursued for now. Alex decided to keep the current margins (IceHistory-PaleoMIST `f_ice > 0.5`,
the backend default); this note records how to make them from the original polygons if we revisit it.


Goal: the margin chronology of the plastic backend from the margin polygons themselves, not from the
PaleoMIST 1.0 output (`f_ice > 0.5` of IceHistory-PaleoMIST). Evidence: ICESHEET 2.0 at 641d172
(`$FESMDATA_SOURCES/gowan2021_icesheet/icesheet-641d172.../`, local copy `~/models/icesheet`), Gowan et al.
(2021, Nat. Commun. 12, 1199, Methods "Margin reconstructions", read from Europe PMC PMC7902671), the
data pages of the other sources. Prototype on Levante: `/work/ba1442/robinson/margin-proto/`.

## 1. How Gowan et al. (2021) built the margins

- Hand-drawn grounded-ice extent polygons, lon-lat WGS84 (GMT/OGR text exported from QGIS, `@Je4326`),
  one file per region and time, every 2.5 kyr from 80 to 0 ka: `global/{North_America,Eurasia,Antarctica,Patagonia}/margins/<t>.gmt`
  (North_America includes Greenland, the Innuitian region and Iceland). 2.5 kyr is the paper's interval;
  the MIS 3-4 timing is uncertain by "several thousand years" (Methods).
- After the LGM, existing reconstructions (Methods refs 13, 36, 38, 59, 63-67, 82):
  - North America: Dyke (2004) as in NAICE (Gowan et al. 2016) and Khosravi (2017, Innuitian), replaced in
    the east by Shaw et al. (2006, Atlantic Canada, 25-15 ka) and Occhietti et al. (2011, Quebec-Labrador,
    15.5-7.5 ka), Cordillera at 25 ka from Clague and Ward (2011) (`North_America/margins/readme` 143-229).
    Not Dalton et al. (2020): only the Hudson Bay MIS 3 compilation of Dalton et al. (2019) is used.
  - Greenland: LGM from Funder et al. (2011), present from the RTopo-2 grounding line (Schaffer et al. 2016),
    Holocene and MIS 3-4 as inward/outward buffers of these (-5, -10, -25 km; readme 103-123, 157-163).
  - Eurasia: DATED-1 most-credible (Hughes et al. 2016) for 25-10 ka, the 2.5 kyr slices in between
    "intermediate" of the 1 kyr DATED slices; 30-40 ka from the DATED-1 38-35 ka min/max and 31/29 ka;
    MIS 3-4 from Helmens (2014) and Ingolfsson and Landvik (2013); 0-7.5 ka modern (GLIMS) (`Eurasia/margins/readme`).
  - Antarctica: RAISED Consortium (2014), scenario B, at 20, 15, 10, 5 ka; 0 ka = Bedmap2 grounding line
    (smoothed by a -5/+5 km buffer); slices in between with Gowan's contour-interpolation tool
    (`Antarctica/margins/readme` 1-25, 67-71).
  - Patagonia: Ehlers et al. (2011) LGM, contour interpolation to the present (`Patagonia/margins/readme` 7).
- Before the LGM (NA, Antarctica): "drawn from scratch" from stratigraphic evidence, tied to Heinrich events
  (local maxima at H events; MIS 4 maximum at 60 ka = H6; retreats at 50, 42.5, 37.5 ka).
- Scenarios: there is no `margins_max`. North America has a second MIS 3 scenario, `margins/a1/`
  (37.5-55 ka, 8 files: Hudson Bay stays ice-covered, the "maximal" one; readme 241). `run.sh` takes
  `margins/<folder>/<t>.gmt` when it exists, else `margins/<t>.gmt` (`prepare_icesheet.sh` 172-177). The
  minimal scenario is PaleoMIST 1.0 (= IceHistory-PaleoMIST), the maximal one is PaleoMIST_max.
- Into ICESHEET (`global/prepare_icesheet.sh`): `gmt mapproject` to the regional projection (NA: Lambert
  azimuthal equal-area at 60N 94W, 5 km; 193-197), split into single polygons, each one an ICESHEET run along
  flowlines from its margin (204-262), then `gmt grdmask` of the projected polygons on the 5 km grid (294;
  boundary nodes count as inside) multiplies the gridded thickness (301). So PaleoMIST's ice extent is the
  polygon set rasterised at 5 km, and IceHistory-PaleoMIST `f_ice` is that, sampled 5x5 per 0.25 deg cell and
  remapped conservatively to the kryos grid: two resamplings of the polygons.

## 2. Original sources

| region | source | span | resolution | format | FesmData | licence |
|---|---|---|---|---|---|---|
| NA+Greenland+Iceland, Eurasia, Antarctica, Patagonia | Gowan ICESHEET 2.0 `margins/` (+ `a1/` for NA MIS 3) | 80-0 ka | 2.5 kyr, one scenario (+a1) | polygons, lon-lat | yes: source `gowan2021_icesheet` (used by PaleoMIST_taub only) | GPL-3.0 |
| Laurentide, Cordilleran, Innuitian (no Greenland) | NADI-1, Dalton et al. (2023), Zenodo 10.5281/zenodo.8161764 | 25-1 ka (cal) | 0.5 kyr (49 steps), min/optimal/max | shapefiles | no | CC BY 4.0 |
| same | Dalton et al. (2020), GHub resource 4878 | 18-1 14C ka (36 isochrones) | ~0.5 kyr, radiocarbon ages | shapefiles | no | superseded by NADI-1 |
| same | Dyke (2004) / GSC Open File 1574 | 18-1 14C ka | ~0.5-1 kyr | shapefiles | no | Open Government Licence Canada (to confirm) |
| Eurasia | DATED-1, Hughes et al. (2016), PANGAEA 10.1594/PANGAEA.848117 | 25-10 ka every 1 kyr (min/mc/max polygons); 27, 28-29, 30-32 ka (lines + polygons); 38-35 ka (min/max lines only) | 1 kyr | shapefiles, north-pole LAEA | no | CC BY 3.0 |
| NH | Batchelor et al. (2019) | MIS 4 (65 ka), 45, 40, 35, 30 ka, LGM (+ older stages) | snapshots, min/best/max | shapefiles | yes: `batchelor2019_nh_ice`, IceHistory-Batchelor2019 | none stated (article CC BY 4.0) |
| Greenland | PaleoGrIS 1.0, Leger et al. (2024) | 13-7 ka (deglaciation) | isochrones, 1 km age buffers | NetCDF + shapefiles | yes: `leger2024_paleogris` | CC BY 4.0 |
| Greenland, Antarctica, present | BedMachine / Bedmap / RTopo-2 grounded mask | 0 ka | raster | raster | yes | various |
| Antarctica | RAISED Consortium (2014) | 20, 15, 10, 5 ka, scenarios A/B | 5 kyr | maps; GIS files not found online (ask the authors) | no | article CC BY |

Gaps:
- MIS 3-4 (80-25 ka): only Gowan is continuous; Batchelor has 5 snapshots with min/max (sparse, nominal
  ages); DATED-1 covers Eurasia back to 38 ka. No independent MIS 3-4 chronology for Greenland
  (Gowan's buffers) or Antarctica (Gowan, drawn from scratch).
- Greenland LGM-deglaciation: NADI-1 leaves it out; PaleoGrIS starts at 13 ka; Gowan uses Funder (2011).
- Antarctica after the LGM: RAISED only at 4 slices; Gowan's polygons are the only digital set at 2.5 kyr.
- Gowan's polygons are "original" in the sense that they are his input, not ICESHEET output; they are
  themselves a compilation of the above (NA before NADI-1, Eurasia from DATED-1, Antarctica from RAISED).

## 3. Workflow

### 3.1 Product: polygon margin chronology (FesmData)

- Prepared file: the polygons themselves, per slice (lon-lat vertices, ring offsets, slice index, source and
  region attributes), from the ICESHEET archive (and later NADI-1, DATED-1, ...). No gridding in prepare.
- A polygon "remap" method (fesmdata.jl) instead of the conservative remap: for each target grid, project the
  vertices with Proj to the grid's projection (`crs:proj_params`) and compute per slice
  - `sdist(x, y, time)` [m, + in ice]: exact Euclidean distance from the cell centre to the boundary of the
    union of the slice's polygons, in the grid projection;
  - `f_ice(x, y, time)` [1]: the area fraction of the union in the cell (scanlines, exact along x).
  `sdist` cannot be remapped (it is a property of the target grid), so it has to be computed there.
- Union, not per-region sums: Gowan's polygons overlap (NA 20 ka: 2022 of 12001 edges lie inside another
  polygon: Laurentide/Innuitian/Greenland). The prototype drops the edges whose midpoint is inside another
  polygon and fills the union along scanlines. `max(sd_A, sd_B)` is only exact for disjoint sets.
- Scenarios as products, as now: minimal (`margins/`) and maximal (`margins/a1/` over `margins/`).
- Licence GPL-3.0-only (ICESHEET 2.0), like PaleoMIST_taub.
- Cost (prototype, brute force over all edges, 16 threads): 0.1-1.6 s per region and slice at 32 km.
  Finer grids scale with cells x edges (8 km: ~16x); a band + EDT far field or a spatial index if needed.

### 3.2 Exact sdist vs the backend's EDT of a thresholded mask

Prototype numbers (section 4): the slice masks are the same to within ~0.3 % of the cells, so the gain is
not at the slice times' masks. It is in

- the sub-cell position: the exact sdist spreads the margin within the cells (5-95 % 0.03-0.88), the EDT
  puts it at the cell faces (0.5) at every slice time;
- the interpolation between slices: halfway between 20 and 17.5 ka the exact sdist gives 18.212 Mkm2
  (linear in area: 18.193), the EDT of PaleoMIST's mask 18.077 (-0.6 %, 178 + 30 cells differ);
- robustness: the EDT of a thresholded mask jumps where a small ice cap or an isolated cell appears or
  disappears with the threshold (20 ka, far from the margin: rms 38 km, 99 % 207 km between the EDTs of
  the PaleoMIST and the polygon masks); exact sdist sees sub-cell polygons (2 at 20 ka, 16 at 0 ka in
  Eurasia) as near-zero distances.

Caveat, in the library: `plastic_margin_frac_from_sdist` takes frac = -sd/h, the Euclidean distance as if
it were along the axis. With an exact sd that is biased towards more ice for oblique margins (mean frac
0.42 instead of 0.5, about 0.08 cell, ~2.6 km at 32 km). The zero crossing of sd between the ice-free
cell and its ice neighbour, frac = sd_f / (sd_f - sd_n) (min over the ice neighbours), is unbiased (mean
0.49-0.56) and gives exactly 0.5 for the EDT sd at the slice times (no change for the current product).
So: exact sdist only together with this change in PlasticIceSheet (Fortran and Julia).

### 3.3 Backend

- No change needed to use the product's `f_ice`: `margin_var = "f_ice"`, `margin_H_min = 0.5` (EDT) as now.
- To use `sdist`: one key saying what the variable is, e.g. `margin_kind = "threshold"` (default, as now:
  var > `margin_H_min`, then EDT) or `"sdist"` (read the slices as signed distances [m], mask = sd > 0).
  `margin_init` dispatches on it; both fill `mg%sd`, so `margin_at`, the footprint and `margin_position`
  are unchanged. With `"sdist"`, `margin_H_min` is not used (error if set?).

### 3.4 Combining regions and sources

- v1: Gowan's polygons only, all regions: the margins PaleoMIST was made with, consistent with Gowan's tau_b
  (tuned to them), and the prototype reproduces IceHistory-PaleoMIST `f_ice` (rms 0.02-0.04).
- Later, per region and time window, a table in prepare (e.g. NA 25-1 ka: NADI-1 optimal; Eurasia
  25-10 ka: DATED-1 most-credible; else Gowan), with:
  - one time axis: the 2.5 kyr axis (NADI-1 has every 2.5 kyr multiple; DATED-1 every 1 kyr from 25 ka:
    20, 15, 10 ka directly, 22.5, 17.5, 12.5 ka in between as Gowan did) or the union of the slices; all
    sources in calendar years (avoid Dyke 2004 and Dalton 2020, in 14C ka);
  - spatial seams: NADI-1 has no Greenland, so Greenland stays Gowan's; clipping polygons by region needs
    polygon boolean operations (LibGEOS or GeometryOps: a new dependency, to ask);
  - a jump at a window boundary becomes a fast migration over one interval (the backend interpolates);
  - tau_b re-calibration (ES-MDA) with any margin other than Gowan's.
- min/max margins (NADI-1, DATED-1, Batchelor) could become margin scenarios for the ensemble.

### 3.5 Marine margins

Gowan's polygons are grounded-ice extents, marine-based ice included (shelf edge at the LGM: Labrador,
Barents-Kara, Antarctic shelf from RAISED grounding lines; present = RTopo-2 / Bedmap2 grounding lines), no
ice shelves. The product stays geometric (no bed, no sea level). Whether imposed ice below flotation stays
grounded (ICESHEET's assumption) or is removed online (VILMA2 bed, sea level) is a backend choice.

### 3.6 tau_b product

- Same archive, same reader. The domain polygons are projected (QGIS, x, y from each region's projection
  centre) while the margins are lon-lat: one polygon module (shared/) reading both, inverse-projecting the
  domains to lon-lat, then projecting to the target grid.
- `Gowan2021_PaleoMIST_taub/prepare.jl` still reads PaleoMIST 1.0 output (its 5 km grids and `ice_any` for
  domain 0 and the 150 km extension, `EXTEND_KM`). With the margin product: rasterise the domains directly
  on the kryos grid (as create_ss_grid: last polygon wins; at the cell centre or by dominant area) and
  extend to the margin footprint (union of all slices) instead of PaleoMIST's ice: the tau_b product then
  depends on ICESHEET only, and `EXTEND_KM` (needed for PaleoMIST's smoothed ice) goes.
- Natural home: one dataset for both (e.g. `Gowan2021_ICESHEET`, products `margins`, `margins_max`, `taub`).

## 4. Prototype (Levante, /work/ba1442/robinson/margin-proto)

`src/margin_proto.jl` (Julia 1.13, a copy of the env of Gowan2021_PaleoMIST_taub: NCDatasets, Proj; no new
package), `src/run_proto.sh` (sbatch, 16 cores, 36 s), log `out/proto_28838090.log`, fields
`out/margin_proto_NH-32KM.nc` (sdist, f_ice of the NA+Eurasia union; NA alone; EDT of PaleoMIST
`f_ice > 0.5`; PaleoMIST f_ice). Inputs copied (cp -L) to `data/`: `margins/{20000,17500,0}.gmt` of
North_America and Eurasia (ICESHEET 641d172), NH-32KM_IceHistory-PaleoMIST.nc (v2.3.0 on disk).
Comparison on the cells nearer to the NA polygons than to the Eurasian ones; areas with the grid file's
cell areas (600-1089 km2; dx dy = 1024).

| | 20 ka | 17.5 ka | 0 ka |
|---|---|---|---|
| area, polygon f_ice / PaleoMIST f_ice [Mkm2] | 18.436 / 18.436 | 17.968 / 17.968 | 1.931 / 1.930 |
| area, polygon centre mask / PaleoMIST f_ice > 0.5 | 18.429 / 18.418 | 17.958 / 17.958 | 1.946 / 1.914 |
| cells: centre mask / PaleoMIST > 0.5 | 19837 / 19827 | 19257 / 19256 | 1885 / 1852 |
| mismatched cells, centre mask vs PaleoMIST (only poly / only PM) | 58 / 48 | 53 / 52 | 67 / 34 |
| mismatched cells, polygon f > 0.5 vs PaleoMIST > 0.5 | 25 / 25 | 31 / 30 | 19 / 14 |
| f_ice polygon - PaleoMIST, rms (cells with ice) | 0.024 | 0.024 | 0.041 |
| sdist exact - EDT(PaleoMIST), within 2 cells of the margin: mean / rms [km] | +2.2 / 23.5 | +0.2 / 10.0 | +4.5 / 17.4 |
| same, farther (and > 10 cells from the grid edge): mean / rms | +6.1 / 38.0 | -1.3 / 9.8 | +11.4 / 22.0 |
| sdist exact - EDT(polygon centre mask), within 2 cells: mean / rms | +0.5 / 10.7 | +0.3 / 8.8 | +0.6 / 12.2 |
| margin frac: -sd/h (exact) / zero crossing (exact) / EDT | 0.42 / 0.50 / 0.50 | 0.42 / 0.49 / 0.50 | 0.40 / 0.56 / 0.50 |

- PaleoMIST's `f_ice` is Gowan's polygons to within resampling (two resamplings: 5 km grdmask, 0.25 deg,
  conservative remap): switching to the polygons changes the slice masks by ~0.3-0.5 % of the ice cells.
- The sdist offsets near the margin are the half-cell staircase of the EDT (rms ~10 km on the same mask)
  plus the mask differences (sign differs on ~100 cells per slice); far from it the 20 and 0 ka numbers are
  dominated by small ice caps present in one mask only.
- The grid projection inflates areas by ~10 % at 20 ka (polygon area 20.41 Mkm2 in the stereographic
  plane, 18.53 equal-area); f_ice times the grid file's cell areas gives the true area.
- Edges: median 0.1-2.6 km, but up to 125-147 km in NA (straight lines across the shelves). GMT projected
  the vertices only (straight in the regional LAEA); here straight in the polar stereographic plane; the
  difference is far below a cell. Densify in lon-lat if finer grids need it.

## 5. Open questions

1. Source policy: v1 Gowan-only (exactly PaleoMIST's margins, consistent with his tau_b), or replace
   NA 25-1 ka by NADI-1 and Eurasia 25-10 ka by DATED-1 now (needs Greenland seams, polygon clipping
   dependency, tau_b re-calibration)?
2. FesmData structure: polygons as the prepared file plus a polygon remap method (sdist, f_ice) in
   fesmdata.jl, and one dataset for margins and tau_b (dropping taub's PaleoMIST dependency)?
3. Backend: add `margin_kind = "sdist"`, or stay with `f_ice > 0.5` + EDT on the new product?
4. Library: change `plastic_margin_frac_from_sdist` to the zero crossing (needed for exact sdist; identical
   for the EDT at slice times)?
5. Marine ice below flotation online: keep grounded as imposed, or drop?
6. Antarctica: RAISED GIS files are not online as far as found; ask the authors or keep Gowan's?
