# The bed of the Coulomb law's bed factor lambda (taub_method = "coulomb", taub_lambda_bed = "file", taub_lambda_file;
# docs/backend-plastic.md), one file per grid:
#   julia --project=tools/levante tools/levante/plastic_lambda_bed.jl [<out dir>]
#       writes <out dir>/<grid>_lambda_bed.nc for NH-32KM (GEBCO 2026), ANT-16KM (BedMachine v4) and PAT-PAL-8KM
#       (GEBCO 2026) from the FesmData topography of each grid (default $KRYOS_DATA/plastic; never into $FESMDATA)
#
# z_lambda = z_bed + (rho_i/rho_m) max(H_ice, 0) [m], rho_i 910, rho_m 3300 kg m-3: the present-day bed relaxed by local
# isostasy (the bed without its ice load), relative to present sea level. lambda is a proxy for the geology of the bed
# (marine sediments), so it is one fixed field per grid, the same for the coupled runs on any Earth and for the
# calibration programs (tests/plastic_coulomb_{calib,scan,pd}.f90, lambda_file): until lambda-file the runs took it on
# the reference topography of their Earth (refcorr: + r), the calibration on GEBCO / BedMachine. The field is what the
# calibration programs computed with lambda_bed = "relaxed" (missing values and beds below -9000 m: 0; H_ice < 0: 0).
#
# File: xc, yc (km, from the topography file), z_lambda(xc, yc) (Float64, no fill value), attributes source, rho_i, rho_m.
using NCDatasets, Printf

const RHO_I = 910.0
const RHO_M = 3300.0
const GRIDS = [("North", "NH-32KM", "GEBCO2026"), ("Antarctica", "ANT-16KM", "BedMachine-v4"),
               ("PatagoniaPaleo", "PAT-PAL-8KM", "GEBCO2026")]

# the relaxed present-day bed of a FesmData topography file, as the calibration programs read it
function lambda_bed(topo_file)
    NCDataset(topo_file) do ds
        zb = Float64.(coalesce.(ds["z_bed"][:, :], 0.0))
        H = Float64.(coalesce.(ds["H_ice"][:, :], 0.0))
        nz = count(<(-9000.0), zb)
        zb[zb .< -9000.0] .= 0.0
        z = zb .+ (RHO_I / RHO_M) .* max.(H, 0.0)
        z, Float64.(ds["xc"][:]), Float64.(ds["yc"][:]), nz + count(ismissing, ds["z_bed"][:, :])
    end
end

function write_lambda_bed(file, z, xc, yc, src)
    isfile(file) && rm(file)
    NCDataset(file, "c") do ds
        defVar(ds, "xc", xc, ("xc",); attrib = ["units" => "km", "axis" => "X"])
        defVar(ds, "yc", yc, ("yc",); attrib = ["units" => "km", "axis" => "Y"])
        defVar(ds, "z_lambda", z, ("xc", "yc");
               attrib = ["units" => "m", "long_name" => "bed of the Coulomb law's lambda, relative to present sea level: " *
                         "present-day bed relaxed by local isostasy, z_bed + (rho_i/rho_m) H_ice"])
        ds.attrib["source"] = "kryos tools/levante/plastic_lambda_bed.jl from $src"
        ds.attrib["rho_i"] = RHO_I
        ds.attrib["rho_m"] = RHO_M
    end
    file
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || (println("usage: plastic_lambda_bed.jl [<out dir>]"); exit(1))
    fesmdata = get(ENV, "FESMDATA", "fesmdata")
    out = length(ARGS) == 1 ? ARGS[1] : joinpath(get(ENV, "KRYOS_DATA", "kryos_data"), "plastic")
    isdir(fesmdata) && startswith(realpath(out) * "/", realpath(fesmdata) * "/") && error("$out is in FesmData: not written")
    for (dom, g, topo) in GRIDS
        tf = joinpath(fesmdata, dom, g, "$(g)_TOPO-$(topo).nc")
        z, xc, yc, nz = lambda_bed(tf)
        f = write_lambda_bed(joinpath(out, "$(g)_lambda_bed.nc"), z, xc, yc, basename(tf))
        @printf("%s: %d x %d, z_lambda %.0f to %.0f m, %d cells below -300 m, %d at or above 200 m; %d missing beds set to 0\n",
                f, size(z)..., minimum(z), maximum(z), count(<=(-300.0), z), count(>=(200.0), z), nz)
    end
end
