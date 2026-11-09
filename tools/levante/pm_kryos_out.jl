# The VILMA2 stream of a kryos run (geo = vilma2, out_vilma: <instance>_vilma.nc) as the vilma_out.nc of a standalone pm_run.sh run,
# so that pm_compare.jl scores both the same way:
#   julia --project=tools tools/levante/pm_kryos_out.jl <run folder> <out folder> [<instance> = geo_glob]
# Writes <out folder>/vilma_out.nc with lon, lat, time, rsl, dz_bed, bsl, bsl_ive. The time of a record is moved by +1 yr except the
# first (pm_compare.jl's trec: record t_k of a standalone run is the state at t_k + 1 yr, its 80 ka record the spun-up start); a kryos
# record at t_k is the state after the step ending at t_k, the same state. plastic_esmda.jl includes it (kryos_vilma_out).
using NCDatasets

function kryos_vilma_out(run, out; inst = "geo_glob")
    src = joinpath(run, inst * "_vilma.nc"); mkpath(out); dst = joinpath(out, "vilma_out.nc"); isfile(dst) && rm(dst)
    NCDataset(src) do s
        t = Float64.(s["time"][:]); t2 = copy(t); t2[2:end] .+= 1.0
        NCDataset(dst, "c") do d
            defDim(d, "lon", s.dim["lon"]); defDim(d, "lat", s.dim["lat"]); defDim(d, "time", length(t))
            defVar(d, "lon", Float64.(s["lon"][:]), ("lon",)); defVar(d, "lat", Float64.(s["lat"][:]), ("lat",))
            defVar(d, "time", t2, ("time",))
            for v in ("rsl", "dz_bed"); defVar(d, v, Float64.(nomissing(s[v][:, :, :], NaN)), ("lon", "lat", "time")); end
            for v in ("bsl", "bsl_ive"); defVar(d, v, Float64.(s[v][:]), ("time",)); end
            d.attrib["source"] = src
        end
    end
    dst
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("wrote ", kryos_vilma_out(ARGS[1], ARGS[2]; inst = length(ARGS) >= 3 ? ARGS[3] : "geo_glob"))
end
