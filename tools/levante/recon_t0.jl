# t0 check of the PaleoMIST start (par/kryos_NH-ANT_vilma2.nml): per domain, the ice of the run at its first output record against the
# domain's PaleoMIST history at record 1 (-80 ka), on the ice grid: max |dH|, and the volumes as sea-level equivalent
# (rho_i/rho_sw V / A_sea of the constants); and the gathered load of the log (domain ESL against the reference) for reference.
#   julia --project=tools tools/levante/recon_t0.jl <run> <hub_areas.x> [--rho-i 910 --rho-sw 1028 --asea 3.618e14]
using NCDatasets, Printf
function main(args)
    rdir, areasx = args[1:2]
    ri = 910.0; rsw = 1028.0; asea = 3.618e14
    for (f, v) in (("--rho-i", :ri), ("--rho-sw", :rsw), ("--asea", :asea))
        i = findfirst(==(f), args); i === nothing && continue
        x = parse(Float64, args[i+1]); v == :ri ? (ri = x) : v == :rsw ? (rsw = x) : (asea = x)
    end
    maps = realpath(joinpath(rdir, "maps"))
    fd = "/work/ba1442/fesmdata"
    doms = [("North", "NH-32KM", "ice_nor"), ("Antarctica", "ANT-32KM", "ice_ant")]
    areas = tempname() * ".nc"
    run(`$areasx $maps $areas $([d[2] for d in doms])`)
    for (d, g, ice) in doms
        A = NCDataset(areas) do a; Float64.(a["area_$g"][:, :]); end
        pm = NCDataset(joinpath(fd, d, g, "$(g)_IceHistory-PaleoMIST.nc")) do p
            t = Float64.(Array(p["time"])); k = argmax(t)            # years ago: the largest is -80 ka
            (k, t[k], Float64.(p["H_ice"][:, :, k]))
        end
        y = NCDataset(joinpath(rdir, ice * ".nc")) do s
            t = Float64.(Array(s["time"])); (t[1], Float64.(s["H_ice"][:, :, 1]))
        end
        Hp = replace(pm[3], NaN => 0.0); Hy = y[2]
        dH = Hy .- Hp
        Vp = sum(Hp .* A); Vy = sum(Hy .* A)
        esl(V) = -(ri / rsw) * V / asea
        @printf("%s (%s): PaleoMIST record %d (%g years ago) vs the run's first record (t = %g yr): max |dH| %.3f m (rms %.3f m); volume %.4e vs %.4e m3; ESL %.4f vs %.4f m, difference %.4f m\n",
                d, g, pm[1], pm[2], y[1], maximum(abs.(dH)), sqrt(sum(dH .^ 2) / length(dH)), Vp, Vy, esl(Vp), esl(Vy), esl(Vy) - esl(Vp))
    end
    for l in eachline(joinpath(rdir, "out.out"))
        occursin("sea-level equivalent of the ice change", l) && (println(strip(l)); break)
    end
end
main(ARGS)
