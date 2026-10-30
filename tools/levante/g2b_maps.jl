# G2b diagnosis maps from g2b_diag.jl's g2b_diag_<tag>.nc (any machine with CairoMakie and NCDatasets).
#   julia g2b_maps.jl <g2b_diag_<tag>.nc> <out folder> [label]
# Writes g2b_<label>_NH.png and g2b_<label>_ANT.png: per time (columns) the difference kryos - standalone of dz_bed [m], of the ocean
# function C, of the load L [m of sea water], its ice part (rho_i/rho_w) dI_g and the water part; polar azimuthal equidistant
# projection, the host region (w > 0.5) outlined.
using NCDatasets, CairoMakie, Printf

function polar(lon, lat, hemi)
    c = hemi == :N ? 90.0 .- lat : 90.0 .+ lat              # colatitude from the pole [deg]
    s = hemi == :N ? 1.0 : -1.0
    X = [c[j] * cosd(lon[i]) for i in eachindex(lon), j in eachindex(lat)]
    Y = [s * c[j] * sind(lon[i]) for i in eachindex(lon), j in eachindex(lat)]
    X, Y
end

function main(args)
    f, out = args[1], args[2]; label = length(args) >= 3 ? args[3] : splitext(basename(f))[1]
    mkpath(out)
    ds = NCDataset(f)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); t = Float64.(ds["time"][:]); w = Float64.(ds["host_w"][:, :])
    g(n) = Float64.(ds[n][:, :, :])
    # close the longitude circle
    cyc(a) = cat(a, a[1:1, :, :]; dims = 1)
    lonc = vcat(lon, lon[1] + 360)
    fields = [("dz_bed", cyc(g("dz_bed_k") .- g("dz_bed_s")), 150.0, "m"),
              ("C", cyc(g("C_k") .- g("C_s")), 1.0, "1"),
              ("load L", cyc(g("L_k") .- g("L_s")), 150.0, "m w.e."),
              ("ice part", cyc((910 / 1028) .* (g("dIg_k") .- g("dIg_s"))), 150.0, "m w.e."),
              ("water part", cyc(g("W_k") .- g("W_s")), 150.0, "m w.e."),
              ("int dL dt", cyc(Float64.(ds["int_dL"][:, :, :])) ./ 1000, 150.0, "m w.e. kyr")]
    wc = cat(w, w[1:1, :]; dims = 1)
    for (hemi, cmax, nm) in ((:N, 60.0, "NH"), (:S, 40.0, "ANT"))
        js = hemi == :N ? findall(lat .> 90 - cmax) : findall(lat .< -90 + cmax)
        X, Y = polar(lonc, lat[js], hemi)
        fig = Figure(size = (300 * length(t) + 120, 290 * length(fields)), fontsize = 13)
        for (row, (fname, F, lim, unit)) in enumerate(fields)
            local hm
            for (col, tt) in enumerate(t)
                ax = Axis(fig[row, col], aspect = DataAspect(), title = row == 1 ? @sprintf("t = %.0f yr", tt) : "")
                hidedecorations!(ax); hidespines!(ax)
                hm = surface!(ax, X, Y, zeros(size(X)); color = F[:, js, col], colormap = :balance, colorrange = (-lim, lim),
                              shading = NoShading, rasterize = 2)
                contour!(ax, X, Y, wc[:, js]; levels = [0.5], color = :black, linewidth = 0.6)
                for c in 10:10:Int(cmax)                     # latitude circles every 10 degrees
                    lines!(ax, c .* cosd.(0:2:360), c .* sind.(0:2:360), color = (:grey50, 0.5), linewidth = 0.4)
                end
                limits!(ax, -cmax, cmax, -cmax, cmax)
            end
            Colorbar(fig[row, length(t)+1], hm, label = fname * " [" * unit * "]")
        end
        Label(fig[0, :], "G2b " * label * ", " * nm * ": kryos - standalone", fontsize = 16)
        save(joinpath(out, "g2b_$(label)_$(nm).png"), fig, px_per_unit = 1)
    end
    close(ds)
end
main(ARGS)
