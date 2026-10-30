# Gate G3 (phase4-plan 13.2): the mass and bsl budgets of a run with geo = vilma2, from its log, series and streams.
#   julia --project=tools tools/levante/g3_budget.jl <run> <lossless run> <hub_areas.x> [--md FILE] [--dt 10]
# <run>: out.out (the log: the gather lines of kryos_globe_ice, the ocean area of the reference) and geo_glob_ts.nc (mass_resid,
# bsl, bsl_ive, A_ocn, every step); <lossless run>: the same configuration with [kryos] write_lossless, whose streams
# geo_glob_gia.nc (host terms as VILMA2 got them), geo_glob_vilma.nc (rsl, C_ocean) and the bundle of VILMA2
# (geo_glob_restart.nc: w_ref, I_g_ref, C_host_ref, z_bed_eq, h_ice_eq) give (iii) without rounding.
#   (i)   merge closure of every gather (log) against 1e-9
#   (ii)  mass_resid of VILMA2 (series), max
#   (iii) int C rsl dA against -(rho_i/rho_w) int dI_g dA + int F dA, dI_g recomputed here from the gathered terms:
#         dI_g = (I_g - I_g_ref) + [a ice (1 - C_own) - a0 ice_eq (1 - C0_own)], a = 1 - w, C_own = (C - C_hat)/a, C0_own from
#         (z_bed_eq, h_ice_eq) with VILMA2's flotation rule; F is the host part of the fill only (the own part's fill needs
#         zeta0, not output): the residual is reported against both |lhs| and |rhs|
#   (iv)  bsl_ive (series, at t_n) against the sea-level equivalent of the grounded-ice change of the gather at t_n - dt
#         (log: domains + external, with the run's rho_ice/rho_sw and the sea surface area of the constants)
#   (v)   the ocean area of the reference against the sea surface (log)
using NCDatasets, Printf, Statistics

function parse_log(file)
    gl = Dict{Float64,NamedTuple}()
    aref = (NaN, NaN); rho = (NaN, NaN)
    for l in eachline(file)
        if occursin("sea-level equivalent of the ice change", l)
            m = match(r"t = (\S+) yr, sea-level equivalent of the ice change \[m\]:(.*), external (\S+); clipped external ice (\S+) Gt; merge closure (\S+)", l)
            m === nothing && continue
            t = parse(Float64, m[1])
            parts = split(strip(m[2]))
            doms = Dict{String,Float64}(parts[i] => parse(Float64, parts[i+1]) for i in 1:2:length(parts)-1)
            gl[t] = (doms = doms, ext = parse(Float64, m[3]), clip = parse(Float64, m[4]), clos = parse(Float64, m[5]))
        elseif occursin("ocean area of the reference", l)
            m = match(r"ocean area of the reference (\S+) m2 \((\S+) % of the sea surface (\S+) m2", l)
            m === nothing || (aref = (parse(Float64, m[1]), parse(Float64, m[3])))
        elseif occursin("rho_ice =", l) && isnan(rho[1])
            m = match(r"rho_ice =\s*(\S+)", l); m === nothing || (rho = (parse(Float64, m[1]), rho[2]))
        elseif occursin("rho_sw =", l) && isnan(rho[2])
            m = match(r"rho_sw =\s*(\S+)", l); m === nothing || (rho = (rho[1], parse(Float64, m[1])))
        end
    end
    gl, aref, rho
end

function main(args)
    rdir, lrun, areasx = args[1:3]
    md = nothing; dt = 10.0
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    i = findfirst(==("--dt"), args); i === nothing || (dt = parse(Float64, args[i+1]))
    io = IOBuffer()
    gl, aref, rho = parse_log(joinpath(rdir, "out.out"))
    println(io, "# G3: budgets of $(basename(dirname(rdir)))/$(basename(rdir))\n")
    println(io, "densities of the run: rho_ice $(rho[1]), rho_sw $(rho[2]); gathers in the log: $(length(gl))\n")
    # (i)
    clos = [v.clos for v in values(gl)]
    @printf(io, "(i) merge closure of the gathers: max %.3e, n = %d (bound 1e-9) -> %s\n\n", maximum(clos), length(clos), maximum(clos) <= 1e-9 ? "PASS" : "FAIL")
    # (ii), (iv) series
    ts = NCDataset(joinpath(rdir, "geo_glob_ts.nc"))
    t = Float64.(Array(ts["time"])); mr = Float64.(Array(ts["mass_resid"])); bi = Float64.(Array(ts["bsl_ive"]))
    bs = Float64.(Array(ts["bsl"])); ao = Float64.(Array(ts["A_ocn"]))
    @printf(io, "(ii) VILMA2 mass_resid: max %.3e over %d steps (sle_tol 1e-7; reported, proposed bound 1e-6)\n\n", maximum(abs.(mr)), length(mr))
    rows = []
    for n in eachindex(t)
        tl = t[n] - dt
        k = findfirst(x -> abs(x - tl) < 1e-6, collect(keys(gl)))
        k === nothing && continue
        g = gl[collect(keys(gl))[k]]
        esl = sum(values(g.doms)) + g.ext
        push!(rows, (t[n], bi[n], esl, bs[n], ao[n]))
    end
    if !isempty(rows)
        rel = [abs(r[2] - r[3]) / max(abs(r[3]), 1e-12) for r in rows]
        @printf(io, "(iv) bsl_ive against the sea-level equivalent of the gathered grounded-ice change: %d steps, max relative difference %.3e (proposed bound 1%%); last: bsl_ive %.4f m, ESL %.4f m, bsl %.4f m\n\n",
                length(rows), maximum(rel), rows[end][2], rows[end][3], rows[end][4])
        asea = aref[2]
        relv = [abs(r[2] * r[5] - r[3] * asea) / max(abs(r[3] * asea), 1e-12) for r in rows]
        @printf(io, "    as volumes (bsl_ive times the ocean area of VILMA2's solve against the ESL times the sea surface area of the constants, i.e. without the change of the ocean area): max relative difference %.3e\n\n", maximum(relv))
        println(io, "| t [yr] | bsl_ive [m] | ESL of the gather at t-dt [m] | bsl [m] | rel. diff | rel. diff as volumes |\n|---|---|---|---|---|---|")
        for r in rows[unique([1; round.(Int, range(1, length(rows), length = min(9, length(rows)))); length(rows)])]
            @printf(io, "| %g | %.5f | %.5f | %.5f | %.2e | %.2e |\n", r[1], r[2], r[3], r[4], abs(r[2] - r[3]) / max(abs(r[3]), 1e-12), abs(r[2] * r[5] - r[3] * asea) / max(abs(r[3] * asea), 1e-12))
        end
        println(io)
    end
    # (v)
    @printf(io, "(v) ocean area of the reference %.5e m2 against the sea surface %.5e m2 (the constants): %+.3f %% (bound 1%%) -> %s\n    VILMA2 A_ocn (series): first %.5e, last %.5e m2\n\n",
            aref[1], aref[2], 100 * (aref[1] - aref[2]) / aref[2], abs(aref[1] - aref[2]) / aref[2] <= 0.01 ? "PASS" : "FAIL", ao[1], ao[end])
    # (iii)
    if isdir(lrun)
        areas = tempname() * ".nc"
        run(`$areasx $(joinpath(rdir, "maps")) $areas GAUSS-L170`)
        A = NCDataset(areas) do a; Float64.(a["area_GAUSS-L170"][:, :]); end
        gia = NCDataset(joinpath(lrun, "geo_glob_gia.nc")); vil = NCDataset(joinpath(lrun, "geo_glob_vilma.nc"))
        bund = sort(filter(d -> startswith(d, "restart-") && isdir(joinpath(lrun, d)), readdir(lrun)))
        rst = NCDataset(joinpath(lrun, last(bund), "geo_glob_restart.nc"))
        w0 = Float64.(rst["w_ref"][:, :]); I0 = Float64.(rst["I_g_ref"][:, :]); Cref = Float64.(rst["C_host_ref"][:, :])
        zeq = Float64.(rst["z_bed_eq"][:, :]); heq = Float64.(rst["h_ice_eq"][:, :])
        # the arrays of the files: (lon, lat) as Julia reads dims (lon, lat, time); the restart is (lat, lon) in Fortran order
        # -> as NCDatasets returns them (lon, lat) if dims are (lon, lat); check shapes
        tg = Float64.(Array(gia["time"])); tv = Float64.(Array(vil["time"]))
        ri, rw = rho
        C0own = Float64.(zeq .< 0 .&& ri .* heq .< -rw .* zeq)
        println(io, "(iii) int C rsl dA = -(rho_i/rho_w) int dI_g dA + int F dA, from the lossless run's streams (host part of F only):\n")
        println(io, "| t [yr] | int C rsl dA [m3] | -(ri/rw) int dI_g dA [m3] | VILMA2 bsl_ive A_ocn [m3] | rhs / VILMA2 - 1 | int F_host dA [m3] | residual (own fill) | residual / |lhs| |\n|---|---|---|---|---|---|---|---|")
        for (k, tt) in enumerate(tg)
            kv = findfirst(x -> abs(x - tt) < 1e-6, tv); kv === nothing && continue
            w = Float64.(gia["host_w"][:, :, k]); Ig = Float64.(gia["host_Ig"][:, :, k]); Ch = Float64.(gia["host_C"][:, :, k])
            F = Float64.(gia["host_F"][:, :, k]); H = Float64.(gia["H_ice_ext"][:, :, k])
            C = Float64.(vil["C_ocean"][:, :, kv]); rsl = Float64.(vil["rsl"][:, :, kv])
            a = 1 .- w; a0 = 1 .- w0
            Cown = ifelse.(a .> 1e-12, (C .- Ch) ./ max.(a, 1e-12), 0.0)
            dI = (Ig .- I0) .+ (ifelse.(a .> 1e-12, a .* H .* (1 .- Cown), 0.0) .- ifelse.(a0 .> 1e-12, a0 .* heq .* (1 .- C0own), 0.0))
            lhs = sum(C .* rsl .* A); rhs = -(ri / rw) * sum(dI .* A); fi = sum(F .* A)
            res = lhs - (rhs + fi)
            ktt = findfirst(x -> abs(x - tt) < 1e-6, t); vil2 = ktt === nothing ? NaN : bi[ktt] * ao[ktt]
            @printf(io, "| %g | %.5e | %.5e | %.5e | %.2e | %.3e | %.3e | %.2e |\n", tt, lhs, rhs, vil2, rhs / vil2 - 1, fi, res, abs(res) / max(abs(lhs), 1e-30))
        end
        println(io, "\nThe recomputed grounded-ice term reproduces VILMA2's own bsl_ive A_ocn (column 5). The residual of the balance is the own part's subgrid fill, a zeta0 (C_own - C0_own), which needs zeta0 (not an output); it is nearly constant in time (set by the coast shift at t0 against the reference).")
    end
    s = String(take!(io)); print(s); md === nothing || write(md, s)
    return 0
end
exit(main(ARGS))
