# The spread of the global sea-level records that measure the same quantity as the scored stack, the structural part of the ESL
# group's model error (plastic_ensemble.jl [esl] spread_file): per age the sample sd s(t) of the central values of two records of
# the ice mass, each as its authors give it (pm_compare.jl sealevel_at, no re-zeroing):
#   Spratt2016               sl_short, the stack of δ18O-based reconstructions (0 at 5 ka): the scored record
#   Lambeck2014              sl, far-field RSL inverted with GIA (ice-volume equivalent)
# Not Grant et al. 2014 (the Red Sea: site RSL, its GIA not removed; scored at its site by the far field groups) nor Waelbroeck et
# al. 2002 (benthic δ18O regressed on corals, a predecessor of the stack's own inputs): with them s(t) was 10-16 m at 25-15 ka and
# the group no longer told an 83-m LGM from a 110-m one (scores v4). s(t) = |Spratt - Lambeck| / √2, every 0.5 kyr from 0 to 33 ka
# (Lambeck et al. 2014 end there), the variance averaged over ±1 kyr (the stack's own scatter on 1-2 kyr is not the disagreement
# of the records' level). The group's model error is σ_m(t) = sqrt(σ_m^2 + s(t)^2). Computed once; the output is an input of the
# scores, kept with this tool:
#   julia --project=tools/levante tools/levante/esl_spread.jl <fesmdata> [<out.csv>]     default out: tools/levante/esl_spread.csv
using Printf, Statistics
include(joinpath(@__DIR__, "pm_compare.jl"))

const RECORDS = [("Spratt2016", "sl_short"), ("Lambeck2014", "sl")]

function esl_spread(fes, out; amax = 33.0, step = 0.5, half = 1.0)
    ages = collect(0.0:step:amax)
    X = hcat([sealevel_at(fes, p, 1000 .* ages; var = v)[1] for (p, v) in RECORDS]...)       # age x record
    var1 = [(x = filter(isfinite, X[i, :]); length(x) >= 2 ? var(x) : NaN) for i in eachindex(ages)]
    k = round(Int, half / step)
    s = [sqrt(mean(filter(isfinite, var1[max(1, i - k):min(end, i + k)]))) for i in eachindex(ages)]
    open(out, "w") do io
        println(io, "# The spread of the data-based ESL records (tools/levante/esl_spread.jl, ", Dates_now(), "; FesmData ", fes, "): s_m, the sd",
                " of the records' central values as given, the variance averaged over ±", half, " kyr; the records' values [m] at the age")
        println(io, "age_ka,n,s_m,", join(first.(RECORDS), ","))
        for i in eachindex(ages)
            println(io, @sprintf("%.1f,%d,%.2f,", ages[i], count(isfinite, X[i, :]), s[i]), join([@sprintf("%.2f", x) for x in X[i, :]], ","))
        end
    end
    println("wrote ", out)
    for a in (25.0, 22.5, 20.0, 17.5, 15.0, 12.5, 10.0)
        i = findfirst(==(a), ages); println(@sprintf("%5.1f ka: s %.1f m, sigma_m(t) with 5 m %.1f m", a, s[i], hypot(5.0, s[i])))
    end
end
Dates_now() = (t = read(`date -I`, String); strip(t))

if abspath(PROGRAM_FILE) == @__FILE__
    esl_spread(ARGS[1], length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "esl_spread.csv"))
end
