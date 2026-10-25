"""
    KryosTools

Julia tooling of kryos phase 1 (regression-design.md §2-§4):

- `KryosTools.Compare`: two run directories compared bit for bit through a
  mapping table (`tools/gate/map_phase1.toml`) or as two runs of one model.
- `KryosTools.Convert`: a yelmox par file converted into a kryos par file
  (core-design.md §5), fully expanded, with its accounting.
- `KryosTools.Gate`: the stack file, the gate's directories and `-p` lists,
  the `docs/gate.md` row.

Scripts: `tools/bin/kryos-compare.jl`, `tools/bin/yelmox2kryos.jl`. Plots
(`Compare.plot_variable`) need `using CairoMakie` (a package extension), so
the comparison itself runs without Makie.
"""
module KryosTools

const TOOLS_DIR = normpath(joinpath(@__DIR__, ".."))

include("Compare.jl")
include("Convert.jl")
include("Gate.jl")

using .Compare, .Convert, .Gate

end # module
