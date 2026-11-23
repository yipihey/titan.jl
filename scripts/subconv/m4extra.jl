# m4extra.jl <run dir> <out dir>: runrec.jl export of a run directory into a separate output directory (lane dirs stay untouched).
# Usage (from julia/): julia -t 2 --project=<env> ../scripts/subconv/m4extra.jl <rundir> <outdir>
include(joinpath(@__DIR__, "..", "..", "julia", "runrec.jl"))
length(ARGS) == 2 || error("usage: m4extra.jl <rundir> <outdir>")
if isempty(filter(f -> occursin(r"^snaps_\d+_\d+\.jls$", f), readdir(ARGS[1])))
    println("no snapshots yet in $(ARGS[1]); nothing exported"); exit(0)
end
export_frames(ARGS[1]; out = ARGS[2])
