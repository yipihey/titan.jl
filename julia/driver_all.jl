# Run a benchmark through the Julia step driver (all-Fortran routines) and
# byte-compare TITAN's output files with the standalone reference runs.
# usage: julia driver_all.jl NAME DECK [predump_from_reference_dir]
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver; const TD = TitanDriver; const TF = TD.TF
TitanDriver.bitwise!(true)     # parity/byte-compare script: reproduce TITAN 1.2 bitwise
name, deck = ARGS[1], ARGS[2]
pre = length(ARGS) >= 3 ? ARGS[3] : nothing
ref = joinpath(TF.ROOT, "reference_runs", name)
dir = joinpath(TF.ROOT, "lib", "p3_runs", "drv_$name")
stem = name == "RLST_leg2" ? "RLST" : name
if name == "COOL"            # COOL reads its restart model from COOL.dmp = HEAT's dump file
    mkpath(dir); cp(joinpath(TF.ROOT, "reference_runs", "HEAT", "HEAT.dmp"), joinpath(dir, "COOL.dmp"); force = true)
    pre = nothing
end
TF.init(deck; rundir = dir, predump = pre === nothing ? nothing : joinpath(TF.ROOT, "reference_runs", pre, "$stem.dmp"))
t0 = time(); n = 0
while true; global n += 1; TD.titan_step!() != 0 && break; end
wall = time() - t0; TF.finalize()
res = String[]
for e in ("dmp", "out", "nrg", "sty", "tty")
    f = joinpath(dir, "$stem.$e"); g = joinpath(ref, "$stem.$e")
    isfile(f) && isfile(g) && push!(res, "$e:" * (read(f) == read(g) ? "identical" : "DIFF"))
end
println("$name: steps driven = $n, wall = $(round(wall, digits = 2)) s; ", join(res, " "))
