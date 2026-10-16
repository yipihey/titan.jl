# Phase 4, item 1: run a benchmark leg with the Julia driver in MODE = orig (bitwise
# TITAN 1.2) or fixed (clock + energy-booking fixes), store per-step snapshots.
# usage: julia legs.jl NAME DECK MODE [PREDUMP_NAME]
#   PREDUMP_NAME: leg whose dump (from the same MODE) is the restart model
#   (COOL <- HEAT, RLST_leg2 <- RLST).
using Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF
name, deck, mode = ARGS[1], ARGS[2], ARGS[3]
pre = length(ARGS) >= 4 ? ARGS[4] : nothing
TD.bitwise!(mode == "orig")
root = joinpath(TF.ROOT, "lib", "p4_legs"); dir = joinpath(root, "$(name)_$mode"); mkpath(dir)
if pre !== nothing
    src = pre == "HEAT" ? joinpath(root, "HEAT_$mode", "HEAT.dmp") : joinpath(root, "$(pre)_$mode", "RLST.dmp")
    cp(src, joinpath(dir, name == "COOL" ? "COOL.dmp" : "RLST.dmp"); force = true)
end
const STEPDT = Ref(0.0)
TD.ROUTINES[:timstp] = () -> (STEPDT[] = TF.getscalar(:dtime); TD.fcall(:timstp_))
TF.init(deck; rundir = dir)
function snapnt()
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); I = a:b+1; Z = a:b
    v(s) = TF.var(s)
    (; t = TF.getscalar(:timen), dt = STEPDT[], jstep = Int(TF.getscalar(:jstep)), mu = Int(TF.getscalar(:mu)),
       r = v(:rn)[I], u = v(:un)[I], F = v(:frn)[I], rho = v(:dn)[Z], T = v(:tn)[Z], Er = v(:ern)[Z],
       p = v(:pgn)[Z], tetot = TF.getscalar(:tetot), tee = TF.getscalar(:tee), tes = TF.getscalar(:tes),
       tew = TF.getscalar(:tew), tel = TF.getscalar(:tel), teq = TF.getscalar(:teq))
end
S = [snapnt()]
t0 = time(); n = 0
while true
    global n += 1
    st = TD.titan_step!(); push!(S, snapnt()); st != 0 && break
end
wall = time() - t0
TF.finalize()
serialize(joinpath(dir, "snaps.jls"), S)
println("$name/$mode: steps = $n, t_end = $(S[end].t), wall = $(round(wall, digits = 2)) s")
