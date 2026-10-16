# usage: julia --project=env phase2_run.jl NAME NN MODE
#   MODE = fortran : Julia step driver, all-Fortran routines (scaffold check)
#          julia   : Julia driver + Julia-native eos/opac hooked in from t = 0
#          tierb   : Fortran trajectory; at checkpoints do tier-(a) routine checks
#                    and tier-(b) one-step comparisons (Fortran vs Julia eos/opac)
# Writes lib/p2_runs/NAME_MODE/{snaps.jls, log.txt} (+ the usual TITAN files).
using Serialization, Printf
include(joinpath(@__DIR__, "TitanDriver.jl"))
include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "Metrics.jl"))
using .TitanDriver, .JuliaPhysics, .Metrics
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics
TitanDriver.bitwise!(true)     # parity/byte-compare script: reproduce TITAN 1.2 bitwise

name, nn, mode = ARGS[1], parse(Int, ARGS[2]), ARGS[3]
dir = joinpath(TF.ROOT, "lib", "p2_runs", "$(name)_$mode")
logio = open(joinpath(@__DIR__, "..", "lib", "p2_runs", "$(name)_$mode.log"), "w")
say(x...) = (println(x...); println(logio, x...); flush(logio))

mode == "julia" && TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
TF.init(nn; rundir = dir)

# ---------------------------------------------------------------- tier (a)/(b)
const OUT_EOS  = (:pgn, :egn, :xnen, :delsn, :cpn, :betan, :gamman, :asn)
const OUT_OPAC = (:chifn, :xken, :xkpn, :xnen, :dlcfdldn, :dlcfdltn)
krange() = (Int(TF.getscalar(:ngrs)) - 1):(Int(TF.getscalar(:ngre)) + 2)
relmax(a, b) = maximum(abs(x - y) / max(abs(x), abs(y), floatmin()) for (x, y) in zip(a, b); init = 0.0)

function tier_a(label, ffort, fjl, outs)
    snap0 = TD.snapshot()
    ffort(); A = Dict(n => TF.var(n)[krange()] for n in outs)
    TD.restore!(snap0)
    fjl();   B = Dict(n => TF.var(n)[krange()] for n in outs)
    TD.restore!(snap0)
    worst = maximum(relmax(A[n], B[n]) for n in outs)
    say(@sprintf("  tier(a) %-5s %d outputs x %d zones: max rel diff = %.3e  %s", label,
                 length(outs), length(krange()), worst, worst <= 1e-12 ? "PASS (<=1e-12)" : "FAIL"))
    worst
end

const IV = (:rn, :xmen)           # interface, log variables
const ZV = (:dn, :tn, :ern)       # zone, log variables
function stepstate()
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    d = Dict{Symbol,Any}(n => TF.var(n)[a:b+1] for n in (IV..., :un, :frn))
    for n in ZV; d[n] = TF.var(n)[a:b]; end
    d[:iter] = Int(TF.getscalar(:iter)); d[:dtime] = TF.getscalar(:dtime); d[:timen] = TF.getscalar(:timen)
    d
end
function tier_b(nstep)
    snap0 = TD.snapshot()
    TD.step!(); A = stepstate()
    TD.restore!(snap0)
    TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
    TD.step!(); B = stepstate()
    TD.install_hooks!()
    TD.restore!(snap0)
    lg = maximum(relmax(A[n], B[n]) for n in (IV..., ZV...))
    nu = maximum(abs.(A[:un] .- B[:un])) / max(maximum(abs.(A[:un])), floatmin())
    nf = maximum(abs.(A[:frn] .- B[:frn])) / max(maximum(abs.(A[:frn])), floatmin())
    conv = TF.getscalar(:conv)
    ok = lg <= conv && nu <= conv && nf <= conv
    say(@sprintf("  tier(b) step %4d: iter F/J = %d/%d, max rel diff (r,m,ρ,T,E) = %.2e, |Δu|/max|u| = %.2e, |ΔF|/max|F| = %.2e, Δdt/dt = %.2e  %s",
                 nstep, A[:iter], B[:iter], lg, nu, nf, abs(B[:dtime] / A[:dtime] - 1),
                 ok ? "PASS (<= conv=$(conv))" : "FAIL"))
    (; lg, nu, nf, ok)
end

# ---------------------------------------------------------------- run
snaps = Metrics.Snap[]
checkpoints = Set([1, 2, 5, 10, 25, 50, 100, 200, 300, 400, 500, 600, 700])
t0 = time()
n = 0
res_b = []
while true
    global n += 1
    if mode == "tierb" && n in checkpoints
        say("checkpoint before step $n (jstep = $(TF.getscalar(:jstep)), t = $(TF.getscalar(:timen)))")
        tier_a("eos", TD.fortran_eos, JP.eos!, OUT_EOS)
        TF.getscalar(:lrad) == 1 && tier_a("opac", TD.fortran_opac, JP.opac!, OUT_OPAC)
        push!(res_b, tier_b(n))
    end
    st = TD.titan_step!()
    push!(snaps, Metrics.snap(TF))
    st != 0 && break
end
wall = time() - t0
s = snaps[end]
say("$name/$mode: steps = $n, final jstep = $(s.jstep), t = $(s.time), wall = $(round(wall, digits = 2)) s")
mode == "tierb" && say("tier(b) summary: $(count(r -> r.ok, res_b))/$(length(res_b)) checkpoints within conv")
serialize(joinpath(dir, "snaps.jls"), snaps)
TF.finalize()
close(logio)
