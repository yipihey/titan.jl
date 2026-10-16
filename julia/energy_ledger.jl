# Independent energy budget for a TITAN run (phase 3, energy investigation).
# usage: julia energy_ledger.jl NAME NN TAG [key=value ...]
#   keys (applied to COMMON after init, before step 1): stol conv thet lgrid niter
#   ncell: set env TITAN_NCELL before launching.
# Budget (all per 4π, as TITAN stores masses/volumes):
#   E   = Σ_zones m_k (e_g + E_r/ρ)  +  KE
#   KE_z = Σ m_k (u_k² + u_{k+1}²)/4          (TITAN's tee definition)
#   KE_i = Σ_interfaces ½ m̄_k u_k²            (interface-mass definition)
#   W   = ∫ [r²u(p + fE + q)]_in − [ ... ]_out dt        (work done ON the domain)
#   L   = ∫ [r²F]_out − [r²F]_in dt                      (radiation leaving)
#   residual = ΔE − W + L ; integrals use x_w = w·x_new + (1−w)·x_old, w ∈ {½, θ}
using Serialization, Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver
const TD = TitanDriver; const TF = TD.TF

name, nn, tag = ARGS[1], parse(Int, ARGS[2]), ARGS[3]
over = Dict(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[4:end])
outdir = joinpath(TF.ROOT, "lib", "p3_runs"); mkpath(outdir)
TF.init(nn; rundir = joinpath(outdir, "run_$(name)_$tag"))
for (k, v) in over
    T = eltype(TF.var(k)); TF.setscalar!(k, T <: Integer ? T(round(Int, v)) : v)
end
θ = TF.getscalar(:thet)
# dt actually used by each step (TITAN's clock does not advance on divergence
# retries: reset.f sets timen = timeo after step.f set timen = timeo + dtime)
const STEPDT = Ref(0.0)
TD.ROUTINES[:timstp] = () -> (STEPDT[] = TF.getscalar(:dtime); TD.fcall(:timstp_))

v(s) = TF.var(s)
function state()
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    rn, dn, un, egn, ern, dvoln, pgn, fedd, frn, qf = (v(s) for s in
        (:rn, :dn, :un, :egn, :ern, :dvoln, :pgn, :fedd, :frn, :qf))
    m = [dn[k] * dvoln[k] for k in a:b]
    Eint_g = sum(m[i] * egn[k] for (i, k) in enumerate(a:b))
    Eint_r = sum(ern[k] * dvoln[k] for k in a:b)
    KEz = sum(m[i] * 0.25 * (un[k]^2 + un[k+1]^2) for (i, k) in enumerate(a:b))
    mi = [0.5 * ((k > a ? m[k-a] : 0.0) + (k <= b ? m[k-a+1] : 0.0)) for k in a:b+1]
    KEi = sum(0.5 * mi[i] * un[k]^2 for (i, k) in enumerate(a:b+1))
    r2(k) = rn[k]^2 * (TF.getscalar(:mu) == 2 ? 1.0 : 0.0) + (TF.getscalar(:mu) == 0 ? 1.0 : 0.0)
    (; t = TF.getscalar(:timen), dt = STEPDT[], jstep = Int(TF.getscalar(:jstep)), iter = Int(TF.getscalar(:iter)),
       Eg = Eint_g, Er = Eint_r, KEz, KEi,
       Win_p = r2(a) * un[a] * pgn[a], Win_r = r2(a) * un[a] * fedd[a] * ern[a], Win_q = r2(a) * un[a] * qf[a],
       Wout_p = r2(b+1) * un[b+1] * pgn[b], Wout_r = r2(b+1) * un[b+1] * fedd[b] * ern[b],
       Wout_q = r2(b+1) * un[b+1] * qf[b],
       Lin = r2(a) * frn[a], Lout = r2(b+1) * frn[b+1],
       tetot = TF.getscalar(:tetot), tee = TF.getscalar(:tee), tew = TF.getscalar(:tew),
       tel = TF.getscalar(:tel), teq = TF.getscalar(:tel) * 0 + TF.getscalar(:teq), mass = sum(m), rshock = 0.0,
       a_rn = rn[a:b+1], a_un = un[a:b+1], a_m = m, a_qf = qf[a:b], a_dudr = v(:dudr)[a:b],
       a_r1 = v(:r1)[a:b+1], a_r3 = v(:r3)[a:b])
end

S = [state()]                 # initial model
nsteps = 0
while true
    st = TD.titan_step!(); global nsteps += 1
    push!(S, state())
    st != 0 && break
end
TF.finalize()

# ledger
function ledger(S, w; q = +1)
    n = length(S); W = zeros(n); L = zeros(n)
    for i in 2:n
        dt = S[i].dt
        x(f) = w * getfield(S[i], f) + (1 - w) * getfield(S[i-1], f)
        W[i] = W[i-1] + dt * ((x(:Win_p) + x(:Win_r) + q * x(:Win_q)) - (x(:Wout_p) + x(:Wout_r) + q * x(:Wout_q)))
        L[i] = L[i-1] + dt * (x(:Lout) - x(:Lin))
    end
    W, L
end
E(s, ke) = s.Eg + s.Er + getfield(s, ke)
scale = max(maximum(abs, ledger(S, 0.5)[1]), E(S[1], :KEz))   # total piston work (or E0 if closed)
say(io, x...) = (println(x...); println(io, x...))
open(joinpath(outdir, "energy_$(name)_$tag.log"), "w") do io
    say(io, "$name [$tag] overrides=$(over) ncell=$(get(ENV, "TITAN_NCELL", "default")) θ=$θ conv=$(TF.getscalar(:conv)) stol=$(TF.getscalar(:stol)) lgrid=$(TF.getscalar(:lgrid))")
    say(io, @sprintf("  steps=%d  t_end=%.5e  Σdt=%.5e (clock lost %.4e s in %d retried steps)  mass drift=%.2e  scale=%.4e",
                     nsteps, S[end].t, sum(s.dt for s in S[2:end]), sum(s.dt for s in S[2:end]) - S[end].t,
                     count(i -> abs(S[i].t - S[i-1].t - S[i].dt) > 1e-9 * S[i].dt, 2:length(S)),
                     S[end].mass / S[1].mass - 1, scale))
    for (lbl, w, q, ke) in (("midpoint, +q, KE_zone", 0.5, +1, :KEz), ("θ-weighted, +q, KE_zone", θ, +1, :KEz),
                            ("θ-weighted, −q, KE_zone", θ, -1, :KEz), ("θ-weighted, +q, KE_iface", θ, +1, :KEi))
        W, L = ledger(S, w; q)
        res = [E(S[i], ke) - E(S[1], ke) - W[i] + L[i] for i in eachindex(S)]
        say(io, @sprintf("  %-26s residual/W at end = %+.4e   max |res|/W = %.4e", lbl, res[end] / scale,
                         maximum(abs, res) / scale))
    end
    # component comparison with TITAN (both from step 1, where TITAN's ledger starts)
    W, L = ledger(S, 0.5)
    i1 = 2
    say(io, @sprintf("  components from step 1 to end, /scale:  ΔE_mine %+.5e  ΔE_titan(tee) %+.5e | W_mine %+.5e  -tew %+.5e | L_mine %+.5e  tel %+.5e | teq %+.3e",
        (E(S[end], :KEz) - E(S[i1], :KEz)) / scale, (S[end].tee - S[i1].tee) / scale,
        (W[end] - W[i1]) / scale, -(S[end].tew - S[i1].tew) / scale,
        (L[end] - L[i1]) / scale, (S[end].tel - S[i1].tel) / scale, (S[end].teq - S[i1].teq) / scale))
    # TITAN's own diagnostic for comparison (relative to its first value)
    e_t = (S[end].tetot - S[2].tetot) / scale
    say(io, @sprintf("  TITAN tetot drift (from step 1) / W = %+.4e", e_t))
    # temporal localisation of the θ-weighted residual
    W, L = ledger(S, θ)
    res = [E(S[i], :KEz) - E(S[1], :KEz) - W[i] + L[i] for i in eachindex(S)] ./ scale
    inc = diff(res)
    ord = sortperm(abs.(inc), rev = true)
    say(io, "  largest per-step residual increments (step, t, Δres/W, dt):")
    for i in ord[1:6]
        say(io, @sprintf("     step %4d  t=%.3e  %+.3e  dt=%.2e", i, S[i+1].t, inc[i], S[i+1].dt))
    end
    for f in (0.01, 0.05, 0.2, 0.5, 1.0)
        i = clamp(round(Int, f * length(res)), 1, length(res))
        say(io, @sprintf("     after %3d%% of steps (t=%.3e): residual/W = %+.4e", round(Int, 100f), S[i].t, res[i]))
    end
end
serialize(joinpath(outdir, "ledger_$(name)_$tag.jls"), S)
