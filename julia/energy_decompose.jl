# Decompose the independent energy residual of a ledger run into
#   (i)  the θ-scheme kinetic-energy term   −(θ−½) Σ m̄ (u_new − u_old)²   per step
#   (ii) the artificial-viscosity mismatch   ∫ [heat − work] dt with TITAN's discrete forms
#        heat = −Σ_zones qf·dudr·dV_θ ,  work = Σ_interior u_θ r1 (r3_k qf_k − r3_{k−1} qf_{k−1})
#        mismatch = heat − work − (piston viscous work); the's viscous work u_θ r1 r3 qf at the inner boundary
# usage: julia energy_decompose.jl NAME TAG θ
using Serialization, Printf
name, tag, θ = ARGS[1], ARGS[2], parse(Float64, ARGS[3])
S = deserialize(joinpath(@__DIR__, "..", "lib", "p3_runs", "ledger_$(name)_$tag.jls"))
n = length(S)
W = zeros(n); L = zeros(n); Tt = zeros(n); Dv = zeros(n)
E(s) = s.Eg + s.Er + s.KEz
for i in 2:n
    s, o = S[i], S[i-1]; dt = s.dt
    x(f) = θ * getfield(s, f) + (1 - θ) * getfield(o, f)
    W[i] = W[i-1] + dt * ((x(:Win_p) + x(:Win_r) + x(:Win_q)) - (x(:Wout_p) + x(:Wout_r) + x(:Wout_q)))
    L[i] = L[i-1] + dt * (x(:Lout) - x(:Lin))
    if length(s.a_un) == length(o.a_un)
        m = s.a_m; nz = length(m)
        mbar = [0.5 * ((k > 1 ? m[k-1] : 0.0) + (k <= nz ? m[k] : 0.0)) for k in 1:nz+1]
        Tt[i] = Tt[i-1] - (θ - 0.5) * sum(mbar .* (s.a_un .- o.a_un) .^ 2)
        uθ = θ .* s.a_un .+ (1 - θ) .* o.a_un
        rθ = θ .* s.a_rn .+ (1 - θ) .* o.a_rn
        mu = abs(rθ[end] - rθ[1]) > 0 && s.a_r1[1] != 1.0 ? 2 : 0
        dV = [(rθ[k+1]^(mu+1) - rθ[k]^(mu+1)) / (mu + 1) for k in 1:nz]
        heat = -sum(s.a_qf .* s.a_dudr .* dV)
        work = sum(uθ[k] * s.a_r1[k] * (s.a_r3[k] * s.a_qf[k] - s.a_r3[k-1] * s.a_qf[k-1]) for k in 2:nz)
        bnd = uθ[1] * s.a_r1[1] * s.a_r3[1] * s.a_qf[1]
        Dv[i] = Dv[i-1] + dt * (heat - work - bnd)
    else
        Tt[i] = Tt[i-1]; Dv[i] = Dv[i-1]
    end
end
res = [E(S[i]) - E(S[1]) - W[i] + L[i] for i in 1:n]
sc = max(maximum(abs, W), E(S[1]))
@printf("%s [%s] θ=%.2f: residual %+.3e | θ-scheme KE term %+.3e | viscosity mismatch %+.3e | remainder %+.3e   (all / %.3e)\n",
        name, tag, θ, res[end] / sc, Tt[end] / sc, Dv[end] / sc, (res[end] - Tt[end] - Dv[end]) / sc, sc)
