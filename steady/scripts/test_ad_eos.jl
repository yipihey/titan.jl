push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using ForwardDiff, Printf
module V1; include("../src/constants.jl"); include("../src/eos_v1.jl"); end
relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)
pts = [(1e-12, 6000.0), (4e-11, 6241.0), (1e-10, 9000.0), (3e-11, 1.2e4), (1e-10, 2e4), (1e-9, 5e4),
       (1e-8, 3e5), (6e-8, 8e5), (1e-15, 8000.0), (1e-7, 1e7)]
println("== 1. AD vs centred FD (h=1e-6 in ln) for P, e, s, ρ(P,T) ==")
worst = 0.0
for (ρ, T) in pts
    for (nm, f) in (("lnP", (a, b) -> log(Q.Ptot(exp(a), exp(b)))), ("e", (a, b) -> Q.eint(exp(a), exp(b))),
                    ("s", (a, b) -> Q.entropy(exp(a), exp(b))))
        g = ForwardDiff.gradient(v -> f(v[1], v[2]), [log(ρ), log(T)])
        hh = 1e-6
        fd = [(f(log(ρ) + hh, log(T)) - f(log(ρ) - hh, log(T))) / 2hh, (f(log(ρ), log(T) + hh) - f(log(ρ), log(T) - hh)) / 2hh]
        e = maximum(relerr.(g, fd)); global worst = max(worst, e)
        e > 1e-6 && @printf("  large: %s at ρ=%g T=%g: AD=%s FD=%s\n", nm, ρ, T, g, fd)
    end
    P = Q.Ptot(ρ, T)
    g = ForwardDiff.gradient(v -> log(Q.rho_PT(exp(v[1]), exp(v[2]))), [log(P), log(T)])
    hh = 1e-6
    fd = [(log(Q.rho_PT(P * exp(hh), T)) - log(Q.rho_PT(P * exp(-hh), T))) / 2hh,
          (log(Q.rho_PT(P, T * exp(hh))) - log(Q.rho_PT(P, T * exp(-hh)))) / 2hh]
    e = maximum(relerr.(g, fd)); global worst = max(worst, e)
    e > 1e-6 && @printf("  large: rho_PT at ρ=%g T=%g: AD=%s FD=%s (β=%.3g)\n", ρ, T, g, fd, 1 - Q.Prad(T) / P)
end
@printf("  worst relative AD–FD difference: %.2e  (FD truncation/round-off ~1e-8)\n", worst)

println("== 2. Thermodynamic consistency (exact identities) ==")
w1 = 0.0; w2 = 0.0
for (ρ, T) in pts
    st = Q.eos_rhoT(ρ, T)
    # (∂e/∂ρ)_T = (P − T (∂P/∂T)_ρ) / ρ²
    lhs = st.dedρ; rhs = (st.P - T * st.dPdT) / ρ^2
    global w1 = max(w1, abs(lhs - rhs) / (abs(st.P) / ρ^2))
    # (∂s/∂T)_ρ = (1/T)(∂e/∂T)_ρ
    dsdT = ForwardDiff.derivative(t -> Q.entropy(ρ, t), T)
    global w2 = max(w2, relerr(dsdT, st.cV / T))
end
@printf("  max |∂e/∂ρ − (P − T∂P/∂T)/ρ²| / (P/ρ²) = %.2e\n  max rel |∂s/∂T − c_V/T| = %.2e\n", w1, w2)

println("== 3. Eddington standard-model Γ1(β) ==")
for β in (0.5, 0.1, 0.01, 0.007)
    T = 1e7; P = Q.Prad(T) / (1 - β); st = Q.thermo(P, T)
    @printf("  β=%.3f: Γ1=%.7f analytic=%.7f\n", β, st.Γ1, (32 - 24β - 3β^2) / (24 - 21β))
end

println("== 4. Second derivatives (nested Duals): Hessian of ln P and e vs FD of AD gradient ==")
for (ρ, T) in pts[[2, 4, 6]]
    for (nm, f) in (("lnP", v -> log(Q.Ptot(exp(v[1]), exp(v[2])))), ("e", v -> Q.eint(exp(v[1]), exp(v[2]))))
        x = [log(ρ), log(T)]
        H = ForwardDiff.hessian(f, x)
        hh = 1e-5
        Hfd = hcat([(ForwardDiff.gradient(f, x .+ hh .* (1:2 .== k)) .- ForwardDiff.gradient(f, x .- hh .* (1:2 .== k))) ./ 2hh for k in 1:2]...)
        @printf("  %-3s ρ=%.0e T=%.0e: max |H−H_fd|/max|H| = %.1e\n", nm, ρ, T, maximum(abs.(H .- Hfd)) / maximum(abs.(H)))
    end
end

println("== 5. Regression vs v1 (finite-difference thermo) ==")
for (ρ, T) in pts
    P = Q.Ptot(ρ, T); a = Q.thermo(P, T); b = V1.thermo(P, T)
    @printf("  ρ=%.0e T=%.0e: ρ rel %.1e  ∇ad rel %.1e  Γ1 rel %.1e  cP rel %.1e\n", ρ, T,
            relerr(a.ρ, b.ρ), relerr(a.∇ad, b.∇ad), relerr(a.Γ1, b.Γ1), relerr(a.cP, b.cP))
end

println("== 6. Timing ==")
st = Q.eos_rhoT(4e-11, 6241.0)
t = @elapsed for _ in 1:2000; Q.eos_rhoT(4e-11, 6241.0); end
@printf("  eos_rhoT: %.1f μs/call (with AD first derivatives)\n", 1e6 * t / 2000)
t = @elapsed for _ in 1:2000; Q.saha(4e-11, 6241.0); end
@printf("  saha (Float64): %.1f μs/call\n", 1e6 * t / 2000)
