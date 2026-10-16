push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using ForwardDiff, Printf
relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)
pts = [(1e-12, 6000.0), (4e-11, 6241.0), (1e-10, 9000.0), (3e-11, 1.2e4), (1e-9, 5e4), (1e-8, 3e5), (1e-15, 8000.0), (1e-7, 1e7)]
println("== fast eos_rhoT vs reference (3 gradients) ==")
w = 0.0
for (ρ, T) in pts
    a = Q.eos_rhoT(ρ, T); b = Q.eos_rhoT_ref(ρ, T)
    for k in (:P, :e, :s, :dPdρ, :dPdT, :dedρ, :cV, :Γ1, :∇ad)
        global w = max(w, relerr(getfield(a, k), getfield(b, k)))
    end
end
@printf("  max relative difference over 9 fields × %d points: %.1e\n", length(pts), w)
println("== nested AD: d(Γ1)/dT and d(c_V)/dρ through eos_rhoT vs FD ==")
for (ρ, T) in pts[[2, 4, 5]]
    d1 = ForwardDiff.derivative(t -> Q.eos_rhoT(ρ, t).Γ1, T)
    hh = 1e-6T; f1 = (Q.eos_rhoT(ρ, T + hh).Γ1 - Q.eos_rhoT(ρ, T - hh).Γ1) / 2hh
    d2 = ForwardDiff.derivative(r -> Q.eos_rhoT(r, T).cV, ρ)
    hr = 1e-6ρ; f2 = (Q.eos_rhoT(ρ + hr, T).cV - Q.eos_rhoT(ρ - hr, T).cV) / 2hr
    @printf("  ρ=%.0e T=%.1e: dΓ1/dT AD=%.6e FD=%.6e | dcV/dρ AD=%.6e FD=%.6e\n", ρ, T, d1, f1, d2, f2)
end
println("== timing ==")
Q.eos_rhoT(4e-11, 6241.0); Q.eos_rhoT_ref(4e-11, 6241.0)
t1 = @elapsed for _ in 1:5000; Q.eos_rhoT(4e-11, 6241.0); end
t2 = @elapsed for _ in 1:5000; Q.eos_rhoT_ref(4e-11, 6241.0); end
@printf("  eos_rhoT fast: %.1f μs   reference: %.1f μs\n", 1e6t1/5000, 1e6t2/5000)
println("== thermodynamic identities with fast path ==")
w1 = 0.0
for (ρ, T) in pts
    st = Q.eos_rhoT(ρ, T)
    global w1 = max(w1, abs(st.dedρ - (st.P - T * st.dPdT) / ρ^2) / (st.P / ρ^2))
end
@printf("  max |∂e/∂ρ − (P − T∂P/∂T)/ρ²|/(P/ρ²) = %.1e\n", w1)
