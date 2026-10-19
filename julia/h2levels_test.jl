# Tests of H2Levels (G6b): non-LTE H2 level populations.  julia --project=env h2levels_test.jl
using Printf
include(joinpath(@__DIR__, "PrimordialChem.jl")); include(joinpath(@__DIR__, "H2Lines.jl")); include(joinpath(@__DIR__, "H2Levels.jl"))
const PC = PrimordialChem; const HL = H2Lines; const HV = H2Levels
L = HV.load(); LL = HL.load()
println("H2Levels: $(length(L.E)) levels, $(length(L.A)) lines; collision coverage:"); HV.coverage(L)
ok = Bool[]
check(name, c) = (push!(ok, c); println(c ? "  PASS  " : "  FAIL  ", name))

# ---- GA08 low-density fits (Glover & Abel 2008 Tables 1-6, eqs. 27-34, separate ortho/para; combined at 3:1 with eqs. 30-32;
# erg cm³ s⁻¹ per H₂ per collider; coefficients: docs/popiii_rates_extraction.md §D.4) and the 3:1 fit of Table 8
poly(a, x) = sum(a[i] * x^(i - 1) for i in eachindex(a))
const GA = (
    Ho_lo = [-24.330855, 4.4404496, -4.0460989, -1.1390725, 9.8094223, 8.6273872], Ho_hi = [-24.329086, 4.6105087, -3.9505350, 12.363818, -32.403165, 48.853562, -38.542008, 12.066770],
    Hp_lo = [-24.216387, 3.3237480, -11.642384, -35.553366, -35.105689, -10.922078], Hp_hi = [-24.216387, 4.2046488, -1.3155285, -1.6552763, 4.1780102, -0.56949697, -3.3824407, 1.0904027],
    pp = [-23.889798, 1.8550774, -0.55593388, 0.28429361, -0.20581113, 0.13112378], po = [-23.748534, 1.76676480, -0.58634325, 0.31074159, -0.17455629, 0.18530758],
    op = [-24.126177, 2.3258217, -1.0082491, 0.54823768, -0.33679759, 0.20771406], oo = [-24.020047, 2.2687566, -1.0200304, 0.83561432, -0.40772247, 0.096025713],
    Hep = [-23.489029, 1.8210825, -0.59110559, 0.42280623, -0.30171138, 0.12872839], Heo = [-23.7749, 2.40654, -1.23449, 0.739874, -0.258940, 0.120573],
    Hpp = [-21.757160, 1.3998367, -0.37209530, 0.061554519, -0.37238286, 0.23314157], Hpo = [-21.706641, 1.3901283, -0.34993699, 0.075402398, -0.23170723, 0.068938876])
function ga08(coll, T; xo = 0.75)
    x = log10(T / 1000); xp = 1 - xo; f(a) = 10^poly(a, x)
    coll === :H  ? xo * f(T <= 1000 ? GA.Ho_lo : GA.Ho_hi) + xp * f(T <= 1000 ? GA.Hp_lo : GA.Hp_hi) :
    coll === :H2 ? xp * xp * f(GA.pp) + xp * xo * f(GA.po) + xo * xp * f(GA.op) + xo * xo * f(GA.oo) :
    coll === :He ? xo * f(GA.Heo) + xp * f(GA.Hep) :
                   xo * f(GA.Hpo) + xp * f(GA.Hpp)
end
function ga08_t8(coll, T)       # Table 8 (3:1 fits)
    x = log10(T / 1000)
    coll === :H  ? 10^(T <= 1000 ? poly([-24.311209, 3.5692468, -11.332860, -27.850082, -21.328264, -4.2519023], x) :
                                    poly([-24.311209, 4.6450521, -3.7209846, 5.9369081, -5.5108047, 1.5538288], x)) :
    coll === :H2 ? 10^poly([-23.962112, 2.09433740, -0.77151436, 0.43693353, -0.14913216, -0.033638326], x) :
    coll === :He ? 10^poly([-23.689237, 2.1892372, -0.81520438, 0.29036281, -0.16596184, 0.19191375], x) :
                   10^poly([-21.716699, 1.3865783, -0.37915285, 0.11453688, -0.23214154, 0.058538864], x)
end

# ---- a. LTE limit
println("\na. LTE limit (β = 1): max relative deviation from Boltzmann over the 63 levels")
for T in (300.0, 1000.0, 3000.0)
    q = HV.lte(L, T)
    # reference check against H2Lines.populations (partition function over all 302 levels), restricted and renormalised
    qH = HL.populations(LL, T); kk = [findfirst(x -> abs(x - L.E[i]) < 1e-16 * 1e3, LL.E) for i in eachindex(L.E)]
    for nH in (1e12, 1e14)
        p = HV.populations(L, T, nH, 0.0, 0.0, 0.0)
        m = q .> 1e-4; e1 = maximum(abs.(p[m] ./ q[m] .- 1))
        pa = HV.populations(L, T, nH, 0.0, 0.1nH, 0.0)         # He (all 63 levels linked) added: whole level set
        m2 = q .> 1e-7; e2 = maximum(abs.(pa[m2] ./ q[m2] .- 1))
        @printf("  T = %5.0f K n_H = %.0e: H only, levels with q > 1e-4 (%d): %.2e;  H + 0.1 He, q > 1e-7 (%d): %.2e;  Σ ortho = %.4f\n",
                T, nH, count(m), e1, count(m2), e2, sum(p[L.ortho]))
        check(@sprintf("LTE T=%g n_H=%g < 1%%", T, nH), e1 < 0.01 && e2 < 0.01)
    end
    # ExoMol 302-level Boltzmann restricted to the 63 levels vs ours (truncation of the partition function)
    sel = sort!([k for k in eachindex(LL.E) if any(abs(LL.E[k] - e) < 1e-25 for e in L.E)])
    @printf("  (63-level fraction of the 302-level partition function at %g K: %.5f)\n", T, sum(qH[sel]))
end
# radiation in equilibrium with the gas, no collisions: B_ν(T) mean intensity -> Boltzmann with the ortho:para ratio fixed
for T in (300.0, 1000.0)
    q = HV.lte(L, T); Jb = [HL.Bν(L.ν[i], T) for i in eachindex(L.A)]
    p = HV.populations(L, T, 0.0, 0.0, 0.0, 0.0; Jbar = Jb, opr = sum(q[L.ortho]) / sum(q[.!L.ortho]))
    m = q .> 1e-7; e = maximum(abs.(p[m] ./ q[m] .- 1))
    @printf("  J̄ = B_ν(T), no collisions, T = %g K: max |p/q − 1| = %.2e\n", T, e)
    check(@sprintf("Planck J̄ -> Boltzmann T=%g", T), e < 1e-6)
end

# ---- b. low-density limit vs GA08
println("\nb. low-density cooling Λ/(n_H2 n_coll) [erg cm³ s⁻¹] (fixed o:p = 3:1, n_coll = 1e-4 cm⁻³) vs Glover & Abel 2008 fits (Tables 1-6 at 3:1)")
@printf("  %6s | %-22s | %-22s | %-22s | %-22s\n", "T [K]", "H: ours  GA08  ratio", "H2: ours  GA08  ratio", "He: ours  GA08  ratio", "H+: ours  GA08  ratio")
rat = Dict(c => Float64[] for c in (:H, :H2, :He, :Hp))
for T in (100.0, 200.0, 300.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 6000.0)
    row = Float64[]
    for (c, args) in ((:H, (1e-4, 0.0, 0.0, 0.0)), (:H2, (0.0, 1e-4, 0.0, 0.0)), (:He, (0.0, 0.0, 1e-4, 0.0)), (:Hp, (0.0, 0.0, 0.0, 1e-4)))
        p = HV.populations(L, T, args...; opr = 3.0); lam = HV.cooling(L, p, 1.0) / 1e-4
        push!(row, lam, ga08(c, T), lam / ga08(c, T)); push!(rat[c], lam / ga08(c, T))
    end
    @printf("  %6.0f | %.2e %.2e %5.2f | %.2e %.2e %5.2f | %.2e %.2e %5.2f | %.2e %.2e %5.2f   (Table 8: H %5.2f)\n", T, row...,
            row[1] / ga08_t8(:H, T))
end
Ts = [100.0, 200.0, 300.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 6000.0]
for c in (:H, :H2, :He, :Hp)
    r = rat[c]; @printf("  %-3s ratio range %.2f – %.2f; median %.2f; (T ≤ 2000 K: %.2f – %.2f)\n", c, minimum(r), maximum(r),
                        sort(r)[6], minimum(r[Ts .<= 2000]), maximum(r[Ts .<= 2000]))
end
println("  note: H at 200-500 K differs from the GA08 (Wrathmall & Flower 2007) fits by up to a factor 1.8 (0.4 at 100 K); H+ is 15-87 % of GA08 (Gerlich table: v = 0, J <= 9, no o<->p in this comparison)")
check("H ratio within 0.75–1.35 for 700–6000 K", all(0.75 .< rat[:H][Ts .>= 700] .< 1.35))
check("H2, He ratio within 0.8–1.25 for 100–6000 K", all(0.8 .< rat[:H2] .< 1.25) && all(0.8 .< rat[:He] .< 1.25))

# ---- c. Λ(n_H)/Λ_LTE vs the GP98/HM79 combination of PrimordialChem.cool_H2_thin
println("\nc. Λ(n_H)/Λ_LTE, H colliders only (free o:p via reactive H), vs 1/(1 + Λ_LTE^HM/(n_H Λ0^GP)) of cool_H2_thin")
rows = Tuple{Float64,Float64,Float64,Float64}[]
for T in (300.0, 1000.0, 3000.0)
    q = HV.lte(L, T); Llte = HV.cooling(L, q, 1.0); Lhm = PC.Λ_H2_LTE_HM(T); L0 = PC.Λ_H2_0_GP(T)
    @printf("  T = %g K: Λ_LTE(63 levels) = %.3e, HM79 = %.3e, Λ0^GP = %.3e, Λ0^ours (n_H=1e-3)/n = %.3e, n_crit^GP = Λ_LTE/Λ0 = %.2e, n_crit^ours (ratio 0.5) = ",
            T, Llte, Lhm, L0, HV.cooling(L, HV.populations(L, T, 1e-3, 0.0, 0.0, 0.0), 1.0) / 1e-3, Lhm / L0)
    nh = 10.0 .^ (0:0.5:12); rs = Float64[]
    for n in nh
        p = HV.populations(L, T, n, 0.0, 0.0, 0.0); push!(rs, HV.cooling(L, p, 1.0) / Llte)
    end
    i = findfirst(>(0.5), rs); ncr = 10^(log10(nh[i-1]) + (0.5 - rs[i-1]) / (rs[i] - rs[i-1]) * 0.5)
    @printf("%.2e\n", ncr)
    @printf("    n_H     : %s\n", join([@sprintf("%8.0e", n) for n in nh[1:2:end]], " "))
    @printf("    ours    : %s\n", join([@sprintf("%8.3f", r) for r in rs[1:2:end]], " "))
    @printf("    GP/HM   : %s\n", join([@sprintf("%8.3f", 1 / (1 + Lhm / (n * L0))) for n in nh[1:2:end]], " "))
    @printf("    ours/GP : %s\n", join([@sprintf("%8.2f", r * (1 + Lhm / (n * L0))) for (n, r) in zip(nh[1:2:end], rs[1:2:end])], " "))
    check(@sprintf("T=%g: Λ -> Λ_LTE at n_H = 1e12 (within 2%%)", T), abs(rs[end-0] - 1) < 0.02)
end

# ---- d. trapping
println("\nd. trapping: scale all A by β (escape probability) at T = 1000 K, H colliders; n_eff = n_H where Λ = ½ Λ_LTE(β)")
T = 1000.0; q = HV.lte(L, T)
function ncrit(β)
    Ll = HV.cooling(L, q, 1.0; β = β)
    f(n) = HV.cooling(L, HV.populations(L, T, n, 0.0, 0.0, 0.0; β = β), 1.0; β = β) / Ll - 0.5
    lo, hi = -3.0, 14.0
    for _ in 1:60; mid = (lo + hi) / 2; f(10^mid) > 0 ? (hi = mid) : (lo = mid); end
    10^((lo + hi) / 2)
end
n1 = ncrit(1.0); nc = Float64[]
for β in (1.0, 0.3, 0.1, 0.03, 0.01, 1e-3)
    n = ncrit(β); push!(nc, n)
    p = HV.populations(L, T, 1e6, 0.0, 0.0, 0.0; β = β); lr = HV.cooling(L, p, 1.0; β = β) / HV.cooling(L, q, 1.0; β = β)
    @printf("  β = %-6g n_eff = %.3e  n_eff/n_eff(1) = %.4f  (β: %.4f)   Λ/Λ_LTE(β) at n_H = 1e6: %.3f\n", β, n, n / n1, β, lr)
end
βs = [1.0, 0.3, 0.1, 0.03, 0.01, 1e-3]
sl = (log(nc[end]) - log(nc[1])) / (log(βs[end]) - log(βs[1]))
@printf("  log-log slope of n_eff vs β over 1–1e-3: %.3f\n", sl)
check("n_eff decreases monotonically with β", issorted(nc; rev = true))
check("n_eff ∝ β^s with s in 0.8–1.05", 0.8 < sl < 1.05)
println(all(ok) ? "\nALL PASS ($(length(ok)) checks)" : "\nFAILED: $(count(!, ok)) of $(length(ok))")
