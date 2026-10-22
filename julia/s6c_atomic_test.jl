# S6c acceptance tests of the atomic H/He chemistry and cooling (ChemNet He species, AtomicHHe). Plain Julia, no TITAN:
#   julia --project=env s6c_atomic_test.jl            (one process, ≈ 2–3 min)
# Output: PASS/FAIL per check with the numbers; last line "S6C ATOMIC: ALL PASS" / "S6C ATOMIC: n FAIL".
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl"))
const CN = ChemNet; const A = AtomicHHe; const yHe = CN.yHe
const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-34s %s\n", ok ? "PASS" : "FAIL", name, detail)
end
Bν(ν, T) = 2CN.hP * ν^3 / CN.c^2 / expm1(CN.hP * ν / (CN.kB * T))
envz = (Trad = 0.0, J = ν -> 0.0, static = true)

# ====================================================================== a. detailed balance: J = B(T) gives Saha
"Saha composition (independent of the network's rates): H, H⁺, H⁻, H₂, H₂⁺, He, He⁺, He²⁺, Li, Li⁺, e at (n_H, T) from ChemNet.lnKeq, charge balance by bisection."
function saha(nH, T)
    K(R, P) = CN.lnKeq(R, P, T)
    lKh = K([:H], [:Hp, :e]); lK1 = K([:He], [:Hep, :e]); lK2 = K([:Hep], [:Hepp, :e]); lKL = K([:Li], [:Lip, :e])
    lKm = K([:Hm], [:H, :e]); lKH2 = K([:H2], [:H, :H]); lKp = K([:H2p], [:H, :Hp])
    function comp(ne)
        rH = exp(lKh - log(ne)); r1 = exp(lK1 - log(ne)); r2 = exp(lK2 - log(ne)); rL = exp(lKL - log(ne))
        a = 2nH * exp(-lKH2) + 2nH * rH * exp(-lKp); b = 1 + rH + ne * exp(-lKm)
        xH = 2 / (b + sqrt(b^2 + 4a))
        xHp = xH * rH; xHm = xH * ne * exp(-lKm); xH2 = nH * xH^2 * exp(-lKH2); xH2p = nH * xH^2 * rH * exp(-lKp)
        xHe = yHe / (1 + r1 + r1 * r2); xHep = xHe * r1; xHepp = xHep * r2
        xLi = CN.yLi / (1 + rL); xLip = xLi * rL
        Dict(:H => xH, :Hp => xHp, :Hm => xHm, :H2 => xH2, :H2p => xH2p, :He => xHe, :Hep => xHep, :Hepp => xHepp, :Li => xLi, :Lip => xLip, :e => ne / nH),
        xHp + xH2p + xHep + 2xHepp + xLip - xHm
    end
    lo, hi = 1e-300, 1e30
    for _ in 1:400; m = sqrt(lo * hi); d, q = comp(m); (q * nH > m) ? (lo = m) : (hi = m); end
    comp(sqrt(lo * hi))[1]
end
# Network with the inverses of its one-way radiative processes (valid in a thermal photon bath only, see AtomicHHe.atomic)
netB, phB = A.atomic(li = true, jb_inverse = true)
"Largest |forward − reverse|/max flux over all reaction pairs and photo pairs at the composition d."
function fluxcheck(net, ph, d, nH, T, env)
    rx = net.reactions; flux = [(rate = r.k(T, env) * nH^(length(r.R) - 1); for s in r.R; rate *= d[s]; end; rate) for r in rx]
    worst = 0.0; wn = ""; unpaired = String[]
    idx = Dict((sort(r.R), sort(r.P)) => i for (i, r) in enumerate(rx))
    for (i, r) in enumerate(rx)
        j = get(idx, (sort(r.P), sort(r.R)), 0)
        if j == 0; push!(unpaired, r.name); continue; end
        j < i && continue
        f, b = flux[i], flux[j]; m = max(f, b); m < 1e-250 && continue
        e = abs(f - b) / m; e > worst && (worst = e; wn = r.name)
    end
    for p in ph
        R, α = CN.photorates(p, T, env); α *= p.esc(T, env)
        fw = R * d[p.R]; bw = α * nH^(length(p.P) - 1); for s in p.P; bw *= d[s]; end
        m = max(fw, bw); m < 1e-250 && continue
        e = abs(fw - bw) / m; e > worst && (worst = e; wn = p.name)
    end
    worst, wn, unpaired
end
let w0 = 0.0, wn0 = "", np = 0, nun = 0, count = 0
    for T in (1e3, 2e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6), nH in (1e0, 1e4, 1e8, 1e12, 1e16, 1e20)
        env = (Trad = T, J = ν -> Bν(ν, T), static = false)
        w, n, up = fluxcheck(netB, phB, saha(nH, T), nH, T, env)
        count += 1; nun += length(up)
        w > w0 && (w0 = w; wn0 = n)
    end
    check("a1 flux balance at Saha", w0 <= 1e-6 && nun == 0, @sprintf("max |net|/flux = %.1e over %d (T = 1e3–1e6 K, n_H = 1e0–1e20) states, all %d reactions + %d photo pairs paired (worst: %s)", w0, count, length(netB.reactions), length(phB), wn0))
end
let worst = 0.0, ws = :none, wcase = (0.0, 0.0), count = 0, nfail = 0
    for T in (3e3, 1e4, 3e4, 1e5, 3e5, 1e6), nH in (1e4, 1e12, 1e20)
        env = (Trad = T, J = ν -> Bν(ν, T), static = false)
        ex = saha(nH, T)
        x0 = [ex[s] * (iseven(i) ? 0.97 : 1.03) for (i, s) in enumerate(netB.free)]           # then caps so that H, He, Li stay positive
        for (i, s) in enumerate(netB.free)
            s === :Hp && (x0[i] = min(x0[i], 1 - 1e-9 - 2 * ex[:H2] - ex[:H2p] - ex[:Hm]))
            s === :Hepp && (x0[i] = min(x0[i], yHe * 0.999999))
            s === :Hep && (x0[i] = min(x0[i], yHe * 0.999999 - x0[findfirst(==(:Hepp), netB.free)]))
            s === :Lip && (x0[i] = min(x0[i], CN.yLi * 0.999999))
        end
        x = try CN.steady(x0, nH * CN.μH(), T, netB; env, photos = phB, xfloor = 1e-14, tol = 1e-11, xclamp = 1e-14) catch e; (nfail += 1; println("   a2: no convergence at T = $T, n_H = $nH"); continue) end
        d = CN.full(x, netB); count += 1
        for s in (:H, :Hp, :Hm, :H2, :H2p, :He, :Hep, :Hepp, :Li, :Lip)
            thr = 1e-6       # eliminated species (1 − Σ…) are resolved to 1e-16 absolute
            ex[s] < thr && continue
            r = abs(d[s] / ex[s] - 1); r > worst && (worst = r; ws = s; wcase = (T, nH))
        end
    end
    check("a2 relaxation to Saha", worst <= 5e-6 && nfail == 0, @sprintf("pseudo-transient steady state from a Saha state with every free species scaled by 0.97 or 1.03, vs Saha for the species with x ≥ 1e-6 (smaller ones are fixed by a1, the pseudo-time stop rule limits this to ~1e-6; tolerance 5e-6): max rel. dev. %.1e (%s at T = %.0e, n_H = %.0e), %d states, T = 3e3–1e6 K, n_H = 1e4–1e20 (1e3 K: rates ≲ 1e-30 s⁻¹, not relaxable; covered by a1)", worst, ws, wcase[1], wcase[2], count))
end
let Tt = 2e4, K = exp(CN.lnKeq([:He], [:Hep, :e], Tt))
    # independent hand evaluation of the He Saha constant: (2π m_e kT/h²)^{3/2} (m_He+/m_He)^{3/2} g_e Z(He+)/Z(He) e^{-χ/kT} with Z from explicit level sums
    lev = [(0.0, 1), (19.8196, 3), (20.6158, 1), (20.9642, 9), (21.218, 3)]; Zhe = sum(g * exp(-E * CN.eV / (CN.kB * Tt)) for (E, g) in lev) +
          sum(let n = n; E = (CN.χHe1 - CN.Ry_eV * CN.eV / n^2) / CN.eV; (g = n <= 3 ? 36.0 : 4.0n^2; g * exp(-E * CN.eV / (CN.kB * Tt))) end for n in 3:CN.NH_LEVELS)
    Zhep = 2 + 8exp(-40.81 * CN.eV / (CN.kB * Tt))
    Kh = (2π * CN.me * CN.kB * Tt / CN.hP^2)^1.5 * 2 * Zhep / Zhe * (1 - CN.me / (4.0026 * 1.66053907e-24))^1.5 * exp(-CN.χHe1 / (CN.kB * Tt))
    check("a3 He Saha constant", abs(K / Kh - 1) < 5e-3, @sprintf("Keq(He → He⁺ + e) at 2e4 K = %.4e cm⁻³ vs hand evaluation %.4e (ratio %.5f; hand sum drops n ≥ 3 term detail and n_max shell of He⁺)", K, Kh, K / Kh))
end

# ====================================================================== rate comparisons (information + radiative recombination totals)
println("\nRate comparisons (ours / reference):")
@printf("  %9s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s\n", "T [K]", "αA(H)/HG97", "αA(He)/VF96", "αA(He+)/VF96", "Voron/Jan H", "He", "He+", "αB(H)/HG97", "αB(He+)/GA08", "")
env0 = (Trad = 1.0, J = ν -> 0.0, static = true)
net, ph = A.atomic()
αtot(RR) = (T -> sum(CN.photorates(p, T, env0)[2] for p in ph if p.R === RR))
αexc(RR) = (T -> sum(CN.photorates(p, T, env0)[2] for p in ph if p.R === RR && endswith(p.name, "excited")))
hg_A(T) = 1.269e-13 * (315614 / T)^1.503 * (1 + (604625 / T)^0.470)^-1.923
hg_B(T) = 2.753e-14 * (315614 / T)^1.5 * (1 + (115188 / T)^0.407)^-2.242
ga_B2(T) = 5.506e-14 * (1262456 / T)^1.5 * (1 + (460752 / T)^0.407)^-2.242
let eH = 0.0, eHe = 0.0, eHe2 = 0.0, eB = 0.0, eB2 = 0.0
    for T in (3e3, 5e3, 1e4, 3e4, 1e5, 3e5, 1e6)
        aH = αtot(:H)(T) / hg_A(T); aHe = αtot(:He)(T) / A.vf96_HeI(T); aHe2 = αtot(:Hep)(T) / A.vf96_HeII(T)
        bH = αexc(:H)(T) / hg_B(T); bHe2 = αexc(:Hep)(T) / ga_B2(T)
        @printf("  %9.0f | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f\n", T, aH, aHe, aHe2, A.voronov_H(T) / PrimordialChem.k1_aazn(T), A.voronov_He0(T) / A.janev_He0(T), A.voronov_Hep(T) / A.janev_Hep(T), bH, bHe2)
        if T >= 1e4; eH = max(eH, abs(aH - 1)); eHe2 = max(eHe2, abs(aHe2 - 1)); eB = max(eB, abs(bH - 1)); eB2 = max(eB2, abs(bHe2 - 1)); end
        if 3e3 <= T <= 3e4; eHe = max(eHe, abs(aHe - 1)); end
    end
    check("c2 radiative recombination", eH < 0.06 && eHe < 0.06 && eHe2 < 0.15 && eB < 0.06 && eB2 < 0.17,
          @sprintf("Milne totals vs Hui & Gnedin 1997 / Verner & Ferland 1996 (n ≤ 15 levels, which costs 5 %% (H) and 12 %% (He²⁺) at 1e4 K): H case A ≤ %.1f %% (1e4–1e6 K), He⁺→He ≤ %.1f %% (3e3–3e4 K; excited-level σ scale calibrated), He²⁺→He⁺ ≤ %.1f %%; case B (excited levels only): H ≤ %.1f %%, He²⁺ ≤ %.1f %% (GA08 k20,B)", 100eH, 100eHe, 100eHe2, 100eB, 100eB2))
end
let T = 1.0e4
    for (RR, name, χ) in ((:H, "H", CN.χH), (:He, "He", CN.χHe1), (:Hep, "He+", CN.χHe2))
        E = sum(CN.photoenergy(p, T, env0)[2] + CN.photoenergy(p, T, env0)[3] for p in ph if p.R === RR); a = αtot(RR)(T)
        kT = CN.kB * T
        @printf("  energy radiated per recombination to %-3s at 1e4 K: %.3f eV = chi + %.2f kT (chi %.3f eV)\n", name, E / a / CN.eV, (E / a - χ) / kT, χ / CN.eV)
    end
    Ep = [(let E = sum(CN.photoenergy(p, T, env0)[2] + CN.photoenergy(p, T, env0)[3] for p in ph if p.R === RR); (E / αtot(RR)(T) - χ) / (CN.kB * T) end) for (RR, χ) in ((:H, CN.χH), (:He, CN.χHe1), (:Hep, CN.χHe2))]
    check("c3 recombination energy", all(0.2 .< Ep .< 1.5), @sprintf("energy radiated per recombination = χ + (%.2f, %.2f, %.2f) kT for H, He, He⁺ at 1e4 K (kinetic part 0.2–1.5 kT; χ counted once: the explicit loss is the whole photon energy, the ionisation energy leaves the gas through E0 only at the ionisation side)", Ep...))
end

let Trad = 2.725 * 20, T = 3000.0, dd = Dict(:H => 0.0, :Hp => 1.0, :He => yHe, :Hep => 0.0, :Hepp => 0.0, :e => 1.0, :Hm => 0.0), nH = 1.0
    envc = (Trad = Trad, J = ν -> Bν(ν, Trad), static = false)
    Λc = A.parts(nH * CN.μH(), T, dd; photos = CN.Photo[], env = envc).comp / (nH * 1.0)
    cen = 5.65e-36 * 20^4 * (T - 2.7 * 20)                   # Cen 1992 (z = 19)
    check("c4 Compton", abs(Λc / cen - 1) < 0.02, @sprintf("Compton cooling per electron at z = 19, T = 3000 K: %.4e vs Cen 1992 5.65e-36 (1+z)⁴ (T − 2.7(1+z)) = %.4e (ratio %.4f); sign: heats for T < T_CMB", Λc, cen, Λc / cen))
end

# ====================================================================== b, c. collisional ionisation equilibrium and the cooling curve
T5(T) = T / 1e5; g1(T) = 1 / (1 + sqrt(T5(T)))
"Cen (1992) / Katz, Weinberg & Hernquist (1996) CIE (independent comparator): ionisation rates ΓH, ΓHe, ΓHe+ and recombinations (A for H, rad+DR for He⁺)."
function cen_cie(T)
    ΓH = 5.85e-11 * sqrt(T) * g1(T) * exp(-157809.1 / T); ΓHe = 2.38e-11 * sqrt(T) * g1(T) * exp(-285335.4 / T); ΓHep = 5.68e-12 * sqrt(T) * g1(T) * exp(-631515 / T)
    aH = 8.4e-11 / sqrt(T) * (T / 1e3)^-0.2 / (1 + (T / 1e6)^0.7)
    aHep = 1.5e-10 * T^-0.6353 + 1.9e-3 * T^-1.5 * exp(-470000 / T) * (1 + 0.3 * exp(-94000 / T))
    aHepp = 3.36e-10 / sqrt(T) * (T / 1e3)^-0.2 / (1 + (T / 1e6)^0.7)
    rH = ΓH / aH; xHp = rH / (1 + rH); r1 = ΓHe / aHep; r2 = ΓHep / aHepp
    xHe = yHe / (1 + r1 + r1 * r2); xHep = xHe * r1; xHepp = xHep * r2
    (xHp = xHp, xHe = xHe, xHep = xHep, xHepp = xHepp, ne = xHp + xHep + 2xHepp)
end
"Cen (1992) / KWH96 cooling per n_H² [erg cm³ s⁻¹] at composition c: excitation, ionisation, recombination (case B + dielectronic), bremsstrahlung."
function cen_cool(T, c)
    ne = c.ne; nHI = 1 - c.xHp; g = 1.1 + 0.34 * exp(-(5.5 - log10(T))^2 / 3)
    ce = 7.5e-19 * g1(T) * exp(-118348 / T) * ne * nHI + 5.54e-17 * T^-0.397 * g1(T) * exp(-473638 / T) * ne * c.xHep + 9.1e-27 * T^-0.1687 * g1(T) * exp(-13179 / T) * ne^2 * c.xHep
    ci = 1.27e-21 * sqrt(T) * g1(T) * exp(-157809.1 / T) * ne * nHI + 9.38e-22 * sqrt(T) * g1(T) * exp(-285335.4 / T) * ne * c.xHe + 4.95e-22 * sqrt(T) * g1(T) * exp(-631515 / T) * ne * c.xHep + 5.01e-27 * T^-0.1687 * g1(T) * exp(-55338 / T) * ne^2 * c.xHep
    re = 8.70e-27 * sqrt(T) * (T / 1e3)^-0.2 / (1 + (T / 1e6)^0.7) * ne * c.xHp + 1.55e-26 * T^0.3647 * ne * c.xHep + 1.24e-13 * T^-1.5 * exp(-470000 / T) * (1 + 0.3 * exp(-94000 / T)) * ne * c.xHep + 3.48e-26 * sqrt(T) * (T / 1e3)^-0.2 / (1 + (T / 1e6)^0.7) * ne * c.xHepp
    br = 1.42e-27 * g * sqrt(T) * ne * (c.xHp + c.xHep + 4c.xHepp)
    ce + ci + re + br
end
nHc = 1e-2; ρc = nHc * CN.μH()
Tgrid = 10.0 .^ range(4.0, 6.0, length = 41)
cie = Dict{Float64,Any}()
for T in Tgrid
    c0 = cen_cie(T); x0 = [s === :Hp ? c0.xHp : s === :Hep ? c0.xHep : s === :Hepp ? c0.xHepp : 0.0 for s in net.free]
    cie[T] = CN.full(CN.steady(x0, ρc, T, net; env = env0, photos = ph), net)
end
let dH = 0.0, dHe = 0.0, dHe2 = 0.0, fHep = Float64[], fHepp = Float64[]
    for T in Tgrid
        d = cie[T]; c = cen_cie(T)
        T >= 1.5e4 && (dH = max(dH, abs(d[:Hp] - c.xHp)))
        dHe = max(dHe, abs(d[:Hep] / yHe - c.xHep / yHe)); dHe2 = max(dHe2, abs(d[:Hepp] / yHe - c.xHepp / yHe))
        push!(fHep, d[:Hep] / yHe); push!(fHepp, d[:Hepp] / yHe)
    end
    Tpk = Tgrid[argmax(fHep)]; Tx = Tgrid[findfirst(>(0.5), fHepp)]
    c0 = [Tgrid[argmax([cen_cie(T).xHep for T in Tgrid])], Tgrid[findfirst(T -> cen_cie(T).xHepp / yHe > 0.5, Tgrid)]]
    check("b CIE ionisation fractions", dH < 0.1 && dHe < 0.2 && dHe2 < 0.2 && 3e4 <= Tpk <= 6e4 && 6e4 <= Tx <= 1.3e5,
          @sprintf("vs Cen 1992 (J = 0, n_H = 1e-2): max |Δx_H⁺| = %.3f (T ≥ 1.5e4), max |Δ(He⁺/He)| = %.3f, |Δ(He²⁺/He)| = %.3f; He⁺ peak at %.2e K (Cen %.2e), He²⁺/He = 0.5 at %.2e K (Cen %.2e)", dH, dHe, dHe2, Tpk, c0[1], Tx, c0[2]))
end
println("\nCIE cooling Λ/n_H² [erg cm³ s⁻¹] (ours at own CIE / ours at Cen's composition / Cen's total at Cen's composition):")
let rown = Float64[], rowc = Float64[], Λo = Float64[]
    for T in Tgrid
        d = cie[T]; pt = A.parts(ρc, T, d; photos = ph, env = envz); Λ = (pt.rad + pt.exc + pt.dr + pt.ff) / nHc^2
        c = cen_cie(T); dC = Dict(:H => 1 - c.xHp, :Hp => c.xHp, :He => c.xHe, :Hep => c.xHep, :Hepp => c.xHepp, :e => c.ne, :Hm => 0.0)
        pc = A.parts(ρc, T, dC; photos = ph, env = envz); ΛC = (pc.rad + pc.exc + pc.dr + pc.ff) / nHc^2; cc = cen_cool(T, c)
        push!(rown, Λ / cc); push!(rowc, ΛC / cc); push!(Λo, Λ)
        abs(log10(T) * 4 - round(log10(T) * 4)) < 1e-6 && @printf("  T = %8.0f  %.3e  %.3e  %.3e   (%.2f, %.2f)\n", T, Λ, ΛC, cc, Λ / cc, ΛC / cc)
    end
    sel(lo, hi) = findall(T -> lo <= T <= hi, Tgrid)
    pk(lo, hi) = (i = sel(lo, hi); Tgrid[i[argmax(Λo[i])]])
    Th = pk(1.2e4, 3.5e4); Thep = pk(5e4, 2.5e5)
    check("c1 cooling curve", all(0.75 .< rowc .< 1.25) && all(0.65 .< rown .< 1.35) && 1.5e4 <= Th <= 2.5e4 && 6e4 <= Thep <= 1.5e5,
          @sprintf("1e4–1e6 K: Λ(ours, Cen's composition)/Λ(Cen 1992) = %.2f–%.2f (±25 %%), Λ(own CIE)/Λ(Cen) = %.2f–%.2f (±35 %%); H peak at %.2e K, He⁺ peak at %.2e K. Reference: Cen 1992/KWH96 fits coded independently (SD93 tables not obtainable offline)", minimum(rowc), maximum(rowc), minimum(rown), maximum(rown), Th, Thep))
end

# ====================================================================== d. energy conservation of the network step
coolz = A.cooling(ph, envz)
let nH = 1e8, ρ = nH * CN.μH(), T = 3e4
    x = CN.steady([0.05, 0, 0, 0, 0, 0.0], ρ, 1.0e4, net; env = env0, photos = ph)          # partly ionised gas at 1e4 K, then heated: out of equilibrium
    u0 = CN.eos(ρ, 1.0e4, x, net)[1] + 3.0e12; T = 1.0e4; T0 = 1e4
    # closed box, no cooling: u conserved (T drops as the gas ionises)
    xs, Ts = x, 2.5e4; us = CN.eos(ρ, Ts, xs, net)[1]; err = 0.0
    for dt in (1e2, 1e4, 1e6, 1e8)
        xs, Ts = CN.step(xs, ρ, us, Ts, dt, net; env = env0, photos = ph)
        err = max(err, abs(CN.eos(ρ, Ts, xs, net)[1] - us) / abs(us))
    end
    check("d1 closed box, no cooling", err < 1e-10, @sprintf("4 steps (1e2–1e8 s) from T = 2.5e4 K (x_H⁺ = 0.05 → %.3f, T → %.0f K): max |Δu|/|u| = %.1e", CN.full(xs, net)[:Hp], Ts, err))
    # with cooling: Δu = −dt Λ(T_new, x_new)/ρ per step (backward Euler), and Σ Δu = −∫Λ dt/ρ
    xs, Ts = CN.steady([1.0, 0, 0, 0, 0, yHe * 0.5], ρ, 1e5, net; env = env0, photos = ph), 1e5
    err = 0.0; Esum = 0.0; us0 = CN.eos(ρ, Ts, xs, net)[1]; nst = 0
    while Ts > 2e4 && nst < 400
        u = CN.eos(ρ, Ts, xs, net)[1]; Λ0 = coolz(ρ, Ts, CN.full(xs, net)); dt = 0.02 * u * ρ / Λ0
        xn, Tn = CN.step(xs, ρ, u, Ts, dt, net; env = env0, photos = ph, cooling = coolz)
        Λ1 = coolz(ρ, Tn, CN.full(xn, net)); un = CN.eos(ρ, Tn, xn, net)[1]
        err = max(err, abs(un - u + dt * Λ1 / ρ) / abs(u)); Esum += 0.5 * (Λ0 + Λ1) * dt / ρ
        xs, Ts = xn, Tn; nst += 1
    end
    usf = CN.eos(ρ, Ts, xs, net)[1]
    check("d2 energy with cooling", err < 1e-10 && abs((us0 - usf) / Esum - 1) < 0.02,
          @sprintf("%d steps 1e5 → %.0f K: per-step |Δu + dt Λ/ρ|/|u| ≤ %.1e; Σ trapezoid(Λ dt)/ρ = %.4e vs u drop %.4e erg/g (ratio %.4f)", nst, Ts, err, Esum, us0 - usf, (us0 - usf) / Esum))
end

# ====================================================================== e. one-zone isochoric cooling from 1e5 K
function coolrun(nH, T0; Tstops = (3e4, 2e4, 1e4))
    ρ = nH * CN.μH(); x = cie_start(nH, T0)
    T = T0; t = 0.0; out = Dict{Float64,Any}(); nst = 0
    while T > minimum(Tstops) && nst < 5000
        d = CN.full(x, net); u = CN.eos(ρ, T, x, net)[1]; Λ = coolz(ρ, T, d); dt = 0.03 * u * ρ / Λ
        Tp = T; x, T = CN.step(x, ρ, u, T, dt, net; env = env0, photos = ph, cooling = coolz); t += dt; nst += 1
        for Ts in Tstops
            if Tp > Ts >= T; w = (Tp - Ts) / (Tp - T); out[Ts] = (t = t - dt * (1 - w), d = CN.full(x, net)); end
        end
    end
    out
end
cie_start(nH, T) = CN.steady([1.0, 0, 0, 0, 0, yHe * 0.9], nH * CN.μH(), T, net; env = env0, photos = ph)
"Isochoric cooling time from T0 to T1 along the CIE cooling curve (Cen-composition independent estimate): ∫ (3/2) k n_tot dT / (n_H² Λ)."
function t_cie(nH, T0, T1)
    n = 400; Ts = 10.0 .^ range(log10(T0), log10(T1), length = n + 1); t = 0.0
    for i in 1:n
        Tm = sqrt(Ts[i] * Ts[i+1]); c = cen_cie(Tm); ntot = nH * (1 + yHe + c.xHep + c.xHepp + c.ne)
        t += 1.5 * CN.kB * ntot * (Ts[i] - Ts[i+1]) / (nH^2 * cen_cool(Tm, c))
    end
    t
end
println("\nIsochoric cooling from 1e5 K (thin, J = 0):")
let ok = true
    res = Dict(nH => coolrun(nH, 1e5) for nH in (1e3, 1e10))
    for nH in (1e3, 1e10)
        r = res[nH]
        for Ts in (3e4, 2e4, 1e4)
            haskey(r, Ts) || (ok = false; continue)
            d = r[Ts].d; dc = cie[Tgrid[argmin(abs.(log.(Tgrid) .- log(Ts)))]]
            @printf("  n_H = %.0e: T = %.0e K reached at t = %.3e yr (CIE-curve estimate %.3e yr); x_H⁺ = %.3f (CIE %.3f), He⁺/He = %.3f (CIE %.3f), He²⁺/He = %.3f (CIE %.3f)\n", nH, Ts, r[Ts].t / 3.156e7, t_cie(nH, 1e5, Ts) / 3.156e7, d[:Hp], dc[:Hp], d[:Hep] / yHe, dc[:Hep] / yHe, d[:Hepp] / yHe, dc[:Hepp] / yHe)
        end
    end
    r = res[1e3]; lag1 = r[1e4].d[:Hp] > 3 * cie[Tgrid[argmin(abs.(log.(Tgrid) .- log(1e4)))]][:Hp]
    r2 = res[1e3]; lag2 = r2[3e4].d[:Hepp] / yHe > 1.5 * cie[Tgrid[argmin(abs.(log.(Tgrid) .- log(3e4)))]][:Hepp] / yHe
    tr2 = res[1e3][2e4].t / t_cie(1e3, 1e5, 2e4); tr1 = res[1e3][1e4].t / t_cie(1e3, 1e5, 1e4); sc = res[1e3][1e4].t / res[1e10][1e4].t
    check("e isochoric cooling", ok && lag1 && lag2 && 0.5 < tr2 < 3 && tr1 < 0.5 && 1e6 < sc < 1e8,
          @sprintf("n = 1e3: t(1e5 → 2e4 K) = %.1f yr = %.2f × the CIE-curve time (the He²⁺ left over from 1e5 K cools less efficiently than He⁺ in CIE), t(→ 1e4 K) = %.1f yr = %.2f × (the gas stays ionised, x_H⁺ = %.3f at 1e4 K vs CIE %.4f, so Lyα and recombination keep cooling it: Gnat & Sternberg 2007 find the same, non-equilibrium cooling being faster below ~2e4 K); He²⁺/He at 3e4 K %.3f vs CIE %.4f; t(n = 1e3)/t(n = 1e10) = %.1e (thin: ∝ 1/n)", res[1e3][2e4].t / 3.156e7, tr2, res[1e3][1e4].t / 3.156e7, tr1, r[1e4].d[:Hp], cie[Tgrid[argmin(abs.(log.(Tgrid) .- log(1e4)))]][:Hp], r[3e4].d[:Hepp] / yHe, cie[Tgrid[argmin(abs.(log.(Tgrid) .- log(3e4)))]][:Hepp] / yHe, sc))
end

# ====================================================================== f. escape probability
let dd = Dict(:H => 0.99, :Hp => 0.01, :He => yHe, :Hep => 0.0, :Hepp => 0.0, :e => 0.01, :Hm => 0.0), T = 1.0e4, nH = 1e6, ρ = nH * CN.μH()
    A.set_zone!(0.0, Inf, nH, dd); thin = A.parts(ρ, T, dd; photos = ph, env = envz).exc
    rs = Float64[]; τs = Float64[]
    for r in 10.0 .^ (0.0:1.0:18.0)
        A.set_zone!(r * 1e0, Inf, nH, dd)               # static column ℓ = r
        τ = A.πe2mc * 0.4162 * (nH * dd[:H]) * r / (sqrt(π) * (CN.eV * 10.2 / CN.hP) * sqrt(2CN.kB * T / CN.mH) / CN.c); push!(τs, τ)
        push!(rs, A.parts(ρ, T, dd; photos = ph, env = envz).exc / thin)
    end
    A.set_zone!(0.0, Inf, nH, dd)
    # Lyα dominates the H excitation: thick slope between the last decades
    sl = log(rs[end] / rs[end-1]) / log(10.0); thinerr = abs(rs[1] - 1)
    check("f escape probability", thinerr < 1e-6 && abs(sl + 1) < 0.02 && 0.1 < rs[findfirst(>(1.0), τs)] < 0.4,
          @sprintf("Λ_exc(τ)/Λ_exc(thin): %.6f at τ_Lyα = %.1e, %.3f at %.1e, %.2e at %.1e, %.2e at %.1e; log-slope %.3f (thick, → −1); the weaker lines (Lyβ…) turn over later, so the slope is that of the sum", rs[1], τs[1], rs[findfirst(>(1.0), τs)], τs[findfirst(>(1.0), τs)], rs[end-1], τs[end-1], rs[end], τs[end], sl))
    # on-the-spot case B: ground-level recombination continuum absorbed in a thick neutral medium
    A.set_zone!(1e14, Inf, 1e10, Dict(:H => 1.0, :He => yHe, :Hep => 0.0)); pg = [p for p in ph if p.R === :H && endswith(p.name, "ground")][1]
    eA = pg.esc(T, envz); A.set_zone!(0.0, Inf, 1e10, dd); eB = pg.esc(T, envz)
    check("f2 on-the-spot (case B)", eA < 1e-6 && eB == 1.0, @sprintf("ground-level recombination continuum escape fraction %.1e in a neutral column of 1e24 cm⁻² (case B), %.1f when thin (case A)", eA, eB))
end

# ====================================================================== g. equilibrium table (PopIIIEqEOS) and He
let
    function hefrac(p, T)
        nH = p / (CN.kB * T) / 1.1
        for _ in 1:200
            K1 = exp(CN.lnKeq([:H], [:Hp, :e], T)); xHp = 0.5; ne = nH * xHp
            for _ in 1:200; xHp = 1 / (1 + ne / K1); ne = xHp * nH; end
            r1 = exp(CN.lnKeq([:He], [:Hep, :e], T)) / max(ne, 1e-30); r2 = exp(CN.lnKeq([:Hep], [:Hepp, :e], T)) / max(ne, 1e-30)
            xHe = yHe / (1 + r1 + r1 * r2); xHep = xHe * r1; xHepp = xHep * r2
            nn = p / (CN.kB * T) / (1 + 2xHp + yHe + xHep + 2xHepp)
            abs(nn / nH - 1) < 1e-10 && return (xHp, (xHep + xHepp) / yHe, (yHe * (xHep / yHe) * CN.χHe1 + xHepp * (CN.χHe1 + CN.χHe2)) / (CN.χH * xHp + 1e-300))
            nH = nn
        end
        (NaN, NaN, NaN)
    end
    worst = 0.0; Eworst = 0.0
    for p in (1e4, 1e5, 1e6, 1e7, 1e8)
        _, f, E = hefrac(p, 1.2e4); worst = max(worst, f); Eworst = max(Eworst, E)
    end
    f20 = hefrac(1e4, 2e4)
    check("g eq. table He", worst < 1e-3 && Eworst < 3e-3, @sprintf("Saha He ionisation at the table's pressures (1e4–1e8 dyn cm⁻²) and T = 1.2e4 K (NetGas.Teqmax of the atomic net): ≤ %.1e of He, energy ≤ %.1e of the H ionisation energy: the table's neutral He is exact below; at 2e4 K and p = 1e4: He ionised %.2f (He ionisation energy %.2f of H's) — those zones stay in the network", worst, Eworst, f20[2], f20[3]))
end

# ====================================================================== h. stiff cooling step (hot ionised infall zone, t_cool ≈ 4 s, dt = 9e6 s)
let ρ = 5.791258057811794e-14, T0 = 50393.782336134165, u = 9.739868491145969e12, dt = 8.8847562265625e6
    netL, phL = A.atomic(li = false); coolL = A.cooling(phL, envz); CN.XFLOOR[] = 1e-12
    x0 = [s === :Hp ? 0.9997238235077323 : s === :Hm ? 1.6462444167136207e-13 : s === :H2 ? 7.930953103803236e-17 : s === :H2p ? 4.9467028716105e-10 : 0.0 for s in netL.free]
    plain = try CN.step(x0, ρ, u, T0, dt, netL; env = env0, photos = phL, cooling = coolL); true catch; false end
    CN.NESTED[] = true
    xn, Tn = CN.step(x0, ρ, u, T0, dt, netL; env = env0, photos = phL, cooling = coolL)
    CN.NESTED[] = false; CN.XFLOOR[] = 1e-20
    un = CN.eos(ρ, Tn, xn, netL)[1]; res = abs(un - u + dt * coolL(ρ, Tn, CN.full(xn, netL)) / ρ) / abs(u)
    check("h stiff cooling step", res < 1e-9 && 4e3 < Tn < 1.5e4,
          @sprintf("ionised gas at 5.0e4 K, n_H = 2.6e10, Λ t_cool ≈ 4 s, step of 9e6 s: coupled Newton %s; with ChemNet.NESTED (bisection in T over the chemistry): T → %.0f K, energy residual %.1e", plain ? "converges" : "fails", Tn, res))
end

nf = count(!last, RESULTS)
println(nf == 0 ? "S6C ATOMIC: ALL PASS ($(length(RESULTS)) checks)" : "S6C ATOMIC: $nf FAIL of $(length(RESULTS))")
