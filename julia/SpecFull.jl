"""
    SpecFull

Full emergent spectrum of a recorded snapshot: continuum from 1000 μm to 0.01 μm (past the He II edge at 228 Å) plus the atomic and
molecular lines in emission AND absorption, in the observer frame with the snapshot's velocity field. Post-processing only:
nothing here is loaded by a run. Driver: `snapfull.jl`. Plan: scratchpad spec_lines_plan.md (13 Nov 2026).

Structure (what is computed, and the approximations):

1. Points (`points`). Zone centres of the snapshot outside the run's grey optical depth TAUCUT (as snapspec), with the run's
   non-equilibrium composition (H, H⁺, H⁻, H₂, e⁻ and, when recorded, He⁺/He²⁺). Zones that are optically thick in the run's grey
   opacity (Δτ > DTMAX) and lie in its outer TREF are split into sub-points (resolved T(τ), the prerequisite of the plan, option a):
   on the run's own grey τ scale (τ linear in r inside the zone, the zone's χ and ρ), T⁴ = fT ¾ F/σ_SB (τ + ⅔), F the run's radiative
   flux at the zone's outer interface, capped at the zone's T; the composition there is LTE at the sub-point's T (Z0Physics.
   opacity_state: H, He ionisation with occupation probabilities, H⁻ by Saha; H₂ only below 4000 K). fT is iterated (≤ 5 continuum
   solutions) until the frequency-dependent continuum carries the run's luminosity through the zone's outer interface: a one-parameter
   Unsöld–Lucy correction. MODEL ASSUMPTION, used only where the run's mesh puts the whole photosphere into one zone. c128 at 240 M☉
   (13 Nov): one 21,700 K zone of 5,600–7,500 R☉, its grey χ 0.47 cm²/g (mostly Thomson at the zone's mean T); without the sub-points
   the continuum radiated the zone-mean T (10³ × the run's L); with the run's τ scale and fT = 1 it was 7 × too bright (the real
   opacity of the cool sub-layers is far below the zone's χ); with fT = 0.0077 it carries 1.04 × the run's L. A fixed-ρ sub-zone with
   its own Rosseland opacity does not work: at the zone's mean ρ the cool gas is transparent far below the zone (tried, 13 Nov).
   Where T or ρ jumps by more than ×3 between
   unrefined zones, two points carry the two sides (Spectrum.atmosphere's rule).
2. Level populations (`populations!`). H I (n ≤ NH), He I (ModelAtom.heI, nmax NHE1) and He II (n ≤ NHE2) in statistical equilibrium
   at the point's T and n_e with the ionisation balance of the run (ground state and next ion fixed to the snapshot's densities:
   ModelAtom.solve_upper), in the continuum mean intensity J_ν of this module's own continuum transfer, with Sobolev escape
   probabilities in the snapshot's velocity field (thermal-width floor, as snapnlte.jl), Λ-iterated on β. Two passes: first with
   the run's radiation (B(T) where the grey τ > 1, diluted B(T_eff) outside, snapnlte's prescription), then with the J_ν of the first
   continuum solution. Points deeper than grey τ = TLTE, and the sub-points of 1., are in LTE relative to the next ion
   (n_i = Φ_i n_e n₊). H₂: the run's ALI level populations (snapspec h2lines → H2Lines.LAST); none in the merged core. HD: LTE at
   the gas T (`hdload`: Amaral et al. 2019, ExoMol ADJSAAM; no HD collision data in the code), with the run's HD abundance.
3. Continuum opacity and emission (`coeffs3` for `FOpac`): ContinuumOpacity's CIA (H₂–H₂, H₂–He, H₂–H), H⁻ free–free, Rayleigh and
   Thomson; H⁻ bound–free with attachment emission; bound–free of every level of H I, He I, He II from the populations of 2. with
   their own recombination emission (Milne: j_ν = σ n_i* 2hν³/c² e^{−hν/kT}, n_i* = Φ_i n_e n₊; stimulated recombination subtracted
   from κ); free–free on H⁺, He⁺, He²⁺ with Gaunt factors (Z0Physics.gaunt_ff), thermal. Coherent isotropic scattering. Transfer:
   Spectrum.monochromatic (tangent rays, Feautrier ALI), static (continuum features are wider than v/c).
4. Lines (`windows`). Each line of LINES gets a velocity window ±VWIN at DV resolution. Opacity κ = (πe²/mc) f (n_l − g_l n_u/g_u) φ,
   emissivity η = (hν/4π) n_u A φ, with φ the Voigt profile of Z0Physics.line_profile at each point: Doppler (thermal + XI
   microturbulence), radiative damping, electron impact; for H I and He II the exact linear Stark pattern of the ions with the
   Holtsmark distribution. Overlapping lines in a window add. Continuum inside the window: κ, σ, emission at the line frequency, the
   scattering source σ J_c from 3. Formal solution in the OBSERVER frame along the tangent rays (far side → near side, core rays from
   the thermal inner boundary): the comoving frequency ν(1 − v_z/c), v_z = v_r z/r, and every segment subdivided so that v_z changes by
   less than SUBV of the local thermal width, the profiles of both end points evaluated at each sub-point's shift. Infall therefore
   gives the redshifted absorption of an inverse P Cygni profile. Not included: the line J in the populations beyond Sobolev,
   continuum scattering out of the line (electron scattering broadening), line photons of one window feeding another.
5. Output per snapshot: the continuum (L_ν at the edge and through the photosphere, r(τ_ν = ⅔)), per line the observed L_ν profile on
   the window, the continuum level, the line luminosity ∫(L_ν − L_c) dν and equivalent width (negative = emission), the radial
   optical depth at line centre (comoving, static) on a log-r grid.

Validation (snapfull.jl --test): LTE slab (J = B, β arbitrary, isothermal): no line features (|L/L_c − 1| below the solver error).
"""
module SpecFull

import ..Spectrum, ..ContinuumOpacity, ..ModelAtom, ..Z0Physics, ..H2Lines
using Printf, LinearAlgebra
const SP = Spectrum; const CO = ContinuumOpacity; const MA = ModelAtom; const Z0 = Z0Physics
const c = CO.c; const h = CO.h; const kB = CO.kB
const yHe = 1 / 12; const mH = 1.6735575e-24; const μH = 1.6735e-24 * (1 + 4yHe)
const Rs = 6.957e10; const Ls = 3.83e33; const σSB = 5.670374e-5
const TAUCUT = 300.0; const DTMAX = 0.5; const TREF = 30.0; const TDEEP = 200.0; const TLTE = 50.0
const NH = 20; const NHE1 = 10; const NHE2 = 20; const XI = 2e5          # microturbulence 2 km/s
const VWIN = 6e7; const DV = 6e5; const SUBV = 0.3; const NRT = 80          # windows ±600 km/s at 6 km/s; log-r grid of τ_line

const ATOMS = Ref{Any}(nothing)
"The three model atoms (built once): H I and He II with CHIANTI collisions, He I with CHIANTI terms + hydrogenic shells."
function atoms()
    ATOMS[] === nothing || return ATOMS[]
    ad = joinpath(@__DIR__, "atomdata", "chianti")
    H = MA.hydrogenic(1; nmax = NH, coll = :chianti, chianti = joinpath(ad, "h_1"))
    He2 = MA.hydrogenic(2; nmax = NHE2, coll = :chianti, chianti = joinpath(ad, "he_2"))
    He1 = MA.heI(; nmax = NHE1)
    ATOMS[] = (; H, He1, He2, LH = MA.lines(H), L1 = MA.lines(He1), L2 = MA.lines(He2))
end

"Continuum frequency grid [Hz], ascending: log-uniform 1000 μm … 0.01 μm plus pairs ±0.1 % at the H n = 1–3, He I and He II edges."
function nugrid(; n = 250)
    a = atoms(); ν = [c / (1e-4 * λ) for λ in exp10.(range(3, -2, length = n))]
    for e in (a.H.levels[1].χ, a.H.levels[2].χ, a.H.levels[3].χ, a.He1.levels[1].χ, a.He2.levels[1].χ)
        ν0 = e / h; push!(ν, ν0 * 0.999, ν0 * 1.001)
    end
    sort!(ν)
end

"""
The line list: (name, vacuum λ [Å] label, atom, lower, upper) with levels by index (H I, He II: n) or by term label (He I);
H₂ lines by their vacuum wavelength in μm (matched to the H2Lines list within 0.1 %).
"""
const LINEDEF = [
    ("Lyα", :H, 1, 2), ("Lyβ", :H, 1, 3), ("Hα", :H, 2, 3), ("Hβ", :H, 2, 4), ("Hγ", :H, 2, 5), ("Hδ", :H, 2, 6),
    ("Paα", :H, 3, 4), ("Paβ", :H, 3, 5), ("Paγ", :H, 3, 6), ("Brα", :H, 4, 5), ("Brγ", :H, 4, 7),
    ("He II 304", :He2, 1, 2), ("He II 256", :He2, 1, 3), ("He II 1640", :He2, 2, 3), ("He II 1215", :He2, 2, 4),
    ("He II 4686", :He2, 3, 4), ("He II 3203", :He2, 3, 5), ("He II 10124", :He2, 4, 5),
    ("He I 584", :He1, "1s2 1S", "1s.2p 1P"), ("He I 537", :He1, "1s2 1S", "1s.3p 1P"), ("He I 10830", :He1, "1s.2s 3S", "1s.2p 3P"),
    ("He I 3889", :He1, "1s.2s 3S", "1s.3p 3P"), ("He I 5876", :He1, "1s.2p 3P", "1s.3d 3D"), ("He I 4471", :He1, "1s.2p 3P", "1s.4d 3D"),
    ("He I 6678", :He1, "1s.2p 1P", "1s.3d 1D"), ("He I 7065", :He1, "1s.2p 3P", "1s.3s 3S"), ("He I 20581", :He1, "1s.2s 1S", "1s.2p 1P"),
    ("H2 1-0 S(1)", :H2, 2.1218, 0), ("H2 1-0 S(0)", :H2, 2.2233, 0), ("H2 1-0 Q(1)", :H2, 2.4066, 0), ("H2 2-1 S(1)", :H2, 2.2477, 0),
    ("H2 0-0 S(0)", :H2, 28.221, 0), ("H2 0-0 S(1)", :H2, 17.035, 0), ("H2 0-0 S(2)", :H2, 12.279, 0), ("H2 0-0 S(3)", :H2, 9.6649, 0),
    ("H2 0-0 S(5)", :H2, 6.9091, 0), ("H2 0-0 S(7)", :H2, 5.5115, 0),
    ("HD 0-0 R(0)", :HD, (0, 1), (0, 0)), ("HD 0-0 R(1)", :HD, (0, 2), (0, 1)), ("HD 0-0 R(2)", :HD, (0, 3), (0, 2)),
    ("HD 0-0 R(3)", :HD, (0, 4), (0, 3)), ("HD 0-0 R(5)", :HD, (0, 6), (0, 5)), ("HD 1-0 R(0)", :HD, (1, 1), (0, 0)),
    ("HD 1-0 P(1)", :HD, (1, 0), (0, 1)),
]

"""
HD line list of Amaral et al. (2019, ApJ 878, 95; ExoMol dataset ADJSAAM, files 1H-2H__ADJSAAM.states/.trans from exomol.com,
retrieved 13 Nov 2026, in ENV HDLINES or the scratch data folder): states (id, E [cm⁻¹], g incl. nuclear spin, J, lifetime, v,
bound/quasi-bound) and transitions (upper, lower, A [s⁻¹], ν̃ [cm⁻¹]). The code has no HD collision data: the level populations
are LTE at the gas T (Boltzmann over all bound states), which overestimates the lines where n is below their critical density.
"""
const HDL = Ref{Any}(nothing)
function hdload(dir = get(ENV, "HDLINES", "/zpool/nvme/data/tabel_scratch/titan_popiii/lit/HDlines"))
    HDL[] === nothing || return HDL[]
    fs = joinpath(dir, "1H-2H__ADJSAAM.states"); isfile(fs) || return nothing
    st = [split(l) for l in readlines(fs) if !isempty(strip(l))]; n = length(st)
    E = zeros(n); g = zeros(n); J = zeros(Int, n); v = zeros(Int, n); bound = falses(n)
    for q in st
        k = parse(Int, q[1]); E[k] = parse(Float64, q[2]) * h * c; g[k] = parse(Float64, q[3]); J[k] = parse(Int, q[4])
        v[k] = parse(Int, q[6]); bound[k] = q[7] == "b"
    end
    iu = Int[]; il = Int[]; A = Float64[]; ν = Float64[]
    for l in readlines(joinpath(dir, "1H-2H__ADJSAAM.trans"))
        q = split(l); length(q) < 4 && continue
        push!(iu, parse(Int, q[1])); push!(il, parse(Int, q[2])); push!(A, parse(Float64, q[3])); push!(ν, parse(Float64, q[4]) * c)
    end
    HDL[] = (; E, g, J, v, bound, iu, il, A, ν)
end
"LTE level fractions of HD at T (bound states)."
function hdfrac(T)
    L = hdload(); w = [L.bound[k] ? L.g[k] * exp(-L.E[k] / (kB * T)) : 0.0 for k in eachindex(L.E)]; w ./ sum(w)
end

"""
Resolved line list: per line (name, atom, l, u, ν0, f, A, g_l, g_u, m [g], Z, n_l*, n_u* (effective principal quantum numbers for
the impact width)). Lines whose levels are missing from the atoms (He I terms above NHE1, H₂ lines outside the list) are dropped.
"""
function linelist(HL)
    a = atoms(); out = []
    lab1 = [l.label for l in a.He1.levels]
    for d in LINEDEF
        nm, at = d[1], d[2]
        if at === :H || at === :He2
            A_ = at === :H ? a.H : a.He2; l, u = d[3], d[4]; Z = at === :H ? 1 : 2
            u <= MA.nlev(A_) || continue
            ν0 = MA.νline(A_, l, u)
            push!(out, (; name = nm, atom = at, l, u, ν0, f = A_.f[l, u], A = A_.A[u, l], gl = A_.levels[l].g, gu = A_.levels[u].g,
                        m = at === :H ? mH : 4.0026 * 1.66054e-24, Z, nsl = Float64(l), nsu = Float64(u)))
        elseif at === :He1
            l = findfirst(==(d[3]), lab1); u = findfirst(==(d[4]), lab1)
            (l === nothing || u === nothing) && continue
            A_ = a.He1; Aul = A_.A[u, l]; Aul > 0 || continue
            ν0 = MA.νline(A_, l, u); Ry = 13.605693 * 1.602176634e-12
            push!(out, (; name = nm, atom = at, l, u, ν0, f = A_.f[l, u], A = Aul, gl = A_.levels[l].g, gu = A_.levels[u].g,
                        m = 4.0026 * 1.66054e-24, Z = 1, nsl = sqrt(Ry / A_.levels[l].χ), nsu = sqrt(Ry / A_.levels[u].χ)))
        elseif at === :HD
            HD = hdload(); HD === nothing && continue
            ku = findfirst(k -> HD.v[k] == d[3][1] && HD.J[k] == d[3][2] && HD.bound[k], eachindex(HD.E))
            kl = findfirst(k -> HD.v[k] == d[4][1] && HD.J[k] == d[4][2] && HD.bound[k], eachindex(HD.E))
            (ku === nothing || kl === nothing) && continue
            i = findfirst(t -> HD.iu[t] == ku && HD.il[t] == kl, eachindex(HD.iu)); i === nothing && continue
            ν0 = HD.ν[i]; A = HD.A[i]; gu = HD.g[ku]; gl = HD.g[kl]; f = A * gu / (MA.ACOEF * ν0^2 * gl)
            push!(out, (; name = nm, atom = :HD, l = kl, u = ku, ν0, f, A, gl, gu, m = 3mH, Z = 0, nsl = 0.0, nsu = 0.0))
        else
            HL === nothing && continue
            λ = d[3]; i = argmin(abs.(1e4 .* c ./ HL.ν .- λ)); abs(1e4 * c / HL.ν[i] / λ - 1) < 1e-3 || continue
            ν0 = HL.ν[i]; gu = HL.g[HL.iu[i]]; gl = HL.g[HL.il[i]]; A = HL.A[i]
            f = A * gu / (MA.ACOEF * ν0^2 * gl)
            push!(out, (; name = nm, atom = :H2, l = HL.il[i], u = HL.iu[i], ν0, f, A, gl, gu, m = 2mH, Z = 0, nsl = 0.0, nsu = 0.0))
        end
    end
    out
end

# ---------------------------------------------------------------- points
"""
Points of the transfer, inner to outer: r, ρ, T, v_r, dv/dr, zone, the zone width, grey τ, lte flag and the number densities
(nH0, nHp, nHm, nH2, ne, nHe0, nHe1, nHe2), and `jump` per segment. `d` = snapspec densities, `τ` the run's grey τ at interfaces.
"""
function points(x, d, τ; fT = 1.0)
    n = length(x.ρ); kc = findlast(>(TAUCUT), τ); kc = kc === nothing ? 1 : min(kc, n - 8)
    spcv(s) = (i = findfirst(==(s), x.spn); i === nothing ? nothing : max.(x.Y[i, :], 0.0))
    xHe1 = spcv(:Hep); xHe2 = spcv(:Hepp); nHt = x.ρ ./ μH
    xHD = (hasproperty(x, :hd) && x.hd !== nothing && length(x.hd.HD) == n) ? max.(x.hd.HD, 0.0) : zeros(n)
    P = NamedTuple[]; jump = Bool[]
    uat(r, k) = (w = (r - x.r[k]) / (x.r[k+1] - x.r[k]); (1 - w) * x.u[k] + w * x.u[k+1])
    function zonept(k, r)
        if xHe1 === nothing
            s = Z0.saha(x.ρ[k], x.T[k]); f1 = s.nHe1 / s.nHe; f2 = s.nHe2 / s.nHe
        else
            f1 = xHe1[k] / yHe; f2 = xHe2 === nothing ? 0.0 : xHe2[k] / yHe
        end
        nHe = yHe * nHt[k]; τc = 0.5 * (τ[k] + τ[k+1])
        (; r, ρ = x.ρ[k], T = x.T[k], v = uat(r, k), dvdr = (x.u[k+1] - x.u[k]) / (x.r[k+1] - x.r[k]), zone = k, Δr = x.r[k+1] - x.r[k],
         τg = τc, lte = τc > TLTE, nH0 = d.nH[k], nHp = d.nHp[k], nHm = d.nHm[k], nH2 = d.nH2[k], ne = d.ne[k],
         nHe0 = nHe * max(1 - f1 - f2, 0.0), nHe1 = nHe * f1, nHe2 = nHe * f2, nHD = nHt[k] * xHD[k])
    end
    refined = falses(n)
    for k in kc:n
        Δτ = τ[k] - τ[k+1]
        # only a zone through which the emergent light leaves: the run's L at its outer interface ≥ 0.3 × the run's L at the edge
        # (z100, 13 Nov: the grey τ = ⅔ point sat at the shock with 2 L☉ of net flux while 2e5 L☉ left the edge)
        Ltop_k = 4π * x.r[k+1]^2 * x.fr[k+1]; Ledge = 4π * x.r[end-2]^2 * x.fr[end-2]
        if Δτ > DTMAX && τ[k+1] < TREF && isfinite(Ltop_k) && Ltop_k >= 0.3 * Ledge && Ledge > 0
            refined[k] = true
            # sub-layers on the run's own grey τ scale (τ linear in r inside the zone: the run's χ and ρ of the zone), T⁴ = fT ¾ F/σ (τ + ⅔)
            # capped at the zone's T, F the run's flux at the zone's outer interface; fT normalises the frequency-dependent continuum
            # to the run's luminosity there (`spectrum`), a one-parameter Unsöld–Lucy correction. Composition LTE at the sub-point's T.
            τa = τ[k+1]; τb = min(τ[k], TDEEP); χρ = Δτ / (x.r[k+1] - x.r[k])
            F = x.fr[k+1]; (isfinite(F) && F > 0) || (F = σSB * max(x.Teff, 1.0)^4)
            sub = NamedTuple[]
            for Δ in exp.(range(log(1e-3), log(max(τb - τa, 2e-3)), length = 36))
                τj = τa + Δ; rj = x.r[k+1] - Δ / χρ; rj <= x.r[k] && break
                Tj = min((fT * 0.75 * F / σSB * (τj + 2 / 3))^0.25, x.T[k])
                s = Z0.opacity_state(x.ρ[k], Tj); kT = kB * Tj
                nHm = s.pH[1] * s.ne * exp(Z0.χHm / kT) / (4 * Z0.λe3(Tj))
                push!(sub, (; r = rj, ρ = x.ρ[k], T = Tj, v = uat(rj, k), dvdr = (x.u[k+1] - x.u[k]) / (x.r[k+1] - x.r[k]), zone = k,
                            Δr = x.r[k+1] - x.r[k], τg = τj, lte = true, nH0 = s.nH0, nHp = s.nHp, nHm, nH2 = Tj < 4000 ? d.nH2[k] : 0.0,
                            ne = s.ne, nHe0 = s.nHe0, nHe1 = s.nHe1, nHe2 = s.nHe2, nHD = Tj < 4000 ? nHt[k] * xHD[k] : 0.0))
            end
            for q in reverse(sub); push!(P, q); push!(jump, false); end
        else
            rc = 0.5 * (x.r[k] + x.r[k+1])
            # interface pair where T or ρ jumps by more than 3 between this and the previous (unrefined) zone
            if k > kc && !refined[k-1] && !isempty(P)
                q = max(x.T[k] / x.T[k-1], x.T[k-1] / x.T[k], x.ρ[k] / x.ρ[k-1], x.ρ[k-1] / x.ρ[k])
                if q > 3.0
                    ri = x.r[k]; ε = 0.01 * min(ri - P[end].r, rc - ri)
                    push!(P, merge(zonept(k - 1, ri - ε), (; r = ri - ε))); push!(jump, false)
                    push!(P, merge(zonept(k, ri + ε), (; r = ri + ε))); push!(jump, true)
                end
            end
            push!(P, zonept(k, rc)); push!(jump, false)
        end
    end
    # jump[i] refers to the segment ending at point i; Spectrum wants jump per segment i → i+1
    jseg = [jump[i+1] for i in 1:length(P)-1]
    kr = findlast(refined); (; P, jump = jseg, kc, nref = count(refined), kr = kr === nothing ? 0 : kr)
end

# ---------------------------------------------------------------- populations
"Log–log interpolation of J on the ascending grid ν (0 above the grid, ∝ ν² below)."
function Jinterp(ν, J)
    lν = log.(ν); lJ = log.(max.(J, 1e-300))
    function (x)
        x >= ν[end] && return 0.0
        x <= ν[1] && return J[1] * (x / ν[1])^2
        i = searchsortedlast(ν, x); w = (log(x) - lν[i]) / (lν[i+1] - lν[i])
        exp((1 - w) * lJ[i] + w * lJ[i+1])
    end
end
const RY = 13.605693 * 1.602176634e-12

"""
Populations of one atom with the ground state n1 and the next ion n₊ fixed (`MA.solve_upper`), in the continuum J and Sobolev
β (Λ-iterated, ≤ nit). LTE relative to the ion when `lte`; only the ground state below 2500 K; all ionised above 3e5 K.
"""
function atompops(a, L, T, ne, n1, nplus, J, v, r, dvdr, Qmin, lte; nit = 8)
    N = MA.nlev(a); n = zeros(N)
    if T < 2500 || ne <= 0
        n[1] = n1; return n
    elseif T > 3e5
        return n
    end
    if lte
        for i in 1:N; n[i] = MA.Φ(a, i, T) * ne * nplus; end      # (ground included: S = B exactly in bound–free and lines)
        return n
    end
    β = ones(length(L)); idx = Dict((x[1], x[2]) => k for (k, x) in enumerate(L))
    βf = (l, u) -> (k = get(idx, (l, u), 0); k == 0 ? 1.0 : β[k])
    for it in 1:nit
        n = try MA.solve_upper(a, T, ne, nplus; nfix = [n1], J, β = βf) catch; (n .= 0; n[1] = n1; return n) end
        all(isfinite, n) || (n = zeros(N); n[1] = n1; return n)
        n .= max.(n, 0.0)
        βn = [MA.βsob(MA.BCOEF * x[7] * (c / x[3]) * max(n[x[1]] - a.levels[x[1]].g / a.levels[x[2]].g * n[x[2]], 0.0), v, r, dvdr; Qmin) for x in L]
        δ = maximum(abs.(βn .- β); init = 0.0); β = βn
        δ < 1e-3 && break
    end
    n
end

"Level populations of H I, He I, He II at every point (threaded over points)."
function populations(pts, Jf)
    a = atoms(); np = length(pts)
    nH = Vector{Vector{Float64}}(undef, np); n1 = similar(nH); n2 = similar(nH)
    Threads.@threads for k in 1:np
        p = pts[k]; J = Jf(k)
        vth(m) = sqrt(2kB * p.T / m + XI^2)
        nH[k] = atompops(a.H, a.LH, p.T, p.ne, p.nH0, p.nHp, J, p.v, p.r, p.dvdr, vth(mH) / p.Δr, p.lte)
        n1[k] = atompops(a.He1, a.L1, p.T, p.ne, p.nHe0, p.nHe1, J, p.v, p.r, p.dvdr, vth(4mH) / p.Δr, p.lte)
        n2[k] = atompops(a.He2, a.L2, p.T, p.ne, p.nHe1, p.nHe2, J, p.v, p.r, p.dvdr, vth(4mH) / p.Δr, p.lte)
    end
    (; nH, n1, n2)
end

# ---------------------------------------------------------------- continuum coefficients
struct FOpac
    pts::Vector; pop::Any
end
const E0HM = 0.754 * CO.eV
"Bound–free of one atom at ν: (absorption incl. stimulated recombination, emission) [cm⁻¹, erg cm⁻³ s⁻¹ Hz⁻¹ sr⁻¹]."
function bf(a, n, nplus, ne, T, ν)
    ab = 0.0; em = 0.0; kT = kB * T; hν = h * ν; c2 = 2h * ν^3 / c^2; λ3 = Z0.λe3(T)
    for i in 1:MA.nlev(a)
        χ = a.levels[i].χ; hν >= χ || continue
        σ = a.σ[i](ν); σ > 0 || continue
        nst = ne * nplus * a.levels[i].g / (2a.gplus) / λ3 * exp((χ - hν) / kT)     # n_i* e^{−hν/kT}
        ab += σ * (n[i] - nst); em += σ * nst * c2
    end
    ab, em
end
function SP.coeffs3(o::FOpac, k, ν)
    p = o.pts[k]; T = p.T; kT = kB * T; hν = h * ν; stim = -expm1(-hν / kT); a = atoms()
    al, s, q = CO.alpha(ν, T; nH = p.nH0, nH2 = p.nH2, nHe = p.nHe0, nHp = p.nHp, ne = p.ne, nHm = p.nHm)
    B = CO.Bν(ν, T)
    ab = al - q.hm_bf - q.h                                                # CIA, H⁻ free–free
    ab += 3.69e8 * p.ne / (sqrt(T) * ν^3) * stim * ((p.nHp + p.nHe1) * Z0.gaunt_ff(1, T, ν) + 4p.nHe2 * Z0.gaunt_ff(2, T, ν))
    em = ab * B
    if hν > E0HM                                                           # H⁻ bound–free, attachment emission
        σ = CO.σ_Hm_bf(1e4 * c / ν)
        ab += p.nHm * σ * stim
        em += p.nH0 * p.ne / Z0.λe3(T) / 4 * (2h / c^2) * σ * ν^3 * exp(-(hν - E0HM) / kT)    # (Z0.λe3 is λ_e⁻³)
    end
    pop = o.pop
    if pop === nothing                                                     # (first pass: LTE levels of H relative to H⁺)
        nl = [n == 1 ? p.nH0 : (T > 2500 ? MA.Φ(a.H, n, T) * p.ne * p.nHp : 0.0) for n in 1:MA.nlev(a.H)]
        x1, y1 = bf(a.H, nl, p.nHp, p.ne, T, ν); ab += x1; em += y1
        hν > a.He1.levels[1].χ && (ab += p.nHe0 * a.He1.σ[1](ν) * stim; em += p.nHe0 * a.He1.σ[1](ν) * stim * B)
        hν > a.He2.levels[1].χ && (ab += p.nHe1 * a.He2.σ[1](ν) * stim; em += p.nHe1 * a.He2.σ[1](ν) * stim * B)
    else
        x1, y1 = bf(a.H, pop.nH[k], p.nHp, p.ne, T, ν); x2, y2 = bf(a.He1, pop.n1[k], p.nHe1, p.ne, T, ν); x3, y3 = bf(a.He2, pop.n2[k], p.nHe2, p.ne, T, ν)
        ab += x1 + x2 + x3; em += y1 + y2 + y3
    end
    ab = max(ab, 1e-300)
    (ab, s, em / ab)
end

"Continuum on grid ν: L_ν edge, L_ν through r_ph, r(τ_ν = ⅔), J_ν per point (matrix np × nν)."
function continuum(atm, R, op, ν, rph)
    np = length(atm.r); nν = length(ν); ps = [ray.p for ray in R]
    Le = zeros(nν); Lp = fill(NaN, nν); rτ = fill(NaN, nν); Jm = zeros(np, nν)
    ip = (isfinite(rph) && atm.r[1] < rph < atm.r[end]) ? searchsortedlast(atm.r, rph) : 0
    Threads.@threads for i in 1:nν
        L, J, _, m = SP.monochromatic(atm, op, R, ν[i]; maxit = 60, moments = true, maxit_direct = true)
        Le[i] = L; Jm[:, i] .= max.(J, 0.0)
        χ = m.χ; tν = zeros(np)
        for k in np-1:-1:1; tν[k] = tν[k+1] + SP.dtau(χ[k], χ[k+1], atm.r[k+1] - atm.r[k], atm.jump[k]); end
        k3 = findlast(>(2 / 3), tν)
        if k3 !== nothing && k3 < np
            w = (tν[k3] - 2 / 3) / (tν[k3] - tν[k3+1]); rτ[i] = exp((1 - w) * log(atm.r[k3]) + w * log(atm.r[k3+1]))
        end
        if ip > 0
            w = (rph - atm.r[ip]) / (atm.r[ip+1] - atm.r[ip])
            Lp[i] = 16π^2 * ((1 - w) * atm.r[ip]^2 * m.H[ip] + w * atm.r[ip+1]^2 * m.H[ip+1])
        end
    end
    (; Le, Lp, rτ, J = Jm)
end

# ---------------------------------------------------------------- lines
"Line strength (κ per unit φ [cm⁻¹ Hz]) and emissivity per unit φ [erg cm⁻³ s⁻¹ sr⁻¹] of line ln at point k."
function strength(ln, k, pts, pop, h2pop)
    p = pts[k]
    if ln.atom === :H2
        x = h2pop(p.zone); x === nothing && return (0.0, 0.0)
        nl = p.nH2 * x[ln.l]; nu = p.nH2 * x[ln.u]
    elseif ln.atom === :HD
        (p.nHD > 0 && p.T > 10) || return (0.0, 0.0)
        fr_ = hdfrac(p.T); nl = p.nHD * fr_[ln.l]; nu = p.nHD * fr_[ln.u]
    else
        n = ln.atom === :H ? pop.nH[k] : ln.atom === :He1 ? pop.n1[k] : pop.n2[k]
        nl = n[ln.l]; nu = n[ln.u]
    end
    K = MA.BCOEF * ln.f * (nl - ln.gl / ln.gu * nu)
    (max(K, 0.0), h * ln.ν0 / (4π) * nu * ln.A)          # (stimulated emission larger than absorption: maser, taken as thin)
end
"Voigt/Stark profile of line ln at point k (Z0Physics.line_profile, normalised per Hz)."
function profile(ln, p)
    T = p.T; D = ln.ν0 / c * sqrt(2kB * T / ln.m + XI^2)
    if ln.atom === :H2 || ln.atom === :HD
        return Z0.line_profile(ln.ν0, 1.0, D, ln.A / (4π), 0.0, Float64[], Float64[], 0.0, 0.0)
    end
    γr = 0.6670 * ln.ν0^2 / c^2 * ln.f * (ln.gl / ln.gu) / (4π) * 2          # (as Z0Physics.line_list)
    if ln.atom === :He1
        γe = p.ne > 0 ? Z0.impact_width(ln.nsl, ln.nsu, 1, p.ne, T) : 0.0
        return Z0.line_profile(ln.ν0, 1.0, D, γr, γe, Float64[], Float64[], 0.0, 0.0)
    end
    Z = ln.Z; e = 4.80320471e-10; a0 = 5.29177211e-9
    Fi = 2.6031 * e * max(p.nHp + p.nHe1 + 2.828 * p.nHe2, 0.0)^(2 / 3); Fe = 2.6031 * e * max(p.ne, 0.0)^(2 / 3)
    Ci = e * a0 * Fi / (Z * h); Ce = e * a0 * Fe / (Z * h)
    γe = p.ne > 0 ? Z0.impact_width(ln.l, ln.u, Z, p.ne, T) : 0.0
    pat = Z0.stark_pattern(ln.l, ln.u)
    Z0.line_profile(ln.ν0, 1.0, D, γr, γe, pat.X, pat.I, Ci, Ce)
end

"""
Observer-frame formal solution of one window: L_ν on the observed frequencies νo for the lines `lw` (indices into LL) with point
data (χc, ηc at the window's frequency), strengths K, E (np × nlines) and profiles prof[k][j]. `lines = false`: continuum only.

    window_L(atm, R, pts, νo, χc, ηc, K, E, prof, lw, LL; lines = true, rph = nothing)

L_ν leaving the outer edge at the observed frequencies νo (observer frame, tangent rays). With `rph` given (a number) the result is the
pair (L_edge, L_ph); L_ph is NaN where rph is not finite or not inside the solved region. With rph inside the solved region it is
L_ν = 16π² r² H_ν through the sphere r = rph from the same transfer, returned as (L_edge, L_ph): H_ν = ½∫(I⁺ − I⁻)μ dμ with
μ dμ = p dp / r², so L = 8π² ∫₀^rph (I⁺ − I⁻) p dp. I⁺ (outward) is the ray's intensity where it crosses rph outward on the near side;
I⁻ (inward) where a ray crosses rph inward on the far side, and for core rays (no far side) from an inward integration of the near side
from the edge (no incident light) down to rph, with the velocity projected on that direction. Intensities at the crossing are
interpolated linearly in r between the bracketing points; at p = rph (tangent) I⁺ = I⁻, so the integrand vanishes there.
"""
function window_L(atm, R, pts, νo, χc, ηc, K, E, prof, lw, LL; lines = true, rph = nothing)
    np = length(pts); ps = [ray.p for ray in R]; Lν = zeros(length(νo))
    doph = rph !== nothing && isfinite(rph) && atm.r[1] < rph < atm.r[end]; Lph = fill(NaN, length(νo))
    vth = [sqrt(2kB * p.T / mH + XI^2) for p in pts]          # (H thermal width: the sub-step rule, conservative for heavier lines)
    @inline function coef(k, νx, vz)
        χ = χc[k]; η = ηc[k]
        if lines
            for (jj, j) in enumerate(lw)
                Kk = K[k, jj]; Ek = E[k, jj]; (Kk > 0 || Ek > 0) || continue
                φ = Z0.prof_eval(prof[k][jj]::Z0.LineProf, νx * (1 - vz / c) - LL[j].ν0)
                χ += Kk * φ; η += Ek * φ
            end
        end
        χ, η
    end
    Threads.@threads for iν in eachindex(νo)
        νx = νo[iν]; Io = zeros(length(R)); Ip = zeros(length(R)); Im = zeros(length(R))
        # integrate I along path (point index, signed z) in the direction of +z (sgn = 1) or −z (sgn = −1); with jr > 0 record the
        # intensities where the path crosses r = rph outward (Ip) and inward (Im)
        function walk(path, I, sgn, jr)
            for s in 1:length(path)-1
                ka, za = path[s]; kb, zb = path[s+1]; Δz = zb - za; Δz == 0 && continue
                Ia = I
                vza = sgn * pts[ka].v * za / atm.r[ka]; vzb = sgn * pts[kb].v * zb / atm.r[kb]
                m = lines ? clamp(ceil(Int, abs(vzb - vza) / (SUBV * min(vth[ka], vth[kb]))), 1, 400) : 1
                jmp = (ka == kb + 1 || kb == ka + 1) && atm.jump[min(ka, kb)]
                χ0a, η0a = coef(ka, νx, vza)
                χp = χ0a; ηp = η0a
                for q in 1:m
                    w = q / m; vz = (1 - w) * vza + w * vzb
                    χa_, ηa_ = coef(ka, νx, vz); χb_, ηb_ = coef(kb, νx, vz)
                    χn = (1 - w) * χa_ + w * χb_; ηn = (1 - w) * ηa_ + w * ηb_
                    dτ = jmp ? min(χp, χn) * abs(Δz) / m : 0.5 * (χp + χn) * abs(Δz) / m
                    I = SP.segment(I, ηp / χp, ηn / χn, dτ)
                    χp = χn; ηp = ηn
                end
                if jr > 0
                    ra = atm.r[ka]; rb = atm.r[kb]
                    if ra < rph <= rb                     # moving outward through rph
                        w = (rph - ra) / (rb - ra); Ip[jr] = (1 - w) * Ia + w * I
                    elseif ra >= rph > rb                 # moving inward through rph
                        w = (ra - rph) / (ra - rb); Im[jr] = (1 - w) * Ia + w * I
                    end
                end
            end
            I
        end
        for (jr, ray) in enumerate(R)
            idx = ray.idx; z = ray.z; nr = length(idx)
            # path from the far end to the observer: (point index, signed z)
            path = ray.core ? [(idx[m], z[m]) for m in nr:-1:1] : vcat([(idx[m], -z[m]) for m in 1:nr], [(idx[m], z[m]) for m in nr:-1:1])
            I0 = ray.core ? ηc[idx[nr]] / χc[idx[nr]] : 0.0
            rec = doph && ray.p < rph ? jr : 0
            Io[jr] = walk(path, I0, 1, rec)
            # core rays have no far side: the inward intensity at rph from the near side, edge inwards, direction −z
            rec > 0 && ray.core && walk([(idx[m], z[m]) for m in 1:nr], 0.0, -1, rec)
        end
        L = 0.0; for j in 1:length(R)-1; L += 0.5 * (Io[j] * ps[j] + Io[j+1] * ps[j+1]) * (ps[j+1] - ps[j]); end
        Lν[iν] = 8π^2 * L
        if doph
            Lq = 0.0; pl = 0.0; fl = 0.0
            for j in eachindex(R)
                ps[j] < rph || break
                fj = (Ip[j] - Im[j]) * ps[j]; j > 1 && (Lq += 0.5 * (fl + fj) * (ps[j] - pl)); pl = ps[j]; fl = fj
            end
            Lq += 0.5 * fl * (rph - pl)          # to p = rph, where I⁺ = I⁻
            Lph[iν] = 8π^2 * Lq
        end
    end
    rph === nothing ? Lν : (Lν, Lph)
end

"""
    spectrum(x, d, τ, kph; h2pop = k -> nothing, HL = nothing) -> NamedTuple

Everything of one snapshot (see the module docstring). `d`, `τ`, `kph`: snapspec's densities, grey τ at interfaces and photosphere
interface; `h2pop(point)` the H₂ level fractions at a point (or nothing), `HL` the H2Lines list.
"""
function spectrum(x, d, τ, kph; h2pop = k -> nothing, HL = nothing, νgrid = nugrid())
    t0 = time(); CO.tables()
    function build(fT)
        A = points(x, d, τ; fT); pts = A.P
        atm = (r = [p.r for p in pts], ρ = [p.ρ for p in pts], T = [p.T for p in pts], zone = [p.zone for p in pts], jump = A.jump)
        _, R = SP.rays(atm.r; ncore = A.kc > 1 ? 8 : 1)
        A.kc > 1 || (R = [merge(ray, (core = false,)) for ray in R])
        (A, pts, atm, SP.weights!(atm, R))
    end
    Teff = x.Teff
    Jdilf(pts, r0) = k -> (p = pts[k]; p.τg > 1 ? (ν -> MA.Bν(ν, p.T)) :
                           (W = 0.5 * (1 - sqrt(max(1 - (r0 / p.r)^2, 0.0))); ν -> W * MA.Bν(ν, max(Teff, 1.0))))
    ∫(L, ν) = sum(0.5 * (L[i] + L[i+1]) * (ν[i+1] - ν[i]) for i in 1:length(ν)-1 if isfinite(L[i]) && isfinite(L[i+1]); init = 0.0)
    fT = 1.0; A, pts, atm, R = build(fT); fhist = Float64[]
    # refined photosphere: normalise its T⁴ so that the continuum carries the run's luminosity through the refined zone's outer interface
    rtop = A.kr > 0 ? x.r[A.kr+1] : NaN; Ltop = A.kr > 0 ? 4π * rtop^2 * x.fr[A.kr+1] : NaN
    if A.kr > 0 && isfinite(Ltop) && Ltop > 0
        νc = νgrid[1:3:end]
        for it in 1:5
            pp = populations(pts, Jdilf(pts, rtop))
            cc = continuum(atm, R, FOpac(pts, pp), νc, rtop); L = ∫(cc.Lp, νc); push!(fhist, fT)
            (isfinite(L) && L > 0) || break
            abs(L / Ltop - 1) < 0.03 && break
            fT *= (Ltop / L); A, pts, atm, R = build(fT)
        end
    end
    np = length(pts)
    rph = kph === nothing ? NaN : x.r[kph]
    if A.nref > 0                                   # refined: where the sub-points' grey τ crosses ⅔
        i = findlast(p -> p.τg > 2 / 3, pts)
        i !== nothing && i < np && (w = (pts[i].τg - 2 / 3) / (pts[i].τg - pts[i+1].τg); rph = (1 - w) * pts[i].r + w * pts[i+1].r)
    end
    # pass 1: the run's radiation for the populations
    r0 = isfinite(rph) ? rph : atm.r[1]
    pop1 = populations(pts, Jdilf(pts, r0)); t1 = time()
    c1 = continuum(atm, R, FOpac(pts, pop1), νgrid, rph); t2 = time()
    # pass 2: the continuum J of pass 1
    pop = populations(pts, k -> Jinterp(νgrid, c1.J[k, :])); t3 = time()
    op = FOpac(pts, pop)
    cont = continuum(atm, R, op, νgrid, rph); t4 = time()
    # lines
    LL = linelist(HL); nl = length(LL)
    nv = 2 * round(Int, VWIN / DV) + 1; vgrid = collect(range(-VWIN, VWIN, length = nv))
    prof = [Vector{Union{Nothing,Z0.LineProf}}(nothing, nl) for _ in 1:np]; Ks = zeros(np, nl); Es = zeros(np, nl)
    Threads.@threads for k in 1:np
        for j in 1:nl
            Ks[k, j], Es[k, j] = strength(LL[j], k, pts, pop, h2pop)
            prof[k][j] = (Ks[k, j] > 0 || Es[k, j] > 0) ? profile(LL[j], pts[k]) : nothing
        end
    end
    t5 = time()
    Lprof = zeros(nl, nv); Lcont = zeros(nl); Lline = zeros(nl); EW = zeros(nl); τr = zeros(nl, NRT)
    # the same at the photosphere sphere r = rph (L_ν = 16π² r² H_ν through it, observer frame; NaN without a photosphere)
    LprofP = fill(NaN, nl, nv); LcontP = fill(NaN, nl); LlineP = fill(NaN, nl); EWP = fill(NaN, nl)
    rgrid = exp.(range(log(atm.r[1]), log(atm.r[end]), length = NRT))
    Jc = [Jinterp(νgrid, cont.J[k, :]) for k in 1:np]
    for j in 1:nl
        ν0 = LL[j].ν0; νo = ν0 .* (1 .- vgrid ./ c)                      # observed: v > 0 = redshift (away)
        lw = [i for i in 1:nl if abs(LL[i].ν0 / ν0 - 1) < 1.5 * VWIN / c]
        χc = zeros(np); ηc = zeros(np)
        for k in 1:np
            ab, sc, St = SP.coeffs3(op, k, ν0); χc[k] = ab + sc; ηc[k] = ab * St + sc * Jc[k](ν0)
        end
        Kw = Ks[:, lw]; Ew = Es[:, lw]; pw = [[prof[k][i] for i in lw] for k in 1:np]
        Lp, Lpp = window_L(atm, R, pts, νo, χc, ηc, Kw, Ew, pw, lw, LL; rph)
        Lcv, Lcvp = window_L(atm, R, pts, [ν0], χc, ηc, zero(Kw), zero(Ew), pw, lw, LL; rph)      # (same sub-segments as the profile)
        Lc = Lcv[1]
        Lprof[j, :] .= Lp; Lcont[j] = Lc
        dν = ν0 * DV / c; Lline[j] = sum(Lp .- Lc) * dν; EW[j] = -Lline[j] / max(Lc, 1e-300) * (c / ν0^2) * 1e8     # Å; < 0 = emission
        if isfinite(Lcvp[1])
            LprofP[j, :] .= Lpp; LcontP[j] = Lcvp[1]
            LlineP[j] = sum(Lpp .- Lcvp[1]) * dν; EWP[j] = -LlineP[j] / (abs(Lcvp[1]) > 0 ? Lcvp[1] : 1e-300) * (c / ν0^2) * 1e8
        end
        # radial line-centre optical depth (comoving, static), from the outer edge inwards
        tk = zeros(np)
        for k in np-1:-1:1
            φa = prof[k][j] === nothing ? 0.0 : Z0.prof_eval(prof[k][j], 0.0); φb = prof[k+1][j] === nothing ? 0.0 : Z0.prof_eval(prof[k+1][j], 0.0)
            tk[k] = tk[k+1] + 0.5 * (Ks[k, j] * φa + Ks[k+1, j] * φb) * (atm.r[k+1] - atm.r[k])
        end
        for (i, rr) in enumerate(rgrid)
            k = clamp(searchsortedlast(atm.r, rr), 1, np - 1); w = (rr - atm.r[k]) / (atm.r[k+1] - atm.r[k])
            τr[j, i] = exp((1 - w) * log(max(tk[k], 1e-30)) + w * log(max(tk[k+1], 1e-30)))
        end
    end
    t6 = time()
    (; step = x.step, t = x.t, ν = νgrid, cont, np, kc = A.kc, fT, fhist, rtop, Ltop, nrefined = count(p -> p.lte && p.τg <= TDEEP, pts), lines = LL,
       vgrid, Lprof, Lcont, Lline, EW, LprofP, LcontP, LlineP, EWP, rgrid, τr, rph, Teff, Lrun_edge = 4π * x.r[end-2]^2 * x.fr[end-2],
       Lrun_ph = kph === nothing ? NaN : 4π * x.r[kph]^2 * x.fr[kph],
       times = (points = t1 - t0, cont1 = t2 - t1, pop2 = t3 - t2, cont2 = t4 - t3, profiles = t5 - t4, windows = t6 - t5))
end

end # module
