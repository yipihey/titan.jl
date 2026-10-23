# S7c-2: trapped Lyman-series radiation force on the neutral gas of the grid (included by RadTransfer.jl; module RadTransfer).
#
# McKee & Tan (2008, ApJ 681, 771) §4 and Appendices B, C: in an opaque neutral medium the Lyα radiation is nearly isotropic, P_iso = 4πJ/(3c), and the radiative force per unit
# volume is −∇P_iso = ρ κ F/c (pure scattering, ∇·F = 0; C4). The scattered photons stay in the net flux (spherical symmetry: nothing returns to the small star), so the force of a thick
# scattering layer is κ F/c per unit mass with the *conserved* flux F = L_ν e^{−τ_abs}/(4π r²), not the saturating (1 − e^{−τ}) F/c of a beam from which the scattered photons are removed
# (what `transfer`'s direct beam does). Photons injected in the Doppler core and the diffusion regime |x̂| < 1 of a line escape in the damping wings after frequency diffusion (Adams 1972;
# Harrington 1973; Neufeld 1990), which bounds J: J/F = 0.411 (a τ̄)^{1/3} / (min[1, (aτ̄/450)^{1/3}] + Γ) (B7–B8, B13–B14), τ̄ = σ̄ N_eff; far in the wings the scattering is coherent,
# J = (3/4π) τ_ν F (B10). MT08's B11 is the harmonic blend of the two; here the coherent limit uses the full Kramers–Heisenberg cross section (all Lyman lines, the Rayleigh continuum) and
# each frequency takes the diffusion plateau of the line nearest in Doppler widths (`LymanWings.jratio`). P_iso = (4π/3c) Σ_bins min[B*_ν Δν, (J/F) F_bin] (B22), the spherical column enters as
# the slab-equivalent N_eff = r_j² Σ_{k≥j} n_k Δr_k/r_k² (C6–C10 without the power-law assumption: (2/5) n r for ρ ∝ r^{-3/2}, exact for any profile and the finite edge); the force on a zone is
# −ΔP_iso/(ρ Δr) (P_iso = 0 at the outer edge, C4).
#
# The component stored in `ARAD_LYA` is the *excess* over what the direct beam already applies for the same opacity: −ΔP_iso/(ρΔr) minus the beam's single-scattering force from the
# neutral-H scattering (`arad_hi` of `transfer`); where the medium is thin the two agree and the component vanishes; with `lya` the total force of the neutral-H scattering is MT08's.
# Line sources (not in the beam): (i) `Linject`, the line-centre Lyα luminosity entering at the inner boundary (the H II region inside r_s), (ii) Lyα produced in the grid zones by the stellar
# photoionisation (case-B recombination cascades, 0.68 per recombination, limited by the photoionisation rate; thermal collisional Lyα is not a source), half of it assumed to be
# emitted outwards (a source inside the cloud: the rest escapes inwards through the lower column); they enter at line centre with the Lyα plateau and the two-photon destruction.
#
# Calibration: a thin shell has M_F ≡ F_rad c/L_α = (4π/3) J/F = 14 / 29 / 63 for N_H = 1e19 / 1e20 / 1e21 cm⁻² (Δv_D = 12.9 km/s) against M_F ≈ 10 N19^0.44 = 10 / 27.5 / 76 of
# the static-shell Monte-Carlo results of Dijkstra & Loeb (2008, Fig. 6) quoted by Smith, Bromm & Loeb (2017), and t_trap/t_light ≈ 15 (τ_0/10^5.5)^{1/3} (Adams 1975).
# Validity: aτ̄ ≳ 1 (N_eff ≳ 1e16 Δv₆² cm⁻²), pure scattering (no dust; two-photon destruction included, collisional de-excitation ignored), velocity range over the column small compared
# with the escape width (aτ̄)^{1/3} b (hundreds of km/s here: yes for the infall), spherical symmetry (the 2D flashlight effect, a polar cavity, is absent by construction), the neutral
# column of the grid only (the trapping inside r_s, the free-fall column between the star and the first interface, is not on the grid).

const E_LYA = LymanWings.hP * LymanWings.ν_line(2)       # 10.2 eV [erg]

"""
Lyα photon luminosity density [erg cm⁻³ s⁻¹] of zone k that is driven by the stellar radiation: case-B recombination cascades (0.68 Lyα per recombination) of the photoionised gas, limited by
the ground-level photoionisation rate n_H Γ_H of the zone (ionising photons > 13.6 eV; the Balmer-continuum ionisation of the thermally excited n = 2 atoms of a hot zone is a thermal channel, not a driver) (in photoionisation equilibrium the two are equal; in gas ionised by collisions, e.g. a shock-heated first zone, there is no stellar driver and the
term is zero). Collisional excitation of H is *not* a source: it radiates thermal energy of the zone (its rate at n_H = 1e10, T = 6e3 K is 1e49 erg/s per zone — a transient of the
shocked first zone — and the photons are thermalised in that dense gas, ε n_e ≫ n_e,cr), not stellar luminosity.
"""
function lya_production(med::Medium, k, GammaH)
    T = med.T[k]; nH = med.nH[k]
    ne = nH * xs(med, :e, k); np = nH * xs(med, :Hp, k); n1 = nH * xs(med, :H, k)
    rec = 2.59e-13 * (T / 1e4)^-0.7 * ne * np
    0.68 * min(rec, n1 * GammaH) * E_LYA
end

"""
    lya_pressures(g, tab, La, med, o; Linject = 0, Rstar = 0, Lsrc = nothing, TinjK = 1.5e4) -> (P, parts)

Isotropic pressure P_iso [dyn cm⁻²] at the n + 1 interfaces from the stellar spectrum (`La`: (n+1) × nf luminosity per bin attenuated by true absorption only), the injected line and the grid
sources; `parts` = (stellar, injected, grid pressures, N_eff, k_ρ, T). `tab` is the `LymanWings.Table` of the grid.
"""
function lya_pressures(g::FreqGrid, tab, La, med::Medium, o::Opts; Linject = 0.0, Rstar = 0.0, Lsrc = nothing, TinjK = 1.5e4, GammaH = nothing, skipbins = nothing, nolines = false)
    LW = LymanWings
    n = length(med); r = med.r; nf = nfreq(g); Δr = diff(r)
    nHI = [med.nH[k] * xs(med, :H, k) for k in 1:n]
    Nout = zeros(n + 1); Tw = zeros(n + 1); Ew = zeros(n + 1)
    for k in n:-1:1
        dN = nHI[k] * Δr[k]
        Nout[k] = Nout[k+1] + dN; Tw[k] = Tw[k+1] + dN * med.T[k]; Ew[k] = Ew[k+1] + dN * med.nH[k] * xs(med, :e, k)
    end
    for j in 1:n
        if Nout[j] > 0; Tw[j] /= Nout[j]; Ew[j] /= Nout[j]; else Tw[j] = med.T[j]; Ew[j] = 0.0; end
    end
    rc = [r[1] > 0 ? sqrt(r[k] * r[k+1]) : 0.5 * (r[k] + r[k+1]) for k in 1:n]
    # slab-equivalent column (MT08 C6–C10 without the power-law assumption): the coherent-scattering pressure P(r_j) = (1/c) Σ_{k≥j} κ_k Δr_k F(r_k) with F ∝ r⁻², written as
    # τ_eff F(r_j)/c, i.e. N_eff = r_j² Σ_{k≥j} n_k Δr_k / r_k² (= (2/5) n r for ρ ∝ r^{-3/2} to infinity, MT08's value; exact for any profile and a finite edge)
    Neffs = zeros(n + 1)
    for k in n:-1:1; Neffs[k] = Neffs[k+1] + nHI[k] * Δr[k] / rc[k]^2; end
    kρ = fill(1.5, n + 1)                     # (diagnostic only: the local density slope)
    if n >= 2
        for j in 1:n
            ka, kb = j == 1 ? (1, 2) : (j - 1, j)
            s = log(max(med.nH[ka], 1e-300) / max(med.nH[kb], 1e-300)) / log(rc[kb] / rc[ka]); kρ[j] = isfinite(s) ? s : 1.5
        end
    end
    Lcum = zeros(n + 1)          # Σ_{k<j} ½ L_k of the grid sources
    for k in 1:n; Lcum[k+1] = Lcum[k] + 0.5 * lya_production(med, k, GammaH === nothing ? 0.0 : GammaH[k]) * volume(r, k); end
    Pst = zeros(n + 1); Pln = zeros(n + 1); Pgr = zeros(n + 1); Neffv = zeros(n + 1)
    νa = LW.ν_line(2); vt = o.lya_vt
    Threads.@threads for j in 1:n
        Nout[j] > 0 || continue
        Neff = r[j]^2 * Neffs[j]; Neffv[j] = Neff
        T = Tw[j]
        eps = Ew[j] > 0 ? 4.4e-9 / (1 + 1.55e4 / Ew[j]) : 0.0
        A4 = 4π * r[j]^2
        # line sources: the line-centre ratio J/F is the plateau of Lyα (σ → the damping core: coherent J ≫ plateau)
        Jl = LW.line_plateau(2, Neff, T; vt, eps)
        if Jl[1] > 0 && !nolines
            D = LW.Δν_D(2, T; vt)
            for (Fl, TB, dst) in ((Linject, TinjK, Pln), (Lcum[j], max(T, TinjK), Pgr))
                Fl > 0 || continue
                F = Fl / A4
                Bcap = 2 * Jl[2] * D * Bν(νa, TB)
                dst[j] = 4π / (3c) * min(Jl[1] * F, Bcap)
            end
        end
        # the stellar spectrum, bin by bin (structured bins: sub-sampled)
        s = 0.0
        for f in 1:nf
            skipbins !== nothing && skipbins[f] && continue
            Fb = La[j, f] / A4
            Fb > 0 || continue
            Jr = 0.0
            if tab.structured[f]
                νq, wq, σq = LW.sigma_sub(tab, f, T)
                for q in eachindex(νq); Jr += wq[q] * LW.jratio(νq[q], σq[q], Neff, T; vt, eps); end
            else
                Jr = LW.jratio(g.ν[f], LW.sigma_bin(tab, f, T), Neff, T; vt, eps)
            end
            jf = Jr * Fb
            if Rstar > 0 && Lsrc !== nothing
                jf = min(jf, Lsrc[f] / (4π^2 * Rstar^2))          # J ≤ B*_ν Δν (B22)
            end
            s += jf
        end
        Pst[j] = 4π / (3c) * s
    end
    (Pst .+ Pln .+ Pgr, (stellar = Pst, injected = Pln, grid = Pgr, Neff = Neffv, kρ = kρ, T = Tw))
end

"Radial acceleration [cm s⁻²] of each zone from the trapped radiation pressure, a_k = max((P_k − P_{k+1})/(ρ_k Δr_k), 0) (the total force of the neutral-H scattering; `transfer` subtracts the beam's share)."
function lya_accel(g::FreqGrid, tab, La, med::Medium, rho, o::Opts; Linject = 0.0, Rstar = 0.0, Lsrc = nothing, GammaH = nothing, skipbins = nothing, nolines = false)
    n = length(med); r = med.r
    P, = lya_pressures(g, tab, La, med, o; Linject, Rstar, Lsrc, GammaH, skipbins, nolines)
    [max((P[k] - P[k+1]) / (rho[k] * (r[k+1] - r[k])), 0.0) for k in 1:n]      # (a pressure that rises outward, from a jump of the density slope, is not a force on the gas)
end

"Mask of the frequency bins within |ν/ν_α − 1| < win of Lyα (S8a: the part of the spectrum that the Monte Carlo transports)."
lya_window(g::FreqGrid, win) = BitVector([abs(g.ν[f] / LymanWings.ν_line(2) - 1) < win for f in 1:nfreq(g)])
