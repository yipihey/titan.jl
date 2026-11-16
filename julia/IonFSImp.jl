# Included in PopIIIMaterial after IonFS.jl. IONFSIMP=1: the IonFS ray (ionizing groups and the Lyman–Werner band) as unknowns of the
# hydro Newton (TitanResidFast.RAYIMP). Nothing below runs, and no code path changes, when IONFSIMP is off. See docs/POPIII_IONFS_IMPLICIT.md.
"""
IONFSIMP=1 (default off; needs IONFS=1, COUPLED=1, STARLIGHT=1): the formal-solution-sourced ray of IonFS, implicit in the Newton solve.

Groups (4, `IFSI_EDGES` and the LW band): g = 1 H I band 13.598–24.587 eV (absorber H), 2 He I band 24.587–54.418 eV (H, He), 3 He II band
54.418–300 eV (H, He, He⁺), 4 Lyman–Werner band 11.2–13.598 eV (H₂ lines, with the H I Lyman-line shielding). The edges are the absorbers'
thresholds, so the set of absorbers is fixed inside every group; the group cross sections are photon-weighted over the formal solution's own
spectrum (`IONFSNSUB` frequencies per group). Broad groups lose the hardening inside a band; a front zone (τ_912 ~ 1e3 per zone) absorbs all
of a group anyway, so the number of ionizations per zone (photon conservation) is unaffected, only how the last few photons are shared.

Unknowns: y_g(k) = ln N_g, the photon rate [s⁻¹] of group g entering zone k (TitanResidFast rows 10–13). Ray row of zone k > k₀:
y_g(k) = y_g(k−1) − Δτ_g(k−1), with Δτ_g of zone k−1 a function of its iterate (ρ, T, x_H⁺, x_H₂ and through the CK expansion He, He⁺),
and y_g(k₀) = ln N_g,src. In coupled zone k the photo-rates, the photo-heating and (LW) the H₂ dissociation depend on y_g(k) and on the
zone's own composition (photon-conserving: absorbed N_g (1 − e^{−Δτ_g}), shared by Δτ_g,s/Δτ_g), all with exact derivatives in the
Newton (the zone's local Jacobian is 8 wide: ρ, T, x_H⁺, x_H₂, y_1…y_4). He I/He II photo-rates per atom are N_g fτ(Δτ_g) σ_g,s Δr/V,
fτ = (1 − e^{−τ})/τ, with Δτ_g from the zone's own He stages: a fixed point u = expansion(rates(u)) in (x_H, x_He, x_He⁺), solved by
Newton inside the residual and differentiated by the implicit-function theorem (`ifsi_ufix`).

Lyman–Werner (photon-conserving self-shielding consistent with SLSHIELD's fits): the dissociation rate per molecule at H₂ column N is
k(N) = σ_d φ f_sh(N, T) f_HI(N_HI) (Wolcott-Green, Haiman & Bryan 2011 eq. 10 with α = 1.1, and eq. 13 of Wolcott-Green & Haiman 2011:
`fsh_h2`, `fsh_hi`), φ the unattenuated LW photon flux, σ_d = 1.1e8 L_ν(12.87 eV)/N_LW from the k27 normalisation of Abel et al.
(1997) (k = 1.1e8 F_ν(12.87 eV) s⁻¹, as STARLIGHT) and the formal solution's spectrum. Draine & Bertoldi (1996, eq. 38) relate f_sh to the
band's absorption: the pumping (absorptions) per area down to column N is ζ_pump(0)/σ ∫₀ᴺ f_sh dN', and ⟨p_diss⟩ ≈ 0.15 of the
absorptions dissociate (DB96 §4, "⟨p_diss⟩ ≈ 0.15 averaged over all UV absorptions"; `IONFSPD`). So the LW unknown carries the photons
not yet absorbed by H₂, and the photons used so far fix the equivalent static column N₁ of the zone's inner face:
  N_LW(k) = Φ [1 − (σ_d/p) G(N₁)],   G(N) = ∫₀ᴺ f_sh dN' (closed form, `ifsi_G`).
Zone k (ΔN = x_H₂ n_H Δr) absorbs A = Φ (σ_d/p) [G(N₁ + ΔN) − G(N₁)] carrier photons (Sobolev shielding, SLSHIELD=sobolev: the zone's mean
f is at least f_sh(N_sob), as STARLIGHT's min(N_static, N_sob)), written as Δτ_LW = A/N_LW(k) so that the row is the same as the
ionizing groups' and a zone can never take more than enters it; dissociations D = p N_LW (1 − e^{−Δτ_LW}) f_HI(N_HI). For a static column
this reproduces STARLIGHT's rate σ_d φ f_sh(N) f_HI exactly; p cancels except where the carrier runs out (it then stops the dissociation).
N_HI (Lyman-line shielding) stays the start-of-step column, as in STARLIGHT; f_HI does not deplete the carrier.

Rates above the edges come only from the ray: in every LOCALJ zone outside the grey photosphere (k ≥ k₀; `IONFSLJ=all`: every LOCALJ
zone) the diluted blackbody keeps only the levels below 13.598 eV (`ljcut`), with case-B He recombination (`caseB`); the LW rate of
STARLIGHT (W B_ν(T_c) at 12.87 eV) is replaced by the ray's in every zone k ≥ k₀ (beyond the ray's last zone: its start-of-step sweep).
The energy: photo-heating Σ_g A_g E_g/V in the coupled zone's energy row (implicit), its photon energy leaving the grey field through
RAYQ after the step (as IonFS); the LW heat (4.48 + 0.4 eV per dissociation) through GASQ, as STARLIGHT.
Diagnostics: <RUNREC>/ionfsi.csv per step (photon bookkeeping per group, LW dissociations of the ray and of the SLSHIELD path, He²⁺ inside
3e3 R☉, the Newton's ray rows).
"""
const IFSI_EDGES = [13.598, 24.587, 54.418, 300.0]
const IFSI_LW = (11.2, 13.598)
const IFSI_SUB = get(ENV, "IONFSSUB", "1") == "1"                     # sub-zone photosphere in the formal solution (ifsi_source!)
const IFSI_TAUC = parse(Float64, get(ENV, "IONFSTAUC", "3e3"))       # the formal solution starts in the innermost zone with grey τ above it
const IFSI_NZMIN = parse(Int, get(ENV, "IONFSNZMIN", "6"))           # … and at least this many zones inside r_b
const IFSI_PD = parse(Float64, get(ENV, "IONFSPD", "0.15"))          # ⟨p_diss⟩, Draine & Bertoldi (1996) §4
>>>

const IFSI_NSUB = parse(Int, get(ENV, "IONFSNSUB", "8"))
const IFSI_LJ = get(ENV, "IONFSLJ", "out")
const IFSI_LWMIN = parse(Float64, get(ENV, "IONFSLWMIN", "1e-2"))
const IFSI_LWON = get(ENV, "IONFSLW", "1") == "1"
const IFSI_K27 = 1.1e8                                                    # AAZN97 k27 (as STARLIGHT's kLW)
const IFSI_E12 = 12.87                                                    # eV (STARLIGHT's νLW)
const IFSI_TMAX = 700.0                                                   # Δτ of a ray row is capped (N e^{−700} = 0): bounded ln N
const IFSI_HECE = get(ENV, "IONFSHECE", "1") == "1"                      # He⁺ collisional-excitation cooling in the ray zones (ifsi_herec)
mutable struct IFIState
    ysrc::Vector{Float64}; Phi::Float64; σd::Float64; ELW::Float64
    fhi::Vector{Float64}; sob::Vector{Float64}
    kLW::Vector{Float64}; kLWsl::Vector{Float64}; D::Vector{Float64}; Dsl::Vector{Float64}
    kend::Int; emit::Vector{Float64}; absb::Vector{Float64}; esc::Vector{Float64}; rowerr::Float64; niter::Int
    phz::Vector{Any}; lrec::Vector{Float64}
    yint::Matrix{Float64}; qion::Vector{Float64}          # (post) ln N_g at every interface k0…b+1 (NaN inside k0), ionizing heat per zone
end
const IFSI = Ref{Any}(nothing)
ifsi_alloc(mgr) = IFIState(zeros(4), 0.0, 0.0, 0.0, ones(mgr), fill(Inf, mgr), zeros(mgr), zeros(mgr), zeros(mgr), zeros(mgr), 0, zeros(4), zeros(4), zeros(4), 0.0, 0,
                           Any[nothing for _ in 1:mgr], zeros(mgr), fill(NaN, 4, mgr + 1), zeros(mgr))
"""
He radiative-recombination cooling of a ray zone [erg cm⁻³ s⁻¹] at (ρ, T, full composition d), from the network's own recombination rates:
the Milne inverses α of its He I and He II photo-processes (`ph`: tuples (process, R, α), case B in the ray's env), × (χ + 0.75 kT) per
recombination. ChemNet keeps the binding energies in u: a recombination gives χ back to the gas, so without this term a photoionization–
recombination cycle heated the gas by the whole hν of the absorbed photon instead of hν − χ − ⟨E_kin⟩ (T → 9e4 K at the He I front,
t_ifs9). The kinetic-energy loss per recombination is the three-level atom's 0.75 kT (tl3_heat). The emitted energy enters the grey field (RAYQ).
"""
function ifsi_herec(ρ, T, d, ph)
    nH = ρ / CNn.μH(); ne = nH * max(d[:e], zero(d[:e])); q = zero(ne) + zero(T)
    # He⁺ collisional excitation (n = 1 → 2, 40.8 eV; IONFSHECE=1, default): ceHeII of Cen (1992) as coded in EmissionKernels.jl
    # (cooling_atomic.jl, the AAZN97 network's coefficient): 5.54e-17 e^{−473638/T} T^−0.397 / (1 + (T/1e5)^½) erg cm³ s⁻¹ × n_e n_He⁺.
    # The network has no He line cooling (ChemRates: "He⁺ two-photon and line cooling — these come with the model atom"); without it the
    # ray's photo-heated zones went to 0.5–1e5 K (t_ifs10, t_ifs14). H's collisional excitation is the three-level atom's (LYABB).
    IFSI_HECE && (q += 5.54e-17 * exp(-min(69.07755278982137, 473638.0 / T)) * T^(-0.397) / (1 + sqrt(T / 1e5)) * ne * nH * max(d[:Hep], zero(d[:Hep])))
    for (p, _, α) in ph
        p isa CNn.Photo || continue
        if p.R === :He
            q += α * ne * nH * max(d[:Hep], zero(d[:Hep])) * (CNn.χHe1 + 0.75CNn.kB * T)
        elseif p.R === :Hep
            q += α * ne * nH * max(d[:Hepp], zero(d[:Hepp])) * (CNn.χHe2 + 0.75CNn.kB * T)
        end
    end
    q
end
"Ray zone k (on, within k₀…k_end) of the current step."
ifsi_on(k) = (st = IFS[]; st !== nothing && st.k0 > 0 && k <= length(st.on) && st.on[k])

# ---------------------------------------------------------------- self-shielding integral
"G(N) = ∫₀ᴺ f_sh(N', T) dN' of `fsh_h2` (WGHB11 eq. 10, α = 1.1), closed form."
function ifsi_G(N, T)
    x = N / 5e14; b5 = sqrt(CNn.kB * T / 1.6735575e-24) / 1e5; a = 8.5e-4
    5e14 * (9.65 * b5 * (1 - (1 + x / b5)^(-0.1)) + (0.07 / a) * (exp(-a) - exp(-a * sqrt(1 + x))))
end
"Column N with G(N, T) = g (bisection in ln N on the values, one Newton correction for the derivatives in g and T)."
function ifsi_ginv(g, T)
    gv = _fv(g); Tv = _fv(T); z = zero(g) + zero(T)
    gv <= 0 && return z
    lo = gv; hi = 1e30
    ifsi_G(hi, Tv) <= gv && return z + hi
    for _ in 1:90
        mid = sqrt(lo * hi); ifsi_G(mid, Tv) < gv ? (lo = mid) : (hi = mid); hi / lo < 1 + 1e-13 && break
    end
    N0 = sqrt(lo * hi)
    N0 - (ifsi_G(N0, T) - g) / fsh_h2(N0, T)
end
"LW in zone k: (Δτ_LW, dissociations D [s⁻¹]) at the carrier y = ln N_LW, T and x_H₂ (generic in the number types)."
function ifsi_lw(I, k, y, T, xH2, nH, dr)
    z = zero(y) + zero(T) + zero(xH2) + zero(nH)
    (IFSI_LWON && I.Phi > 0) || return (z, z)
    ΔN = max(xH2, zero(xH2)) * nH * dr
    q = -expm1(min(y - I.ysrc[4], zero(y)))                      # fraction of the carrier used inside zone k
    N1 = ifsi_ginv(q * IFSI_PD / I.σd, T)
    fst = _fv(ΔN) < 1e-4 * max(_fv(N1), 5e14) ? fsh_h2(N1 + ΔN / 2, T) : (ifsi_G(N1 + ΔN, T) - ifsi_G(N1, T)) / ΔN
    fe = isfinite(I.sob[k]) ? max(fst, fsh_h2(max(xH2, zero(xH2)) * nH * I.sob[k], T)) : fst
    A = I.σd * I.Phi / IFSI_PD * ΔN * fe
    τ = _fv(A) > 0 ? exp(min(log(A) - y, 69.0)) : A * exp(-min(y, 690.0))   # (capped at e^69: the carrier is gone)
    D = IFSI_PD * exp(y) * (-expm1(-τ)) * I.fhi[k]
    (τ + z, D + z)
end

# ---------------------------------------------------------------- source
"Group frequencies: 24 below 11.2 eV (luminosity normalisation), IONFSNSUB per LW and ionizing group, and 12.87 eV (σ_d) last."
function ifsi_freqs()
    νs = Float64[]; dν = Float64[]; grp = Int[]
    Elo = exp.(range(log(0.05), log(IFSI_LW[1]), length = 25))
    for i in 1:24; push!(νs, sqrt(Elo[i] * Elo[i+1]) * EVERG / CNn.hP); push!(dν, (Elo[i+1] - Elo[i]) * EVERG / CNn.hP); push!(grp, 0); end
    bands = vcat([IFSI_LW], [(IFSI_EDGES[g], IFSI_EDGES[g+1]) for g in 1:3]); ids = (4, 1, 2, 3)
    for (b_, gid) in zip(bands, ids)
        e1, e2 = b_; q = (e2 / e1)^(1 / IFSI_NSUB)
        for i in 1:IFSI_NSUB
            push!(νs, e1 * q^(i - 0.5) * EVERG / CNn.hP); push!(dν, e1 * (q^i - q^(i - 1)) * EVERG / CNn.hP); push!(grp, gid)
        end
    end
    push!(νs, IFSI_E12 * EVERG / CNn.hP); push!(dν, 0.0); push!(grp, -1)
    νs, dν, grp
end
const IFSI_FREQ = ifsi_freqs()
"""
IONFSLDIAG=n (diagnostic, default 0): at the first n source solves, why L_fs ≠ L_run. (1) the same formal solution with the run's grey
opacity (χ = chifn ρ, absorption κ_P ρ, LTE B(T) source): if that gives L_run, the method is right and the difference is the opacity;
(2) the Rosseland and Planck means of the source opacity against the run's χ and κ_P at the points of the solved region.
"""
const IFSI_LDIAG = Ref(parse(Int, get(ENV, "IONFSLDIAG", "0")))
struct IFGrey; χ::Vector{Float64}; κ::Vector{Float64}; T::Vector{Float64}; end
function ifsi_ldiag(atm, op, R, νs, dν, grp, L, Lfs, Lrun, zone, k0, tauz)
    SP_ = IFSP(); chi = TM.v(:chifn); kp = TM.v(:xkpn); dn = TM.v(:dn)
    @eval Main.Spectrum.coeffs3(o::$IFGrey, k, ν) = (o.κ[k], max(o.χ[k] - o.κ[k], 0.0), $(IFCO()).Bν(ν, o.T[k]))
    gop = IFGrey([chi[z] * dn[z] for z in zone], [min(kp[z], chi[z]) * dn[z] for z in zone], copy(atm.T))
    Lg = zeros(length(νs))
    Threads.@threads for i in eachindex(νs)
        grp[i] >= 0 && (Lg[i] = max(Base.invokelatest(SP_.monochromatic, atm, gop, R, νs[i]; maxit = 60, maxit_direct = true)[1], 0.0))
    end
    Lgrey = sum(Lg[i] * dν[i] for i in eachindex(νs) if grp[i] >= 0)
    @printf("IONFSLDIAG step %d: L_run %.4e  L_fs %.4e (ratio %.3f)  L_fs with the run's grey opacity %.4e (ratio to L_run %.3f)  [r_b %.4g R_sun, zones %d..%d]\n",
            Int(TF.getscalar(:jstep)), Lrun / 3.83e33, Lfs / 3.83e33, Lrun / Lfs, Lgrey / 3.83e33, Lgrey / Lrun, TM.v(:rn)[k0] / 6.957e10, zone[1], zone[end])
    # means of the source opacity at every solved point, against the run's
    for j in unique([1, length(zone) ÷ 4, length(zone) ÷ 2, (3 * length(zone)) ÷ 4, length(zone) - 2, length(zone)])
        (1 <= j <= length(zone)) || continue
        z = zone[j]; T = atm.T[j]; sR = 0.0; wR = 0.0; sP = 0.0; wP = 0.0
        for i in eachindex(νs)
            grp[i] >= 0 || continue
            ab, sc, _ = ifs_coeffs3(op, j, νs[i]); x = CNn.hP * νs[i] / (CNn.kB * T)
            B = IFCO().Bν(νs[i], T); dB = x > 600 ? 0.0 : B * x * exp(x) / expm1(x) / T
            sR += dB * dν[i] / (ab + sc); wR += dB * dν[i]; sP += ab * B * dν[i]; wP += B * dν[i]
        end
        @printf("IONFSLDIAG   point %d zone %d T %.0f rho %.3e tau_grey %.3g: Rosseland(FS) %.3e Planck(FS) %.3e per cm | run chi %.3e kappa_P %.3e per cm\n",
                j, z, T, atm.ρ[j], tauz[z], wR / sR, sP / wP, chi[z] * dn[z], kp[z] * dn[z])
    end
    nothing
end
"IonFS's source (ifs_source!) for the implicit ray: three ionizing groups in st (N, G), the LW band in IFSI (Φ, σ_d, mean energy)."
function ifsi_source!(st, m, a, b, k0, tauz)
    t0 = time(); dn, tn, rn, fr = TM.v(:dn), TM.v(:tn), TM.v(:rn), TM.v(:frn)
    I = IFSI[]
    kc = a; for k in a:k0-1; tauz[k] > IFSI_TAUC && (kc = k); end
    kc = max(a, min(kc, k0 - IFSI_NZMIN))
    SP_ = IFSP(); atm0 = SP_.atmosphere((r = rn[kc:k0], d = dn[kc:k0-1], T = tn[kc:k0-1]); split = 3.0)
    zone = copy(atm0.zone .+ (kc - 1)); r = copy(atm0.r); jump = copy(atm0.jump); Tp = copy(atm0.T)
    kf = k0 - 1; rb = rn[k0]; Lrun0 = 4π * rb^2 * fr[k0]; nsub = 0
    χf = TM.v(:chifn)[kf] * dn[kf]; τz = χf * (rb - r[end])           # grey optical depth from the front zone's centre to r_b
    if IFSI_SUB && Lrun0 > 0 && τz > 1
        # IONFSSUB (default with IONFSIMP): the photosphere inside zone k₀−1 (an unresolved transparency front, SUBZONE) as sub-points:
        # grey Eddington T⁴ = ¾ T_eff⁴ (τ + 2/3), T_eff from L_run at r_b, capped at the zone's T, on the zone's grey τ scale from r_b.
        # Without it the last stretch to r_b had the zone-centre T (z300mm: 7302 K against T_eff 6000 K): L_fs = 2.2 L_run, and the
        # Wien tail (LW, Lyman continuum) e^{hν/k (1/T_eff − 1/T)} ~ 100× too bright.
        Te_ = (Lrun0 / (4π * rb^2 * 5.670374419e-5))^0.25; rcf = r[end]
        for τ in (300.0, 100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001, 0.0)
            τ < 0.9τz || continue
            push!(r, τ == 0 ? rb * (1 - 1e-9) : rb - (rb - rcf) * τ / τz); push!(zone, kf); push!(jump, false)
            push!(Tp, min(tn[kf], Te_ * (0.75 * (τ + 2 / 3))^0.25)); nsub += 1
        end
    else
        push!(zone, kf); push!(r, rb * (1 - 1e-9)); push!(jump, false); push!(Tp, tn[kf])
    end
    nH = dn[zone] ./ CNn.μH(); ds = [CNn.full(m.x[k], m.net) for k in zone]
    g(s) = [max(Float64(ForwardDiff.value(d[s])), 0.0) for d in ds] .* nH
    op = IFOpac(Tp, g(:H), g(:Hp), g(:e), g(:Hm), g(:H2), g(:He), g(:Hep), g(:Hepp))
    atm = (r = r, ρ = dn[zone], T = Tp, zone = zone, jump = jump)
    νs, dν, grp = IFSI_FREQ
    if nsub > 0       # the sub-points' source opacity scaled to the zone's grey opacity per length (their Rosseland mean = χ ρ of zone k₀−1)
        for j in length(r)-nsub+1:length(r)
            T = Tp[j]; sR = 0.0; wR = 0.0
            for i in eachindex(νs)
                grp[i] >= 0 || continue
                ab, sc, _ = ifs_coeffs3(op, j, νs[i]); x = CNn.hP * νs[i] / (CNn.kB * T)
                dB = x > 600 ? 0.0 : IFCO().Bν(νs[i], T) * x * exp(x) / expm1(x) / T
                sR += dB * dν[i] / (ab + sc); wR += dB * dν[i]
            end
            sR > 0 && (op.sc[j] = χf / (wR / sR))
        end
    end
    _, R = SP_.rays(atm.r; ncore = 8); R = SP_.weights!(atm, R)
    L = zeros(length(νs))
    IFCO().tables()
    Threads.@threads for i in eachindex(νs)
        L[i] = max(SP_.monochromatic(atm, op, R, νs[i]; maxit = 60, maxit_direct = true)[1], 0.0)
    end
    N = zeros(4); Ec = zeros(4); σ = zeros(3, length(IFS_SP))
    Lfs = sum(L[i] * dν[i] for i in eachindex(νs) if grp[i] >= 0)
    for i in eachindex(νs)
        gi_ = grp[i]; gi_ <= 0 && continue
        n_ = L[i] * dν[i] / (CNn.hP * νs[i]); N[gi_] += n_; Ec[gi_] += n_ * CNn.hP * νs[i]
        gi_ <= 3 && for (j, s) in enumerate(IFS_SP); σ[gi_, j] += n_ * ifs_σ(s, νs[i]); end
    end
    for gg in 1:4
        if N[gg] > 0
            Ec[gg] /= N[gg]; gg <= 3 && (σ[gg, :] ./= N[gg])
        else
            e1, e2 = gg <= 3 ? (IFSI_EDGES[gg], IFSI_EDGES[gg+1]) : IFSI_LW
            Ec[gg] = sqrt(e1 * e2) * EVERG; gg <= 3 && for (j, s) in enumerate(IFS_SP); σ[gg, j] = ifs_σ(s, Ec[gg] / CNn.hP); end
        end
    end
    L12 = L[end]
    Lrun = 4π * rn[k0]^2 * fr[k0]
    IFSI_LDIAG[] > 0 && (IFSI_LDIAG[] -= 1; ifsi_ldiag(atm, op, R, νs, dν, grp, L, Lfs, Lrun, zone, k0, tauz))
    f = (IFSNORM && Lfs > 0 && Lrun > 0) ? Lrun / Lfs : 1.0
    st.N = f .* N[1:3]; dE = [diff(IFSI_EDGES)...] .* EVERG
    st.G = Main.IonRay.Groups(Ec[1:3], dE, IFS_SP, σ); st.k0 = k0; st.rb = rn[k0]; st.step = Int(TF.getscalar(:jstep))
    st.fnorm = f; st.Lfs = Lfs; st.Lrun = Lrun
    st.Qsrc = (sum(st.N), st.N[2] + st.N[3], st.N[3])
    I.Phi = f * N[4]; I.σd = N[4] > 0 ? IFSI_K27 * L12 / N[4] : 0.0; I.ELW = Ec[4]
    st.nsrc += 1; st.tsrc += time() - t0
    (st.nsrc <= 20 || st.nsrc % 50 == 0) && @printf("IONFSIMP source step %d: zones %d-%d (r_b %.1f R☉, T %.0f-%.0f K), Q_H %.3e Q_HeI %.3e Q_HeII %.3e Q_LW %.3e s⁻¹ (σ_d %.3e cm²), L_fs %.3e L_run %.3e L☉ (norm %.3g), %.2f s\n",
                                                  st.step, kc - a + 1, k0 - a, rn[k0] / 6.957e10, tn[kc], tn[k0-1], st.Qsrc..., I.Phi, I.σd, Lfs / 3.83e33, Lrun / 3.83e33, f, time() - t0)
    nothing
end

# ---------------------------------------------------------------- per-step set-up (start of the step, COUPLED prepare pass)
"Ionizing Δτ (3 groups) of zone k at absorbers xa = (H, H₂, He, He⁺, …)."
ifsi_tauion(st, xa, nH, dr) = ntuple(g -> (xa[1] * st.G.σ[g, 1] + xa[2] * st.G.σ[g, 2] + xa[3] * st.G.σ[g, 3] + xa[4] * st.G.σ[g, 4]) * nH * dr, 3)
"Start-of-step sweep in ln N of all four groups from the start-of-step absorbers; the ray's zones; the Newton's ray state (TX.YR0, YRN)."
function ifsi_setup!(st, m, a, b, k0, tauz, slrates, dt)
    dn, tn, rn, un = TM.v(:dn), TM.v(:tn), TM.v(:rn), TM.v(:un); mgr = length(dn)
    I = IFSI[]; (I === nothing || length(I.fhi) != mgr) && (I = IFSI[] = ifsi_alloc(mgr))
    k0 = st.k0
    # H I column (Lyman-line shielding of LW, start of step, to the zone centre from k₀) and the Sobolev length factor (SLSHIELD=sobolev)
    s1 = 0.0; sob = STARLIGHT[].shield === :sobolev
    fill!(I.fhi, 1.0); fill!(I.sob, Inf); fill!(I.phz, nothing)
    for k in k0:b
        nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; xH = st.xabs0[k][1]
        STARLIGHT[].shield === :none || (I.fhi[k] = fsh_hi(s1 + 0.5 * nH * xH * dr))
        s1 += nH * xH * dr
        if sob
            dudr = abs(un[k+1] - un[k]) / dr
            dudr > 0 && (I.sob[k] = sqrt(2CNn.kB * tn[k] / 1.6735575e-24) / sqrt(2) / dudr)
        end
    end
    I.ysrc .= log.(max.(vcat(st.N, I.Phi), 1e-300))
    for g in 1:4; length(TX.YR0[][g]) == mgr || (TX.YR0[][g] = zeros(mgr)); length(TX.YRN[][g]) == mgr || (TX.YRN[][g] = zeros(mgr)); end
    y = copy(I.ysrc); last = k0
    for k in 1:mgr; for g in 1:4; TX.YR0[][g][k] = I.ysrc[g]; end; end
    for k in k0:b
        nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3); xa = st.xabs0[k]
        for g in 1:4; TX.YR0[][g][k] = y[g]; end
        τi = ifsi_tauion(st, xa, nH, dr); τL, D = ifsi_lw(I, k, y[4], tn[k], xa[2], nH, dr)
        I.D[k] = D; I.kLW[k] = xa[2] > 0 ? D / (xa[2] * nH * V) : 0.0
        (sum(exp(y[g]) for g in 1:3) * dt / (nH * V) > IFSXMIN || (xa[2] > 1e-20 && I.kLW[k] * dt > IFSI_LWMIN)) && (last = k)
        for g in 1:3; y[g] -= min(τi[g], IFSI_TMAX); end; y[4] -= min(τL, IFSI_TMAX)
    end
    I.kend = min(b, last + 2)
    for g in 1:4; copyto!(TX.YRN[][g], TX.YR0[][g]); end
    TX.YSRC[] = copy(I.ysrc); TX.RAYK[] = (k0, I.kend)
    for k in k0:I.kend
        st.on[k] = true; for g in 1:3; st.nin[k][g] = exp(TX.YR0[][g][k]); end
        slrates === nothing || slrates(k)
        st.p3[k] = (slrates !== nothing && hasproperty(STARLIGHT[], :tl3) && STARLIGHT[].tl3) ? SLDIAG[].tl3[k] : nothing
        ifs_rates!(st, k, st.xabs0[k], tn[k])
        st.r0[k] = (st.fpHe[k].rate[], st.fpHep[k].rate[]); st.κ20[k] = st.κ2[k]
    end
    for k in a:b; st.nin0[k] .= st.nin[k]; end
    st.niter = 0; st.epsit = 0.0; I.niter = 0
    nothing
end
"STARLIGHT's LW rate of zone k replaced by the ray's (zones k ≥ k₀); the SLSHIELD value kept for the comparison."
function ifsi_klw(k, kLW)
    st = IFS[]; I = IFSI[]
    (IFSI_LWON && st !== nothing && I !== nothing && st.k0 > 0 && st.k0 <= k <= length(I.kLW)) || return kLW
    I.kLWsl[k] = kLW
    I.kLW[k]
end
"LOCALJ env of a zone k ≥ k₀ that is not a ray zone: rates below 13.598 eV only, case-B He recombination (as ifs_terms)."
function ifsi_env(envk)
    hasproperty(envk, :ljW) && return merge(envk, (; caseB = true, ljcut = IFS_CUT))
    (hasproperty(envk, :static) && !envk.static && hasproperty(envk, :J)) || return envk
    J0 = envk.J; merge(envk, (; caseB = true, J = ν -> CNn.hP * ν >= IFS_CUT ? 0.0 : J0(ν)))
end
"Does IONFSIMP's ljcut apply to (non-ray) LOCALJ zone k?"
ifsi_ljzone(k) = (st = IFS[]; st !== nothing && (IFSI_LJ == "all" || (st.k0 > 0 && k >= st.k0)))
"Ray zone of the prepare pass: the zone's own terms with the ray's env (case B, ljcut, He backward Euler); no ray rates (they are the Newton's)."
function ifsi_terms(m, k, exs, cools, envk, phk)
    _, _, env2, _ = ifs_terms(m, k, exs, cools, envk, phk)
    I = IFSI[]; I.phz[k] = Tv -> Any[(p, CNn.photorates_t(p, Tv, env2)...) for p in phk if p isa CNn.Photo]
    exs, cools, env2, phk
end
"A split ray zone after the Newton (post pass): IonFS's Float ray terms, with the He recombination cooling (ifsi_herec)."
function ifsi_terms_split(m, k, exs, cools, envk, phk)
    ex2, c2, env2, ph2 = ifs_terms(m, k, exs, cools, envk, phk)
    I = IFSI[]; phz = Tv -> Any[(p, CNn.photorates_t(p, Tv, env2)...) for p in phk if p isa CNn.Photo]; I.phz[k] = phz
    c3 = (ρ, T, d) -> c2(ρ, T, d) + ifsi_herec(ρ, T, d, phz(_fv(T)))
    ex2, c3, env2, ph2
end

# ---------------------------------------------------------------- the coupled zone with the ray
"""
The CK expansion of ray zone k at z = (ρ, T, x_H⁺, x_H₂, y₁…y₄) with the ray's terms, the He photo-rates per atom taken at the absorber
fractions u = (x_H, x_He, x_He⁺): returns (x, f_slow, d).
"""
function ifsi_core(m::NetGas, k, z, t, u)
    # (ck_expand's composition takes the number type of x_s and T: carry the ray's and u's dual numbers there too)
    o_ = zero(u[1]) + zero(u[2]) + zero(u[3]) + zero(z[5]) + zero(z[6]) + zero(z[7]) + zero(z[8]) + zero(z[1])
    st = IFS[]; I = IFSI[]; ρ = z[1]; T = z[2] + o_; rn = TM.v(:rn)
    nH = ρ / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
    ys = (z[5], z[6], z[7], z[8]); N = ntuple(g -> exp(ys[g]), 3); σ = st.G.σ; κ2 = st.κ2[k]
    fr_ = m.net.free; iHp = findfirst(==(:Hp), fr_); iH2 = findfirst(==(:H2), fr_)
    function exray(x)
        d = CNn.full(x, m.net); o = zero(eltype(x)) + zero(N[1]) + zero(ys[4]) + zero(T)
        xH = max(d[:H], zero(d[:H])); xH2 = max(d[:H2], zero(d[:H2])); xHe = max(d[:He], zero(d[:He])); xHep = max(d[:Hep], zero(d[:Hep]))
        xp = max(d[:Hp], zero(d[:Hp]))
        f = fill(o, length(x))
        for g in 1:3
            τ = (xH * σ[g, 1] + xH2 * σ[g, 2] + xHe * σ[g, 3] + xHep * σ[g, 4]) * nH * dr
            _fv(τ) > 0 || continue
            f[iHp] += N[g] * (-expm1(-τ)) * (xH * σ[g, 1] * dr / τ) / V
        end
        f[iHp] -= κ2 * xp * (1 - xp)
        if IFSI_LWON && iH2 !== nothing
            _, D = ifsi_lw(I, k, ys[4], T, xH2, nH, dr); f[iH2] -= D / (nH * V)
        end
        f
    end
    ex = t.extra === nothing ? exray : (x -> t.extra(x) .+ exray(x))
    RHe = zero(N[1]) + zero(T) + zero(u[1]); RHep = zero(RHe)
    for g in 1:3
        τ = (u[1] * σ[g, 1] + u[2] * σ[g, 3] + u[3] * σ[g, 4]) * nH * dr
        fτ = _fv(τ) > 1e-12 ? -expm1(-τ) / τ : one(τ) - τ / 2
        RHe += N[g] * fτ * σ[g, 3] * dr / V; RHep += N[g] * fτ * σ[g, 4] * dr / V
    end
    ph0 = t.phc
    phc = Tv -> vcat(ph0(Tv), Any[(st.fpHe[k], RHe, zero(RHe)), (st.fpHep[k], RHep, zero(RHep))])
    x, fs = CNn.ck_eval([z[3] + o_, z[4] + o_], ρ, T, m.net; env = t.env, phc, extra = ex)
    x, fs, CNn.full(x, m.net)
end
ifsi_u(d) = (max(d[:H], zero(d[:H])), max(d[:He], zero(d[:He])), max(d[:Hep], zero(d[:Hep])))
const IFSI_NNEWT = Threads.Atomic{Int}(0); const IFSI_NFAIL = Threads.Atomic{Int}(0)
const IFSI_HEDBG = Ref(parse(Int, get(ENV, "IONFSHEDBG", "0")))   # (diagnostic) print the first n unconverged He Newtons
"""
He absorbers of ray zone k consistent with their own photo-rates: u* with u = expansion(rates(u)), Newton in ln u on the Float values
(from the last record; converged to 1e-11, so the result does not depend on the start), then one implicit-function correction with the
inputs' dual numbers (w = w* − J⁻¹ F(w*; z)), so that the residual's derivatives are those of the fixed point. Zones whose He photons are
negligible (N_He Δt < 1e-12 of the zone's He) take the start-of-step fractions.
"""
function ifsi_ufix(m::NetGas, k, z, t)
    st = IFS[]; dn, rn = TM.v(:dn), TM.v(:rn); V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
    x0 = st.xabs0[k]; dt = TF.getscalar(:dtime)
    NHe = exp(_fv(z[6])) + exp(_fv(z[7]))
    NHe * dt < 1e-12 * CNn.yHe * dn[k] / CNn.μH() * V && return (x0[1], x0[3], x0[4])
    zf = [_fv(q) for q in z]
    # (fractions floored at 1e-40 as the start: an absorber below it is absent; x_H < 0 from an iterate x_H+ > 1 is 0)
    F(w) = (u = exp.(w); d = ifsi_core(m, k, zf, t, u)[3]; un = ifsi_u(d); [log(max(un[i], 1e-40)) - w[i] for i in 1:3])
    J = zeros(3, 3)
    # Newton with backtracking on |F| (the photon-limited He zone's map is nearly singular: full steps overshoot and cycle); returns
    # (converged, w, |F|) of the best iterate
    function solve(w, maxit)
        bw = copy(w); bn = Inf
        for it in 1:maxit
            Threads.atomic_add!(IFSI_NNEWT, 1)
            res = ForwardDiff.DiffResults.JacobianResult(zeros(3), w); ForwardDiff.jacobian!(res, F, w)
            r = copy(ForwardDiff.DiffResults.value(res)); J .= ForwardDiff.DiffResults.jacobian(res)
            n0 = maximum(abs, r); n0 < bn && (bn = n0; bw .= w)
            # (tight: a looser |F| < 1e-6 left that noise in the residual, and the hydro Newton's median Δt fell from 1.6e3 s to 13 s,
            # t_ifs24 against t_ifs15; where Newton stalls above it — x_H = 1 − x_H⁺ − … at x_H ~ 1e-8, He⁺ under dominant He²⁺ — the
            # best iterate is taken)
            n0 < 1e-11 && return (true, w, n0)
            δ = try -(J \ r) catch; fill(NaN, 3) end; all(isfinite, δ) || return (false, bw, bn)
            δ .= clamp.(δ, -20.0, 20.0)
            maximum(abs, δ) < 1e-11 && (w .+= δ; return (true, w, n0))
            λ = 1.0; wn = w .+ δ
            for _ in 1:12
                maximum(abs, F(wn)) < (1 - 1e-4λ) * n0 && break
                λ /= 2; wn = w .+ λ .* δ
            end
            w = wn
        end
        (false, bw, bn)
    end
    xr = st.xabs[k]
    ok, w, nrm = solve([log(max(xr[1], 1e-40)), log(max(xr[3], 1e-40)), log(max(xr[4], 1e-40))], 60)
    ok || (J .= ForwardDiff.jacobian(F, w))
    if !ok
        Threads.atomic_add!(IFSI_NFAIL, 1)
        if IFSI_HEDBG[] > 0
            IFSI_HEDBG[] -= 1; r_ = F(w)
            @printf("IONFSHEDBG zone %d: He Newton unconverged, |F| %.3e, u %s, z %s, J %s\n", k - Int(TF.getscalar(:ngrs)) + 1, maximum(abs, r_),
                    string(exp.(w)), string(zf), string(J))
        end
    end
    all(z_ -> z_ isa Float64, z) && return Tuple(exp.(w))
    # implicit function: w(z) = w* − J⁻¹ F(w*; z)
    d = ifsi_core(m, k, z, t, Tuple(exp.(w)))[3]; un = ifsi_u(d)
    Fz = [log(max(un[i], 1e-40 + zero(un[i]))) - w[i] for i in 1:3]
    Ji = try inv(J) catch; zeros(3, 3) end; all(isfinite, Ji) || (Ji = zeros(3, 3))
    wz = [w[i] - sum(Ji[i, j] * Fz[j] for j in 1:3) for i in 1:3]
    Tuple(exp.(wz))
end
"""
(p, e, S_H⁺, S_H₂, Λ, cq, yq, Δτ₁, Δτ₂, Δτ₃, Δτ_LW) of coupled ray zone k at z = (ρ, T, x_H⁺, x_H₂[, y₁…y₄]) (y from TX.YRN if absent).
`wantx`: also the CK-expanded composition.
"""
function ifsi_g(m::NetGas, k, z, t; wantx = false)
    st = IFS[]; I = IFSI[]; rn = TM.v(:rn)
    length(z) >= 8 || (z = vcat(z[1:4], [TX.YRN[][g][k] + zero(z[1]) for g in 1:4]))
    ρ = z[1]; T = z[2]; nH = ρ / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
    u = ifsi_ufix(m, k, z, t)
    x, fs, d = ifsi_core(m, k, z, t, u)
    N = ntuple(g -> exp(z[4+g]), 3)
    uE, p = CNn.eos(ρ, T, x, m.net)
    ifs_record!(k, d)
    xc = (max(d[:H], zero(d[:H])), max(d[:H2], zero(d[:H2])), max(d[:He], zero(d[:He])), max(d[:Hep], zero(d[:Hep])))
    τf = ifsi_tauion(st, xc, nH, dr)
    heat = zero(uE)
    for g in 1:3; heat += N[g] * (-expm1(-τf[g])) * st.G.Ec[g]; end
    Λ = (t.cooling === nothing ? zero(uE) : t.cooling(ρ, T, d)) - heat / V + ifsi_herec(ρ, T, d, t.phc(_fv(T)))
    p3_ = hasproperty(t.env, :tl3) ? t.env.tl3 : nothing
    cq, yq = p3_ === nothing ? (zero(uE), zero(uE)) : tl3_heat(ρ, T, d, merge(p3_, (; ne = ρ / CNn.μH() * max(d[TL3NEK], zero(d[TL3NEK])))))
    τL, _ = ifsi_lw(I, k, z[8], T, xc[2], nH, dr)
    out = [p, uE, fs[1], fs[2], Λ, cq, yq, min(τf[1], IFSI_TMAX), min(τf[2], IFSI_TMAX), min(τf[3], IFSI_TMAX), min(τL, IFSI_TMAX)]
    wantx ? (out, x) : out
end
"(diagnostic, with CPLJCHECK) ray zone k's local 8-wide Jacobian against central differences; the He Newton's statistics."
function ifsi_gcheck(m::NetGas, k)
    t = CPLTERMS[][k]; a = Int(TF.getscalar(:ngrs))
    z0 = [TM.v(:dn)[k], TM.v(:tn)[k], TX.XSN[][1][k], TX.XSN[][2][k], (TX.YRN[][g][k] for g in 1:4)...]
    n0, f0 = IFSI_NNEWT[], IFSI_NFAIL[]
    J = ForwardDiff.jacobian(z -> ifsi_g(m, k, z, t), z0)
    @printf("IFSIGCHECK zone %d: He Newton iterations %d, unconverged %d; u = %s\n", k - a + 1, IFSI_NNEWT[] - n0, IFSI_NFAIL[] - f0, string(ifsi_ufix(m, k, z0, t)))
    names = ("p", "e", "S1", "S2", "Λ", "cq", "yq", "τ1", "τ2", "τ3", "τLW"); vn = ("ρ", "T", "x1", "x2", "y1", "y2", "y3", "y4")
    for j in 1:8
        h = j <= 4 ? 1e-6 * abs(z0[j]) : 1e-6; zp = copy(z0); zm = copy(z0); zp[j] += h; zm[j] -= h
        fd = (ifsi_g(m, k, zp, t) .- ifsi_g(m, k, zm, t)) ./ 2h
        bad = [(names[i], J[i, j], fd[i]) for i in 1:11 if abs(J[i, j] - fd[i]) > 1e-3 * max(abs(J[i, j]), abs(fd[i])) && max(abs(J[i, j]), abs(fd[i])) > 0]
        @printf("IFSIGCHECK zone %d d/d%s: %s\n", k - a + 1, vn[j], isempty(bad) ? "all 11 agree to 1e-3" : join((@sprintf("%s AD %.5e FD %.5e", b_...) for b_ in bad), " | "))
    end
    nothing
end
"Δτ of the four groups in a split (not coupled) zone of the ray: start-of-step absorbers and T; LW in the carrier y₄."
function ifsi_tau(k, y1, y2, y3, y4)
    st = IFS[]; I = IFSI[]; dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn)
    nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; xa = st.xabs0[k]
    τi = ifsi_tauion(st, xa, nH, dr); τL, _ = ifsi_lw(I, k, y4, tn[k], xa[2], nH, dr)
    z = zero(y1) + zero(y2) + zero(y3) + zero(y4)
    (min(τi[1], IFSI_TMAX) + z, min(τi[2], IFSI_TMAX) + z, min(τi[3], IFSI_TMAX) + z, min(τL, IFSI_TMAX) + z)
end
@inline ifsi_close(vJ, v) = vJ !== nothing && length(vJ) == 8 && cpl_close(vJ, v) && all(i -> abs(v[i] - vJ[i]) < CPLJLAG[2], 5:8)
function ifsi_vj(m, k, v)
    t0 = time_ns(); z0 = [v...]; t = CPLTERMS[][k]
    res = ForwardDiff.DiffResults.JacobianResult(zeros(11), z0)
    ForwardDiff.jacobian!(res, z -> ifsi_g(m, k, z, t), z0)
    Threads.atomic_add!(CPLPN, 1); Threads.atomic_add!(CPLPT, Int(time_ns() - t0))
    copy(ForwardDiff.DiffResults.value(res)), copy(ForwardDiff.DiffResults.jacobian(res))
end
function ifsi_lift(y, J, zs::NTuple{8,D}) where {D<:ForwardDiff.Dual}
    P = map(ForwardDiff.partials, zs)
    ntuple(i -> D(y[i], J[i, 1] * P[1] + J[i, 2] * P[2] + J[i, 3] * P[3] + J[i, 4] * P[4] + J[i, 5] * P[5] + J[i, 6] * P[6] + J[i, 7] * P[7] + J[i, 8] * P[8]), 11)
end
ifsi_lift(y, J, zs::NTuple{8,Float64}) = Tuple(y)
ifsi_lift(y, J, zs::Tuple) = ifsi_lift(y, J, promote(zs...))
"TitanResidFast.CPL_RAY: cpl_zone for a ray zone (8-wide local Jacobian, cached and reused as CPLJLAG; y within CPLJLAG[2] in ln N)."
function cpl_rayzone(m::NetGas, k, ρ, T, a, b, y1, y2, y3, y4)
    zs = (ρ, T, a, b, y1, y2, y3, y4)
    v = map(q -> Float64(ForwardDiff.value(q)), zs)
    needJ = !all(q -> q isa Float64, zs)
    c = CPLCACHE[][k]
    if c === nothing || c[1] != v
        if needJ && !(c !== nothing && ifsi_close(c[3], v))
            y, J = ifsi_vj(m, k, v); c = (v, y, v, J)
        else
            y = ifsi_g(m, k, [v...], CPLTERMS[][k]); c = c === nothing ? (v, y, nothing, nothing) : (v, y, c[3], c[4])
            Threads.atomic_add!(CPLPH, 1)
        end
        CPLCACHE[][k] = c
    elseif needJ && !ifsi_close(c[3], v)
        y, J = ifsi_vj(m, k, v); c = (v, y, v, J); CPLCACHE[][k] = c
    else
        Threads.atomic_add!(CPLPH, 1)
    end
    needJ ? ifsi_lift(c[2], c[4], zs) : Tuple(c[2])
end
"cpl_finish! of a ray zone: the CK expansion with the converged ray."
ifsi_finish_x(m, k) = ifsi_g(m, k, [TM.v(:dn)[k], TM.v(:tn)[k], TX.XSN[][1][k], TX.XSN[][2][k]], CPLTERMS[][k]; wantx = true)[2]

# ---------------------------------------------------------------- after the Newton (COUPLED post pass, before the zones)
"The converged ray into IonFS's state (photons per zone, He rates of split zones), the LW rates of every zone k ≥ k₀, photon bookkeeping."
function ifsi_post!(m, a, b)
    st = IFS[]; I = IFSI[]; (st === nothing || st.k0 == 0 || I === nothing) && return nothing
    dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn); k0, ke = TX.RAYK[]
    for k in a:b; fill!(st.nin[k], 0.0); end
    fill!(I.absb, 0.0); I.emit .= exp.(I.ysrc); rowe = 0.0
    fill!(I.yint, NaN); fill!(I.qion, 0.0)
    y = copy(I.ysrc)
    for k in k0:b
        k <= ke && (for g in 1:4; y[g] = TX.YRN[][g][k]; end)
        for g in 1:4; I.yint[g, k] = y[g]; end
        nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
        TX.cplz(k) && (st.xabs[k] = ifs_x(CNn.full(m.x[k], m.net)))
        xa = TX.cplz(k) ? st.xabs[k] : st.xabs0[k]
        for g in 1:3; st.nin[k][g] = exp(y[g]); end
        τi = ifsi_tauion(st, xa, nH, dr); τL, D = ifsi_lw(I, k, y[4], tn[k], xa[2], nH, dr)
        I.D[k] = D; I.kLW[k] = xa[2] > 0 ? D / (xa[2] * nH * V) : 0.0
        τs = (τi..., τL)
        for g in 1:4; I.absb[g] += exp(y[g]) * (-expm1(-τs[g])); end
        for g in 1:3; I.qion[k] += exp(y[g]) * (-expm1(-τs[g])) * st.G.Ec[g] / V; end
        yn = [y[g] - min(τs[g], IFSI_TMAX) for g in 1:4]
        if k + 1 <= ke
            for g in 1:4; rowe = max(rowe, abs(exp(TX.YRN[][g][k+1]) - exp(yn[g])) / I.emit[g]); end
        end
        y = yn
    end
    I.esc .= exp.(y); I.rowerr = rowe
    for g in 1:4; I.yint[g, b+1] = y[g]; end
    for k in k0:ke
        st.on[k] && ifs_rates!(st, k, TX.cplz(k) ? st.xabs[k] : st.xabs0[k], tn[k])
    end
    st.eps = rowe
    nothing
end
"One line per step in <RUNREC>/ionfsi.csv."
function ifsi_log(m, a, b, js)
    st = IFS[]; I = IFSI[]; (st === nothing || st.k0 == 0 || I === nothing) && return nothing
    dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn)
    dir = get(ENV, "RUNREC", get(ENV, "POPIII_OUT", ".")); f = joinpath(dir, "ionfsi.csv"); new = !isfile(f)
    # LW: dissociations of the ray and of STARLIGHT's SLSHIELD path at the same composition, zones k ≥ k₀
    Dr = 0.0; Ds = 0.0
    for k in st.k0:b
        nH = dn[k] / CNn.μH(); V = 4π / 3 * (rn[k+1]^3 - rn[k]^3); xH2 = max(CNn.full(m.x[k], m.net)[:H2], 0.0)
        Dr += I.kLW[k] * xH2 * nH * V; Ds += I.kLWsl[k] * xH2 * nH * V
    end
    # He²⁺ (per He) inside 3e3 R☉, network zones outside the grey photosphere, and every zone
    h2n = 0.0; h2a = 0.0; kh2 = a - 1
    for k in a:b
        0.5 * (rn[k] + rn[k+1]) < 3e3 * 6.957e10 || continue
        x = max(CNn.full(m.x[k], m.net)[:Hepp], 0.0) / CNn.yHe
        h2a = max(h2a, x); (k >= st.k0 && !m.eqzone[k]) && x > h2n && (h2n = x; kh2 = k)
    end
    eb(g) = abs(I.emit[g] - I.absb[g] - I.esc[g]) / max(I.emit[g], 1e-300)
    open(f, "a") do io
        new && println(io, "step,t_yr,dt_s,k0,kend,r_b_Rsun,Q_H,Q_HeI,Q_HeII,Q_LW,abs_H,abs_HeI,abs_HeII,abs_LW,esc_H,esc_HeI,esc_HeII,esc_LW,err_H,err_HeI,err_HeII,err_LW,rowerr,D_ray,D_slshield,sigma_d,xHepp_net_3e3,xHepp_all_3e3,fnorm,L_fs_Lsun,L_run_Lsun,ncoupled_ray,nnewt_He,nfail_He,Tmax_ray,zTmax,khepp")
        nc = count(k -> TX.cplz(k), st.k0:I.kend)
        Tm = 0.0; kTm = 0; for k in st.k0:I.kend; tn[k] > Tm && (Tm = tn[k]; kTm = k); end
        @printf(io, "%d,%.8e,%.4e,%d,%d,%.5g,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.3e,%.3e,%.3e,%.3e,%.3e,%.5e,%.5e,%.4e,%.4e,%.4e,%.5g,%.5e,%.5e,%d,%d,%d,%.1f,%d,%d\n",
                js, TF.getscalar(:timen) / 3.156e7, TF.getscalar(:timen) - TF.getscalar(:timeo), st.k0 - a + 1, I.kend - a + 1, st.rb / 6.957e10,
                I.emit..., I.absb..., I.esc..., eb(1), eb(2), eb(3), eb(4), I.rowerr, Dr, Ds, I.σd, h2n, h2a, st.fnorm, st.Lfs / 3.83e33, st.Lrun / 3.83e33, nc,
                IFSI_NNEWT[], IFSI_NFAIL[], Tm, kTm - a + 1, kh2 - a + 1)
    end
    nothing
end

# ---------------------------------------------------------------- per-snapshot diagnostics (ionfs_snap.csv)
const IFSI_RADII = (1e2, 1e3, 1e4, 1e5)                                   # R☉
"""
<RUNREC>/ionfs_snap.csv, one line per snapshot (IONFSIMP): does the light get ahead of the photosphere, and what does it do there.
Columns (radii in R☉, rates in s⁻¹, heating per mass in erg g⁻¹ s⁻¹; "ph" = the ray's start r_b, the first interface outside the
grey photosphere, τ_grey = 2/3 from outside):
- step, t_yr, dt_s, r_ph_Rsun, T_eff (STARLIGHT's), T_ph (gas T of the zone just inside r_ph), L_run_Lsun, L_fs_Lsun, fnorm
- Qem_LW, Qem_H, Qem_HeI, Qem_HeII: photons per s of the four groups leaving r_ph (the formal solution's emergent spectrum, × fnorm)
- Qesc_*: the same through the outer boundary
- rIF_H, rIF_HeI, rIF_HeII: first interface outside r_ph beyond which x_H⁺ < 0.5, (x_He⁺ + x_He²⁺)/y_He < 0.5, x_He²⁺/y_He < 0.5
  (= r_ph: no ionized region ahead of the photosphere); rH2_01, rH2_05: first zone centre outside r_ph with f_H₂ = 2x_H₂ ≥ 0.1, 0.5
  (NaN: none)
- Dray, Dsl: H₂ dissociations per s outside r_ph, from the ray and from STARLIGHT's SLSHIELD rate at the same composition
- Qion_heat, LW_heat: power of the ray's ionizing photo-heating and of the LW heat (4.88 eV per dissociation) outside r_ph [erg/s]
- per radius R (zone containing R / interface at or above R): Q{LW,H,HeI,HeII}_R (photons crossing it), T_R, u_R (km/s), xe_R, fH2_R,
  kLWray_R, kLWsl_R (per molecule), qray_R (ionizing + LW heat per mass), tflow_theat_R = (R/|u|)/(e/qray): ≳ 1 means the light heats
  the infall noticeably before it arrives
"""
function ifsi_snaplog(m, a, b, js)
    st = IFS[]; I = IFSI[]; (st === nothing || st.k0 == 0 || I === nothing) && return nothing
    dn, tn, rn, un, egn = TM.v(:dn), TM.v(:tn), TM.v(:rn), TM.v(:un), TM.v(:egn); k0 = st.k0; RS = 6.957e10
    fl(k) = CNn.full(m.x[k], m.net)
    xv(d, s) = max(Float64(ForwardDiff.value(d[s])), 0.0)
    front(f) = (for k in k0:b; f(fl(k)) || return rn[k] / RS; end; rn[b+1] / RS)
    rIH = front(d -> xv(d, :Hp) >= 0.5)
    rIHe1 = front(d -> (xv(d, :Hep) + xv(d, :Hepp)) / CNn.yHe >= 0.5)
    rIHe2 = front(d -> xv(d, :Hepp) / CNn.yHe >= 0.5)
    h2f(th) = (for k in k0:b; 2xv(fl(k), :H2) >= th && return 0.5 * (rn[k] + rn[k+1]) / RS; end; NaN)
    Dr = 0.0; Ds = 0.0; Pion = 0.0; Plw = 0.0
    for k in k0:b
        nH = dn[k] / CNn.μH(); V = 4π / 3 * (rn[k+1]^3 - rn[k]^3); xH2 = xv(fl(k), :H2)
        Dr += I.kLW[k] * xH2 * nH * V; Ds += I.kLWsl[k] * xH2 * nH * V; Pion += I.qion[k] * V
        Plw += I.kLW[k] * xH2 * nH * V * (4.48 + 0.4) * EVERG
    end
    sd = SLDIAG[]; Te = sd === nothing ? NaN : Float64(sd.Te)
    Q(g, k) = (y = I.yint[g, k]; isnan(y) ? 0.0 : exp(y))
    cols = Any[js, TF.getscalar(:timen) / 3.156e7, TF.getscalar(:dtime), rn[k0] / RS, Te, tn[k0-1], st.Lrun / 3.83e33, st.Lfs / 3.83e33, st.fnorm,
               (exp(I.ysrc[g]) for g in (4, 1, 2, 3))..., (Q(g, b + 1) for g in (4, 1, 2, 3))..., rIH, rIHe1, rIHe2, h2f(0.1), h2f(0.5),
               Dr, Ds, Pion, Plw]
    for R in IFSI_RADII
        ki = findfirst(r -> r >= R * RS, rn[k0:b+1]); ki = ki === nothing ? b + 1 : ki + k0 - 1      # interface at or above R
        kz = clamp(ki - 1, k0, b); d = fl(kz); nH = dn[kz] / CNn.μH(); xH2 = xv(d, :H2)
        q = I.qion[kz] / dn[kz] + I.kLW[kz] * xH2 * nH * (4.48 + 0.4) * EVERG / dn[kz]
        uR = 0.5 * (un[kz] + un[kz+1]); tfl = 0.5 * (rn[kz] + rn[kz+1]) / max(abs(uR), 1e-30); the = egn[kz] / max(q, 1e-300)
        append!(cols, [(Q(g, ki) for g in (4, 1, 2, 3))..., tn[kz], uR / 1e5, xv(d, :e), 2xH2, I.kLW[kz], I.kLWsl[kz], q, tfl / the])
    end
    dir = get(ENV, "RUNREC", get(ENV, "POPIII_OUT", ".")); f = joinpath(dir, "ionfs_snap.csv"); new = !isfile(f)
    open(f, "a") do io
        if new
            h = ["step", "t_yr", "dt_s", "r_ph_Rsun", "T_eff", "T_ph", "L_run_Lsun", "L_fs_Lsun", "fnorm",
                 "Qem_LW", "Qem_H", "Qem_HeI", "Qem_HeII", "Qesc_LW", "Qesc_H", "Qesc_HeI", "Qesc_HeII",
                 "rIF_H_Rsun", "rIF_HeI_Rsun", "rIF_HeII_Rsun", "rH2_01_Rsun", "rH2_05_Rsun", "Dray", "Dsl", "Qion_heat", "LW_heat"]
            for R in IFSI_RADII
                t = string(Int(R)); append!(h, ["QLW_$t", "QH_$t", "QHeI_$t", "QHeII_$t", "T_$t", "u_kms_$t", "xe_$t", "fH2_$t",
                                                "kLWray_$t", "kLWsl_$t", "qray_$t", "tflow_theat_$t"])
            end
            println(io, join(h, ","))
        end
        println(io, join((c isa Integer ? string(c) : @sprintf("%.6e", c) for c in cols), ","))
    end
    nothing
end
