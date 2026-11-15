# Included in PopIIIMaterial. IONFS=1: the Lyman continuum as a photon-conserving ray sourced by a formal solution.
"""
IONFS=1 (default off; nothing below runs, and no code path changes, when it is off).

The hydrogen- and helium-ionizing continuum (hν ≥ 13.598 eV) as a radial, multi-group, photon-conserving ray (Abel, Norman & Madau
1999, ApJ 523, 66; Mellema et al. 2006, NewA 11, 374), case B, through the network zones. The source is not a T_eff blackbody but
the continuum leaving the run's table (LTE) region. The source is set as follows.

1. Source (`ifs_source!`, every IONFSN steps, or when the first network zone k0 moves). The static spherical transfer equation
   (Spectrum.monochromatic: tangent rays, Feautrier, ALI with coherent Thomson scattering, thermal inner boundary) is solved on the
   zones between the innermost zone deeper than grey τ = 300 and k0, ending at r_b = r(k0). Opacity (`IFOpac`): SnapSpec's
   continuum (ContinuumOpacity.alpha: CIA, H⁻ ff, Rayleigh, Thomson; H⁻ bound–free with attachment emission; H bound–free n ≤ 30
   with Milne emission of the zone's own n_e n_p, excited levels min(Saha, Boltzmann)). On top come H⁺, He⁺, He²⁺ free–free (Gaunt 1)
   and He I 1¹S and He II 1s bound–free with Milne emission of the zone's own n_e n_He⁺ and n_e n_He²⁺. The cross sections of
   H I 1s, He I 1s² and He II 1s are those of the network, Verner et al. (1996, ApJ 465, 487; ChemRates.σ_verner). There are
   11 groups (edges below), 3 frequencies each. N_g = ∫_g L_ν/(hν) dν at r_b. The group cross sections are photon-weighted over
   the solution's own spectrum, and E_g is the mean photon energy. The run's T structure comes from grey diffusion, so ∫L_ν dν
   need not equal the run's L at r_b. With IONFSNORM=1 (default) N_g is scaled by L_run(r_b)/L_FS(r_b): the spectral shape of the
   formal solution, the luminosity of the run. The ratio is reported (SpecFull's T⁴ factor fT plays the same role there).
   Outside the solved region (r > r_b) there is no incident light, and the outer network gas's own diffuse field is left to case B.
2. Ray (`ifs_sweep!`): zone by zone outward from k0, N_g,k+1 = N_g,k e^{−Δτ_g}, Δτ_g = Σ_s x_s n_H σ_g,s Δr,
   s ∈ {H, H₂, He, He⁺}. Each zone's absorption uses its new abundances: Nin_g (1 − e^{−Δτ_g(x)}) split over the absorbers by
   Δτ_g,s/Δτ_g. A zone can never take more photons than enter it, however long the step or wherever the front sits inside it.
   H and H₂: IonRay.source's `extra` (H → H⁺ + e, H₂ → H₂⁺ + e), implicit in the zone's own x. He I, He II: per-atom rates
   R = Σ_g Nin_g ((1 − e^{−Δτ_g})/Δτ_g) σ_g Δr / V (ChemNet.FixedPhoto), because the CK network holds He in ionization equilibrium,
   which `extra` does not reach. They are photon-conserving at the fixed point (R x_He n_H V = the He share of the absorption).
3. Energy: heat = Σ_g Nin_g (1 − e^{−Δτ_g(x)}) E_g / V, implicit in the zone's chemistry–temperature step (the COUPLED energy row,
   or ChemNet.step). ChemNet pays χ from u through the species' E0, so the gas keeps hν − χ. The photon energy leaves the grey field
   through RAYQ (as IonRay, and as the three-level heat).
4. No double counting. The ray carries only the light that left the LTE region (the net outward flux at r_b). In the ray zones the
   diffuse Lyman-continuum field of the gas's own recombinations is not transported: recombinations straight to the ground state are
   left out (case B, on the spot). For H the three-level atom does this (C α_B). For He ChemNet's `caseB` env does it, in the LOCALJ
   tables via `ljcut` + `caseB`. LOCALJ's diluted blackbody keeps its photons below 13.598 eV (`ljcut`): the ionizing part of
   W B(T_c), applied locally without attenuation, would count the ray's photons again. So the network sees J − S of the Lyman
   continuum: the ray's photons in, its on-the-spot recombinations out.
5. Front inside a zone. Photon conservation (2.) fixes how many atoms a zone ionizes. The recombinations of a zone holding a sharp
   front (the neutral mean free path ≪ Δr) are those of the ionized part alone: α n_H² V_ion = α n_H² x_p V, not the zone mean
   α n_H² x_p x_e V. The extra recombination κ₂ x_p (1 − x_p) with κ₂ = w C α_B(T) n_H, w = w_ray (1 − e^{−τ̄_H/IONFSTP}), applies
   where the ray dominates the ionization (w_ray = Γ_ray/(Γ_ray + I_3lev + k1 n_e)) and the neutral column is thick (τ̄_H, photon
   weighted). The front radius is reported from the ionized volume fraction, r_f³ = r₁³ + x_p (r₂³ − r₁³).
6. Coupling with the COUPLED Newton (x_H⁺ a Newton unknown). The zones' terms are prepared from the start-of-step state. After every
   Newton update (TitanCore CPLH.update) the ray is swept again from the iterate's abundances (recorded where the residual evaluates
   each zone, `ifs_record!`), and the zones whose Nin, He rates or κ₂ changed lose their cached Jacobian. The converged step is then
   self-consistent: the Jacobian lacks only the inter-zone coupling (a quasi-Newton). Split (not coupled) zones use the sweep of the
   start of the step (lagged by one step). After the step the sweep is repeated from the final state. ε = Σ_k |A_k(final) − A_k(used)|
   / Q_src (A the photons absorbed per zone) is the convergence check. With ε > IONFSTOL the next Δt is halved (DT hook).
Table zones beyond k0 that the ray would change (Nin Δt/(n_H V) > IONFSXMIN) go to the network, as in IonRay.

Diagnostics: <RUNREC>/ionfs.csv per step, and the snapshot field `ionfs` (Q_H, Q_HeI, Q_HeII at the interfaces, k0, front).
"""
const IONFS = get(ENV, "IONFS", "0") == "1"
const IFSN = parse(Int, get(ENV, "IONFSN", "20"))
const IFSXMIN = parse(Float64, get(ENV, "IONFSXMIN", "1e-6"))
const IFSTP = parse(Float64, get(ENV, "IONFSTP", "10"))
const IFSNORM = get(ENV, "IONFSNORM", "1") == "1"
const IFSTOL = parse(Float64, get(ENV, "IONFSTOL", "0.1"))
# (tests) IONFSITER=0: no re-sweep after the Newton updates (start-of-step ray only); IONFSENV=0: no case-B/ljcut env change in the ray zones;
# IONFSDBG=n: the ray's zone table at the first n set-ups
const IFSITER = get(ENV, "IONFSITER", "1") == "1"; const IFSENV = get(ENV, "IONFSENV", "1") == "1"; const IFSDBG = Ref(parse(Int, get(ENV, "IONFSDBG", "0")))
const IFS_EDGES = [13.598, 15.0, 17.0, 20.0, 24.587, 30.0, 40.0, 54.418, 70.0, 100.0, 150.0, 300.0]
const IFS_SP = [:H, :H2, :He, :Hep]
const IFS_CUT = 13.598 * EVERG
# (ContinuumOpacity, ChemRates, Spectrum and IonRay are looked up in Main at run time: other scripts load PopIIIMaterial without them)
IFCO() = Main.ContinuumOpacity; IFSP() = Main.Spectrum
function ifs_σ(s, ν)
    CR_ = Main.ChemRates
    s === :H ? CR_.σ_verner(ν, 13.598, CR_.VERNER_HI) : s === :He ? CR_.σ_verner(ν, 24.587, CR_.VERNER_HeI) :
    s === :Hep ? CR_.σ_verner(ν, 54.418, CR_.VERNER_HeII) : 0.0      # H₂: trace (≤ 1e-15) where the ray reaches; its photoionization (15.4 eV,
    # 1e7 s⁻¹ per molecule at r_b) made the coupled H₂ row stiff to divergence (t_ion1 B1, first step) and is left out
end
mutable struct IFState
    G::Any; N::Vector{Float64}; k0::Int; rb::Float64; step::Int; fnorm::Float64; Lfs::Float64; Lrun::Float64; Qsrc::NTuple{3,Float64}
    nin::Vector{Vector{Float64}}; nin0::Vector{Vector{Float64}}; on::Vector{Bool}
    xabs::Vector{NTuple{6,Float64}}; xabs0::Vector{NTuple{6,Float64}}
    fpHe::Vector{Any}; fpHep::Vector{Any}; κ2::Vector{Float64}; κ20::Vector{Float64}; r0::Vector{NTuple{2,Float64}}
    absd::Vector{Float64}; p3::Vector{Any}; eps::Float64; epsit::Float64; niter::Int; nsrc::Int; tsrc::Float64
end
const IFS = Ref{Any}(nothing)
function ifs_alloc(mgr, ng)
    IFState(nothing, zeros(ng), 0, 0.0, -10^9, 1.0, 0.0, 0.0, (0.0, 0.0, 0.0), [zeros(ng) for _ in 1:mgr], [zeros(ng) for _ in 1:mgr], falses(mgr) |> collect,
            [ntuple(_ -> 0.0, 6) for _ in 1:mgr], [ntuple(_ -> 0.0, 6) for _ in 1:mgr],
            Any[CNn.FixedPhoto("He + γ → He+ + e (IONFS ray)", :He, [:Hep, :e], Ref(0.0)) for _ in 1:mgr],
            Any[CNn.FixedPhoto("He+ + γ → He++ + e (IONFS ray)", :Hep, [:Hepp, :e], Ref(0.0)) for _ in 1:mgr],
            zeros(mgr), zeros(mgr), [(0.0, 0.0) for _ in 1:mgr], zeros(mgr), Any[nothing for _ in 1:mgr], 0.0, 0.0, 0, 0, 0.0)
end
ifs_x(d) = (Float64(ForwardDiff.value(max(d[:H], zero(d[:H])))), Float64(ForwardDiff.value(max(d[:H2], zero(d[:H2])))),
            Float64(ForwardDiff.value(max(d[:He], zero(d[:He])))), Float64(ForwardDiff.value(max(d[:Hep], zero(d[:Hep])))),
            Float64(ForwardDiff.value(max(d[:Hp], zero(d[:Hp])))), Float64(ForwardDiff.value(max(d[:e], zero(d[:e])))))
"The residual's evaluation of coupled zone k at the iterate (cpl_g): its absorbers, for the next sweep."
ifs_record!(k, d) = (st = IFS[]; (st !== nothing && k <= length(st.xabs)) && (st.xabs[k] = ifs_x(d)); nothing)

# ---------------------------------------------------------------- the source: formal solution of the table region
struct IFOpac
    T::Vector{Float64}; nH0::Vector{Float64}; nHp::Vector{Float64}; ne::Vector{Float64}; nHm::Vector{Float64}; nH2::Vector{Float64}
    nHe0::Vector{Float64}; nHe1::Vector{Float64}; nHe2::Vector{Float64}
end
const IF_E0HM = 0.754 * EVERG; const IF_χH = 13.598 * EVERG; const IF_σ0 = 7.906e-18
const IF_χHe1 = 24.587387 * EVERG; const IF_χHe2 = 54.417763 * EVERG
"(absorption, scattering, thermal source) at point k: SnapSpec's continuum point by point, plus He I/He II ground bound–free and He ff."
function ifs_coeffs3(o::IFOpac, k, ν)
    IFCO_ = IFCO(); T = o.T[k]; kT = CNn.kB * T; hν = CNn.hP * ν; c2 = 2CNn.hP / CNn.c^2
    a, s, p = IFCO_.alpha(ν, T; nH = o.nH0[k], nH2 = o.nH2[k], nHe = o.nHe0[k], nHp = o.nHp[k], ne = o.ne[k], nHm = o.nHm[k])
    B = IFCO_.Bν(ν, T); stim = -expm1(-hν / kT)
    ab = a - p.hm_bf - p.h + 3.692e8 * o.ne[k] * (o.nHp[k] + o.nHe1[k] + 4o.nHe2[k]) / (sqrt(T) * ν^3) * stim
    em = ab * B
    λe3 = (CNn.hP^2 / (2π * CNn.me * kT))^1.5
    if hν > IF_E0HM                                                       # H⁻ bound–free, attachment emission
        σ = IFCO_.σ_Hm_bf(1e4 * CNn.c / ν)
        ab += o.nHm[k] * σ * stim; em += o.nH0[k] * o.ne[k] * λe3 / 4 * c2 * σ * ν^3 * exp(-(hν - IF_E0HM) / kT)
    end
    nep = o.ne[k] * o.nHp[k]; lnep = nep > 0 ? log(nep * λe3) : -Inf; ln1 = o.nH0[k] > 0 ? log(o.nH0[k]) : -Inf
    for n in 1:30                                                          # H bound–free (1s: Verner, n ≥ 2 Kramers)
        χn = IF_χH / n^2; hν >= χn || continue
        σ = n == 1 ? ifs_σ(:H, ν) : IF_σ0 * n * (χn / hν)^3
        ls = log(n^2) + lnep + (χn - hν) / kT; nst = ls > -700 ? exp(ls) : 0.0                  # n_n* e^{−hν/kT}
        nn = n == 1 ? o.nH0[k] : n^2 * exp(min(lnep + χn / kT, ln1 - (IF_χH - χn) / kT))
        ab += σ * max(nn - nst, 0.0); em += σ * nst * c2 * ν^3
    end
    for (n0, ni, χ, gr, σs) in ((o.nHe0[k], o.nHe1[k], IF_χHe1, 1 / 4, :He), (o.nHe1[k], o.nHe2[k], IF_χHe2, 1.0, :Hep))   # g/(2g₊): 1/(2·2), 2/(2·1)
        hν >= χ || continue
        σ = ifs_σ(σs, ν); ls = o.ne[k] * ni > 0 ? log(o.ne[k] * ni * gr * λe3) + (χ - hν) / kT : -Inf; nst = ls > -700 ? exp(ls) : 0.0
        ab += σ * max(n0 - nst, 0.0); em += σ * nst * c2 * ν^3
    end
    ab = max(ab, 1e-300)
    (ab, s, em / ab)
end
"Group sub-frequencies (3 per group, log mid-points) and their widths; below 13.598 eV 24 for the luminosity normalisation."
function ifs_freqs()
    νs = Float64[]; dν = Float64[]; grp = Int[]
    Elo = exp.(range(log(0.05), log(13.598), length = 25))
    for i in 1:24; push!(νs, sqrt(Elo[i] * Elo[i+1]) * EVERG / CNn.hP); push!(dν, (Elo[i+1] - Elo[i]) * EVERG / CNn.hP); push!(grp, 0); end
    for g in 1:length(IFS_EDGES)-1
        e1, e2 = IFS_EDGES[g], IFS_EDGES[g+1]; q = (e2 / e1)^(1 / 3)
        for i in 1:3
            push!(νs, e1 * q^(i - 0.5) * EVERG / CNn.hP); push!(dν, e1 * (q^i - q^(i - 1)) * EVERG / CNn.hP); push!(grp, g)
        end
    end
    νs, dν, grp
end
const IF_FREQ = ifs_freqs()
"IONFS start (popiii_m2c, after Spectrum.jl is loaded): Spectrum's coefficient method for the source opacity."
ifs_init!() = (@eval Main.Spectrum.coeffs3(o::$IFOpac, k, ν) = $ifs_coeffs3(o, k, ν); nothing)
function ifs_source!(st, m, a, b, k0, tauz)
    t0 = time(); dn, tn, rn, fr = TM.v(:dn), TM.v(:tn), TM.v(:rn), TM.v(:frn)
    kc = a; for k in a:k0-1; tauz[k] > 300 && (kc = k); end
    kc = max(a, min(kc, k0 - 3))
    SP_ = IFSP(); atm0 = SP_.atmosphere((r = rn[kc:k0], d = dn[kc:k0-1], T = tn[kc:k0-1]); split = 3.0)
    zone = vcat(atm0.zone .+ (kc - 1), k0 - 1); r = vcat(atm0.r, rn[k0] * (1 - 1e-9)); jump = vcat(atm0.jump, false)
    nH = dn[zone] ./ CNn.μH(); ds = [CNn.full(m.x[k], m.net) for k in zone]
    g(s) = [max(Float64(ForwardDiff.value(d[s])), 0.0) for d in ds] .* nH
    op = IFOpac(tn[zone], g(:H), g(:Hp), g(:e), g(:Hm), g(:H2), g(:He), g(:Hep), g(:Hepp))
    atm = (r = r, ρ = dn[zone], T = tn[zone], zone = zone, jump = jump)
    _, R = SP_.rays(atm.r; ncore = 8); R = SP_.weights!(atm, R)
    νs, dν, grp = IF_FREQ; L = zeros(length(νs))
    IFCO().tables()
    Threads.@threads for i in eachindex(νs)
        L[i] = max(SP_.monochromatic(atm, op, R, νs[i]; maxit = 60, maxit_direct = true)[1], 0.0)
    end
    ng = length(IFS_EDGES) - 1; N = zeros(ng); Ec = zeros(ng); σ = zeros(ng, length(IFS_SP)); dE = diff(IFS_EDGES) .* EVERG
    Lfs = sum(L .* dν)
    for i in eachindex(νs)
        gi_ = grp[i]; gi_ == 0 && continue
        n_ = L[i] * dν[i] / (CNn.hP * νs[i]); N[gi_] += n_; Ec[gi_] += n_ * CNn.hP * νs[i]
        for (j, s) in enumerate(IFS_SP); σ[gi_, j] += n_ * ifs_σ(s, νs[i]); end
    end
    for gg in 1:ng
        if N[gg] > 0
            Ec[gg] /= N[gg]; σ[gg, :] ./= N[gg]
        else                                                   # no photons: the group's centre
            Ec[gg] = sqrt(IFS_EDGES[gg] * IFS_EDGES[gg+1]) * EVERG; for (j, s) in enumerate(IFS_SP); σ[gg, j] = ifs_σ(s, Ec[gg] / CNn.hP); end
        end
    end
    Lrun = 4π * rn[k0]^2 * fr[k0]
    f = (IFSNORM && Lfs > 0 && Lrun > 0) ? Lrun / Lfs : 1.0
    st.N .= f .* N; st.G = Main.IonRay.Groups(Ec, dE, IFS_SP, σ); st.k0 = k0; st.rb = rn[k0]; st.step = Int(TF.getscalar(:jstep))
    st.fnorm = f; st.Lfs = Lfs; st.Lrun = Lrun
    st.Qsrc = (sum(st.N), sum(st.N[gg] for gg in 1:ng if IFS_EDGES[gg] >= 24.587), sum(st.N[gg] for gg in 1:ng if IFS_EDGES[gg] >= 54.418))
    st.nsrc += 1; st.tsrc += time() - t0
    (st.nsrc <= 20 || st.nsrc % 50 == 0) && @printf("IONFS source step %d: zones %d-%d (r_b %.1f R☉, T %.0f-%.0f K), Q_H %.3e Q_HeI %.3e Q_HeII %.3e s⁻¹, L_fs %.3e L_run %.3e L☉ (norm %.3g), %.2f s\n",
                                                  st.step, kc - a + 1, k0 - a, rn[k0] / 6.957e10, tn[kc], tn[k0-1], st.Qsrc..., Lfs / 3.83e33, Lrun / 3.83e33, f, time() - t0)
    nothing
end

# ---------------------------------------------------------------- the ray
"Photons per group entering every zone (zones below k0: none), from the absorbers xs(k)."
function ifs_sweep!(st, xs, a, b)
    dn, rn = TM.v(:dn), TM.v(:rn); N = copy(st.N); G = st.G
    for k in a:b
        if k < st.k0; fill!(st.nin[k], 0.0); continue; end
        st.nin[k] .= N; x = xs(k); nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]
        for gg in eachindex(N)
            τ = 0.0; for j in 1:4; τ += x[j] * nH * G.σ[gg, j] * dr; end
            N[gg] *= exp(-τ)
        end
    end
    N
end
"He rates, two-phase recombination κ₂, and the photons absorbed in zone k at the absorbers x and gas T."
function ifs_rates!(st, k, x, T)
    dn, rn = TM.v(:dn), TM.v(:rn); nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
    G = st.G; Nin = st.nin[k]; RHe = 0.0; RHep = 0.0; AH = 0.0; A = 0.0; τn = 0.0; τd = 0.0
    for gg in eachindex(Nin)
        Nin[gg] > 0 || continue
        τ = 0.0; for j in 1:4; τ += x[j] * nH * G.σ[gg, j] * dr; end
        fτ = τ > 1e-12 ? -expm1(-τ) / τ : 1.0
        RHe += Nin[gg] * fτ * G.σ[gg, 3] * dr / V; RHep += Nin[gg] * fτ * G.σ[gg, 4] * dr / V
        τH = x[1] * nH * G.σ[gg, 1] * dr; AH += Nin[gg] * fτ * τH; A += Nin[gg] * fτ * τ
        τn += Nin[gg] * τH; τd += Nin[gg]
    end
    st.fpHe[k].rate[] = RHe; st.fpHep[k].rate[] = RHep; st.absd[k] = A
    ne = nH * x[6]; p3 = st.p3[k]
    ΓR = AH / (nH * max(x[1], 1e-30) * V)
    Io = PC.k1_aazn(T) * ne
    C = 1.0
    if p3 !== nothing
        C = tl3_C(T, merge(p3, (; ne))); Io += (1 - C) * (ne * q12_H(T) + p3.P12)
    end
    w = (ΓR / (ΓR + Io + 1e-300)) * (τd > 0 ? -expm1(-(τn / τd) / IFSTP) : 0.0)
    st.κ2[k] = w * C * alphaB_recfast(T) * nH
    nothing
end
"Per-step set-up (start of the step: the prepare pass, or the only pass without COUPLED)."
function ifs_setup!(m, a, b, tauz, slrates, dt)
    dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn); h = HOT[]
    st = IFS[]; mgr = length(dn); ng = length(IFS_EDGES) - 1
    (st === nothing || length(st.on) != mgr) && (st = IFS[] = ifs_alloc(mgr, ng))
    fill!(st.on, false); fill!(st.κ2, 0.0); for k in 1:mgr; st.fpHe[k].rate[] = 0.0; st.fpHep[k].rate[] = 0.0; end
    k0 = 0
    for k in a:b; (!m.eqzone[k] && !(h !== nothing && tn[k] > h.T1)) && (k0 = k; break); end
    (k0 <= a + 2) && (st.k0 = 0; return nothing)                     # no table region in front of the network: nothing to source
    js = Int(TF.getscalar(:jstep))
    (st.G === nothing || js - st.step >= IFSN || k0 != st.k0) && ifs_source!(st, m, a, b, k0, tauz)
    for k in a:b; st.xabs0[k] = ifs_x(CNn.full(m.x[k], m.net)); st.xabs[k] = st.xabs0[k]; end
    ifs_sweep!(st, k -> st.xabs0[k], a, b)
    last = 0
    for k in k0:b
        V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
        sum(st.nin[k]) * dt / (dn[k] / CNn.μH() * V) > IFSXMIN && (last = k)
    end
    for k in k0:min(b, max(last + 2, k0))
        st.on[k] = true
        slrates === nothing || slrates(k)
        st.p3[k] = (slrates !== nothing && hasproperty(STARLIGHT[], :tl3) && STARLIGHT[].tl3) ? SLDIAG[].tl3[k] : nothing
        ifs_rates!(st, k, st.xabs0[k], tn[k])
        st.r0[k] = (st.fpHe[k].rate[], st.fpHep[k].rate[]); st.κ20[k] = st.κ2[k]
    end
    for k in a:b; st.nin0[k] .= st.nin[k]; end
    st.niter = 0; st.epsit = 0.0
    if IFSDBG[] > 0
        IFSDBG[] -= 1; ng = length(st.N); i1 = findall(gg -> IFS_EDGES[gg] >= 24.587, 1:ng)
        @printf("IONFSDBG step %d k0 %d dt %.3e: zone r[R☉] T n_H x_H x_He x_He+ x_e | Q_H Q_HeI in | absorbed R_He R_He+ κ₂ on\n", js, k0 - a + 1, dt)
        for k in k0:min(b, k0 + 40)
            x = st.xabs0[k]
            @printf("  %3d %9.3e %8.0f %8.2e %8.2e %8.2e %8.2e %6.3f | %9.2e %9.2e | %9.2e %9.2e %9.2e %9.2e %d\n", k - a + 1, rn[k] / 6.957e10, tn[k], dn[k] / CNn.μH(),
                    x[1], x[3], x[4], x[6], sum(st.nin[k]), sum(st.nin[k][i1]), st.absd[k], st.fpHe[k].rate[], st.fpHep[k].rate[], st.κ2[k], st.on[k])
        end
    end
    nothing
end
"After each Newton update (COUPLED): sweep from the iterate, refresh the ray zones' rates, drop the cached Jacobians that changed."
function ifs_iter!(m)
    st = IFS[]; (st === nothing || st.k0 == 0 || !IFSITER) && return nothing
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); tn = TM.v(:tn)
    old = [copy(st.nin[k]) for k in a:b]
    ifs_sweep!(st, k -> TX.cplz(k) ? st.xabs[k] : st.xabs0[k], a, b)
    dmax = 0.0
    for k in a:b
        st.on[k] || continue
        o = old[k-a+1]; r1 = (st.fpHe[k].rate[], st.fpHep[k].rate[], st.κ2[k])
        ifs_rates!(st, k, TX.cplz(k) ? st.xabs[k] : st.xabs0[k], tn[k])
        r2 = (st.fpHe[k].rate[], st.fpHep[k].rate[], st.κ2[k])
        dn_ = maximum(abs(st.nin[k][g] - o[g]) / (abs(o[g]) + 1e-30 * st.N[1] + 1e-300) for g in eachindex(o))
        dr_ = maximum(abs(r2[i] - r1[i]) / (abs(r1[i]) + 1e-300) for i in 1:3)
        dmax = max(dmax, min(dn_, 1e3))
        (dn_ > 1e-10 || dr_ > 1e-10) && k <= length(CPLCACHE[]) && (CPLCACHE[][k] = nothing)
    end
    st.niter += 1; st.epsit = dmax
    nothing
end
"A failed step's retry (COUPLED reset): back to the start-of-step ray."
function ifs_reset!()
    st = IFS[]; (st === nothing || st.k0 == 0) && return nothing
    for k in eachindex(st.nin)
        st.nin[k] .= st.nin0[k]
        if st.on[k]
            st.fpHe[k].rate[], st.fpHep[k].rate[] = st.r0[k]; st.κ2[k] = st.κ20[k]; st.xabs[k] = st.xabs0[k]
            k <= length(CPLCACHE[]) && (CPLCACHE[][k] = nothing)
        end
    end
    nothing
end
"Post pass, before the zones: the sweep from the solved coupled zones (split zones: start of step), and the convergence check ε."
function ifs_post!(m, a, b)
    st = IFS[]; (st === nothing || st.k0 == 0) && return nothing
    tn = TM.v(:tn); used = copy(st.absd)
    for k in a:b; TX.cplz(k) && (st.xabs[k] = ifs_x(CNn.full(m.x[k], m.net))); end
    ifs_sweep!(st, k -> TX.cplz(k) ? st.xabs[k] : st.xabs0[k], a, b)
    e = 0.0
    for k in a:b
        st.on[k] || continue
        ifs_rates!(st, k, TX.cplz(k) ? st.xabs[k] : st.xabs0[k], tn[k]); e += abs(st.absd[k] - used[k])
    end
    st.eps = st.Qsrc[1] > 0 ? e / st.Qsrc[1] : 0.0
    nothing
end
"A table zone beyond k0 that the ray would change goes to the network (it then takes the ray's terms)."
function ifs_force!(k, dt)
    st = IFS[]; (st === nothing || st.k0 == 0 || k < st.k0) && return false
    dn, rn = TM.v(:dn), TM.v(:rn); V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
    sum(st.nin[k]) * dt / (dn[k] / CNn.μH() * V) > IFSXMIN || return false
    st.on[k] || (st.on[k] = true; st.p3[k] = nothing; ifs_rates!(st, k, st.xabs0[k], TM.v(:tn)[k]))
    true
end
const IFS_ENVB = Ref{Any}(nothing)
"The ray's terms of zone k composed with the zone's own (extra, cooling, env, photos)."
function ifs_terms(m, k, exs, cools, envk, phk)
    st = IFS[]; dn, rn = TM.v(:dn), TM.v(:rn)
    src = Main.IonRay.source(dn[k], rn[k], rn[k+1], st.nin[k], st.G, m.net)
    iHp = findfirst(==(:Hp), m.net.free)
    exr = x -> (f = src.extra(x); d = CNn.full(x, m.net); xp = max(d[:Hp], zero(d[:Hp])); f[iHp] -= st.κ2[k] * xp * (1 - xp); f)
    ex2 = exs === nothing ? exr : (x -> exs(x) .+ exr(x))
    c2 = (ρ, T, d) -> (cools === nothing ? zero(T) : cools(ρ, T, d)) - src.heat(d)
    IFS_ENVB[] === nothing && (IFS_ENVB[] = merge(m.env, (; caseB = true)))
    env2 = !IFSENV ? envk : hasproperty(envk, :ljW) ? merge(envk, (; caseB = true, ljcut = IFS_CUT)) :
           (hasproperty(envk, :static) && envk.static) ? merge(envk, (; caseB = true, tabkey = IFS_ENVB[])) :
           (J0 = envk.J; merge(envk, (; caseB = true, J = ν -> CNn.hP * ν >= IFS_CUT ? 0.0 : J0(ν))))
    # He stages: backward Euler from the start of the step (ChemNet.ck_expand heold), not their instantaneous equilibrium
    x0 = st.xabs0[k]; yHe_ = CNn.yHe; fHe0 = (x0[3] / yHe_, x0[4] / yHe_, max(1 - (x0[3] + x0[4]) / yHe_, 0.0))
    env2 = merge(env2, (; heold = () -> (fHe0..., Float64(TF.getscalar(:dtime)))))
    ph2 = Any[phk..., st.fpHe[k], st.fpHep[k]]
    ex2, c2, env2, ph2
end
"After the zones' chemistry: the ray's energy leaves the grey field (RAYQ), diagnostics (ionfs.csv)."
function ifs_after!(m, a, b, dt)
    st = IFS[]; (st === nothing || st.k0 == 0) && return nothing
    dn, rn, tn = TM.v(:dn), TM.v(:rn), TM.v(:tn); G = st.G
    for k in a:b
        st.on[k] || continue
        x = ifs_x(CNn.full(m.x[k], m.net)); nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]; V = 4π / 3 * (rn[k+1]^3 - rn[k]^3); E = 0.0
        for gg in eachindex(st.nin[k])
            τ = 0.0; for j in 1:4; τ += x[j] * nH * G.σ[gg, j] * dr; end
            E += st.nin[k][gg] * (-expm1(-τ)) * G.Ec[gg]
        end
        isempty(TX.RAYQ[]) || (TX.RAYQ[][k] += E / V)
    end
    nothing
end
"Q_H, Q_HeI, Q_HeII through the interfaces a..b+1 (the last sweep), the front zone and radius."
function ifs_profile(m, a, b)
    st = IFS[]; (st === nothing || st.k0 == 0) && return nothing
    ng = length(st.N); iHe1 = findall(gg -> IFS_EDGES[gg] >= 24.587, 1:ng); iHe2 = findall(gg -> IFS_EDGES[gg] >= 54.418, 1:ng)
    dn, rn = TM.v(:dn), TM.v(:rn)
    Q = zeros(b - a + 2); Q1 = zeros(b - a + 2); Q2 = zeros(b - a + 2)
    Nout = ifs_sweep!(st, k -> ifs_x(CNn.full(m.x[k], m.net)), a, b)
    for k in a:b
        Q[k-a+1] = sum(st.nin[k]); Q1[k-a+1] = sum(st.nin[k][iHe1]); Q2[k-a+1] = sum(st.nin[k][iHe2])
    end
    Q[end] = sum(Nout); Q1[end] = sum(Nout[iHe1]); Q2[end] = sum(Nout[iHe2])
    kf = 0; rf = NaN
    for k in st.k0:b
        if Q[k-a+1] >= 0.5 * st.Qsrc[1] && Q[k-a+2] < 0.5 * st.Qsrc[1]
            kf = k; xp = max(CNn.full(m.x[k], m.net)[:Hp], 0.0); rf = cbrt(rn[k]^3 + min(xp, 1.0) * (rn[k+1]^3 - rn[k]^3)); break
        end
    end
    (; k0 = st.k0, rb = st.rb, Q, QHe1 = Q1, QHe2 = Q2, kf, rf, Qsrc = st.Qsrc, fnorm = st.fnorm, Lfs = st.Lfs, Lrun = st.Lrun, eps = st.eps,
       epsit = st.epsit, niter = st.niter, κ2 = [st.κ2[k] for k in a:b], on = [st.on[k] for k in a:b])
end
