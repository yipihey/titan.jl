# Emergent spectra of recorded snapshots (RUNREC streams, see runrec.jl): post-processing, nothing here feeds back into a run.
#
#   julia -t 8 --project=<env> snapspec.jl <run dir> [stepmod]            # fills <run dir>/spec/<step>_<t>.jls, writes export/spec.{json,f32}
#   julia -t 8 --project=<env> snapspec.jl --old <snaps file> <out dir> [stepmod]
#
# For every snapshot two spectra, each at two places:
#
#   where     the outer mesh boundary ("edge": what leaves the computed volume; gas beyond the mesh is not there)
#             the photosphere of the run (the interface where the run's own grey optical depth reaches 2/3, see `greytau`)
#   what      continuum  L_ν  from a formal solution of the static spherical transfer equation on the snapshot's ρ, T and composition
#             H₂ lines   the luminosity of each of the run's quadrupole lines, from the run's own line transfer (H2Lines, ALI)
#
# Continuum (module Spectrum: tangent rays, Feautrier, accelerated Λ-iteration):
# - opacity: ContinuumOpacity.alpha with the snapshot's own number densities (the network's H₂, H⁺, e⁻, H⁻; not equilibrium):
#   collision-induced absorption H₂–H₂, H₂–He, H₂–H; H⁻ bound–free and free–free; H bound–free and H⁺ free–free (Kramers);
#   Rayleigh scattering on H₂ and H, Thomson scattering. No atomic lines (no Lyα, no Balmer lines), no He opacity.
# - source: S_ν = ε_ν S_th + (1 − ε_ν) J_ν, coherent isotropic scattering; thermal part S_th = B_ν(T_gas) except for the bound–free
#   continua of H⁻ and H, which emit by radiative attachment ∝ n_H n_e and recombination ∝ n_e n_p of the zone's own electrons
#   (Kirchhoff–Milne with Maxwellian electrons; for H⁻ as in the run's κ_P, PopIIIMaterial.nlte_means). With LTE emission instead
#   (κ_ν B_ν with Boltzmann levels on the network's neutral H) the collisionally ionised gas outside the photosphere of the
#   accreting star radiated 16 × the run's luminosity in the Balmer continuum. Absorption by the excited levels of H: the
#   run's three-level populations are not in the snapshot; the smaller of Saha's (zone's n_e n_p) and Boltzmann's (zone's n = 1).
# - static: velocities ignored (|u|/c < 1e-3). Zones deeper than grey optical depth TAUCUT (300) are replaced by a thermal inner
#   boundary I⁺ = B_ν(T) there; with less than that to the centre the whole sphere is solved (rays symmetric about the midplane).
# - the temperature structure is the run's (grey flux-limited diffusion, or no continuum at all before the radiation rebuild): the
#   frequency integral of this spectrum need not equal the run's luminosity. Both are stored, the ratio is the check.
#
# H₂ lines: H2Lines.cooling in ALI mode on the snapshot (level populations and transfer as the run's reference solution; the run
# itself uses the escape-probability solution rescaled to ALI every ALIINT steps), with the run's merging of the zones deeper than
# continuum optical depth LTAUCUT. The line transfer has no continuum extinction: line photons emitted inside the photosphere are
# excluded by the cut, those emitted outside travel freely.
using Serialization, Printf
isdefined(Main, :lineage) || include(joinpath(@__DIR__, "runrec.jl"))
isdefined(Main, :Z0Physics) || include(joinpath(@__DIR__, "Z0Physics.jl"))
isdefined(Main, :Spectrum) || include(joinpath(@__DIR__, "Spectrum.jl"))
isdefined(Main, :ContinuumOpacity) || include(joinpath(@__DIR__, "ContinuumOpacity.jl"))
isdefined(Main, :H2Lines) || include(joinpath(@__DIR__, "H2Lines.jl"))
isdefined(Main, :H2Levels) || include(joinpath(@__DIR__, "H2Levels.jl"))

module SnapSpec
import ..Spectrum, ..ContinuumOpacity, ..H2Lines, ..H2Levels
using Printf
const CO = ContinuumOpacity; const SP = Spectrum
const yHe = 1 / 12; const mH = 1.6726e-24; const μH = 1.6735e-24 * (1 + 4yHe)       # mass per H nucleus as PrimordialChem
const Rs = 6.957e10; const Ls = 3.83e33
const TAUCUT = 300.0; const LTAUCUT = 10.0; const NBIN = 30
"Wavelength grid: log-uniform in wavenumber, 1000 μm to 0.05 μm."
const NU = [CO.c * 10^x for x in range(1.0, log10(2e5), length = 216)]

"""
The run's grey continuum optical depth at the interfaces, as the driver's (`starlight_setup`, `linecool`) and runrec's:
τ(r_k) = Σ_{j ≥ k} χ_j ρ_j Δr_j from the outer mesh boundary (τ = 0 there) inwards, with χ the Rosseland mean of the extinction
(absorption + scattering) the run used in its flux (TITAN `chifn`). Zero where the run had no radiation yet.
"""
function greytau(x)
    n = length(x.ρ); τ = zeros(n + 1)
    for k in n:-1:1; c = x.chi[k]; τ[k] = τ[k+1] + (isfinite(c) && c > 0 ? c : 0.0) * x.ρ[k] * (x.r[k+1] - x.r[k]); end
    τ
end

spc(x, s) = (i = findfirst(==(s), x.spn); i === nothing ? nothing : max.(x.Y[i, :], 0.0))
"Number densities per zone [cm⁻³] from the recorded fractions per H nucleus."
function densities(x)
    nH = x.ρ ./ μH; z = zeros(length(nH)); g(s) = (v = spc(x, s); v === nothing ? z : v)
    xp = g(:Hp); xm = g(:Hm); x2 = g(:H2); x2p = g(:H2p)
    xH = (v = spc(x, :H); v === nothing ? max.(1 .- xp .- xm .- 2 .* x2 .- 2 .* x2p, 0.0) : v)
    xe = (v = spc(x, :e); v === nothing ? max.(xp .- xm .+ x2p, 0.0) : v)
    (; nH = nH .* xH, nH2 = nH .* x2, nHe = yHe .* nH, nHp = nH .* xp, ne = nH .* xe, nHm = nH .* xm)
end

struct Opac
    T::Vector{Float64}; nH::Vector{Float64}; nH2::Vector{Float64}; nHe::Vector{Float64}; nHp::Vector{Float64}; ne::Vector{Float64}; nHm::Vector{Float64}
end
Opac(atm, d) = (z = atm.zone; Opac(copy(atm.T), d.nH[z], d.nH2[z], d.nHe[z], d.nHp[z], d.ne[z], d.nHm[z]))
const E0HM = 0.754 * CO.eV; const χH = 13.598 * CO.eV; const σ0H = 7.906e-18
const QB = sqrt(NU[2] / NU[1])                       # a grid frequency ν stands for the band [ν/QB, ν QB]
"""
(absorption, scattering, thermal source) at point k for the band of frequency ν. Everything of ContinuumOpacity.alpha except
hydrogen's bound–free continua emits B_ν(T_gas) × its absorption. The bound–free continua emit by recombination (attachment for
H⁻) of the zone's own electrons: j_ν = σ_n(ν) n_e n_p n² λ_e³ e^{−(hν − χ_n)/kT} 2hν³/c² (Kirchhoff–Milne, Maxwellian electrons),
whatever the level populations are. They absorb with the zone's own n = 1 (and H⁻); the excited levels, which the snapshot does not
hold, are the smaller of Saha's with the zone's n_e n_p and Boltzmann's with its n = 1 (equal in LTE). Stimulated recombination
is subtracted from the absorption.

Recombination emission falls as e^{−(hν − χ_n)/kT} above each threshold, in cold gas within much less than a band (4.7 % in ν):
bound–free emission and absorption are band averages (exact for Kramers' σ ∝ ν⁻³; H⁻ by midpoint sums that resolve kT).
"""
function SP.coeffs3(o::Opac, k, ν)
    T = o.T[k]; kT = CO.kB * T; hν = CO.h * ν
    a, s, p = CO.alpha(ν, T; nH = o.nH[k], nH2 = o.nH2[k], nHe = o.nHe[k], nHp = o.nHp[k], ne = o.ne[k], nHm = o.nHm[k])
    B = CO.Bν(ν, T); stim = -expm1(-hν / kT); c2 = 2 * CO.h / CO.c^2
    ab = a - p.hm_bf - p.h + 3.692e8 * o.ne[k] * o.nHp[k] / (sqrt(T) * ν^3) * stim       # CIA, H⁻ ff, H⁺ ff
    em = ab * B
    νl = ν / QB; νh = ν * QB; Δν = νh - νl
    λe3 = (CO.h^2 / (2π * CO.me * kT))^1.5
    # H⁻ bound–free
    ν0 = E0HM / CO.h
    if νh > ν0 && (o.nHm[k] > 0 || o.nH[k] * o.ne[k] > 0)
        νa = max(νl, ν0); m = clamp(ceil(Int, 4 * CO.h * (νh - νa) / kT), 1, 64); dν = (νh - νa) / m; sa = 0.0; se = 0.0
        x0 = (CO.h * νa - E0HM) / kT
        if x0 < 200 || o.nHm[k] > 0
            for i in 1:m
                νi = νa + (i - 0.5) * dν; σ = CO.σ_Hm_bf(1e4 * CO.c / νi)
                sa += σ * (-expm1(-CO.h * νi / kT)); se += σ * νi^3 * exp(-(CO.h * νi - E0HM) / kT)
            end
            ab += o.nHm[k] * sa * dν / Δν; em += o.nH[k] * o.ne[k] * λe3 / 4 * c2 * se * dν / Δν
        end
    end
    # H bound–free, levels 1…30
    nep = o.ne[k] * o.nHp[k]; n1 = o.nH[k]
    if nep > 0 || n1 > 0
        lnep = nep > 0 ? log(nep * λe3) : -Inf; ln1 = n1 > 0 ? log(n1) : -Inf
        for n in 1:30
            χn = χH / n^2; νn = χn / CO.h; νh > νn || continue
            νa = max(νl, νn); σν3 = σ0H * n * νn^3                                      # σ_n(ν) ν³
            # ∫ e^{−(hν − χ_n)/kT} dν over the band, with the Saha factor: n² n_e n_p λ_e³ e^{χ_n/kT} e^{−hν/kT}
            ea = lnep - (CO.h * νa - χn) / kT
            I = ea > -700 ? n^2 * exp(ea) * (kT / CO.h) * (-expm1(-CO.h * (νh - νa) / kT)) : 0.0
            nn = n == 1 ? n1 : n^2 * exp(min(lnep + χn / kT, ln1 - (χH - χn) / kT))
            νe = min(νa + kT / CO.h, 0.5 * (νa + νh))                                   # where the stimulated term sits in the band
            ab += σν3 * max(nn * (νa^-2 - νh^-2) / 2 - I / νe^3, 0.0) / Δν; em += σν3 * c2 * I / Δν
        end
    end
    ab > 0 ? (ab, s, em / ab) : (ab, s, B)
end

"""
Continuum of snapshot x on the grid NU: (Le, Lp, rτ, τin, its) — L_ν [erg s⁻¹ Hz⁻¹] leaving the outermost zone, L_ν = 16π² r² H_ν through
the photosphere interface kph (NaN without one), the radius [cm] where the monochromatic radial optical depth (extinction, from
the outer boundary) reaches 2/3 (NaN where the whole solved region is thinner), that depth at the inner end of the solved region,
and the iteration counts (0: scattering neglected, −1: direct solution with the full Λ matrix).
"""
function continuum(x, τ, kph)
    n = length(x.ρ); rc = [0.5 * (x.r[k] + x.r[k+1]) for k in 1:n]
    kc = findlast(>(TAUCUT), τ); kc = kc === nothing ? 1 : min(kc, n - 8)       # first zone kept
    atm = SP.atmosphere((r = x.r, d = x.ρ, T = x.T); rmin = kc > 1 ? rc[kc] : nothing, split = 3.0)
    op = Opac(atm, densities(x))
    _, R = SP.rays(atm.r; ncore = kc > 1 ? 8 : 1)
    kc > 1 || (R = [merge(ray, (core = false,)) for ray in R])                   # whole sphere: no thermal core behind the rays
    R = SP.weights!(atm, R); ps = [ray.p for ray in R]
    np = length(atm.r); nν = length(NU)
    Le = zeros(nν); Lp = fill(NaN, nν); rτ = fill(NaN, nν); τin = zeros(nν); its = zeros(Int, nν)
    rph = kph === nothing ? NaN : x.r[kph]
    ip = (isfinite(rph) && atm.r[1] < rph < atm.r[end]) ? searchsortedlast(atm.r, rph) : 0
    CO.tables()
    Threads.@threads for i in 1:nν
        ν = NU[i]; κ = zeros(np); σ = zeros(np); B = zeros(np)
        for k in 1:np; κ[k], σ[k], B[k] = SP.coeffs3(op, k, ν); end
        χ = max.(κ, 1e-300) .+ σ
        tν = zeros(np); ts = 0.0
        for k in np-1:-1:1
            dr = atm.r[k+1] - atm.r[k]; tν[k] = tν[k+1] + SP.dtau(χ[k], χ[k+1], dr, atm.jump[k]); ts += 0.5 * (σ[k] + σ[k+1]) * dr
        end
        τin[i] = tν[1]
        k3 = findlast(>(2 / 3), tν)
        if k3 !== nothing && k3 < np
            w = (tν[k3] - 2 / 3) / (tν[k3] - tν[k3+1]); rτ[i] = exp((1 - w) * log(atm.r[k3]) + w * log(atm.r[k3+1]))
        end
        if maximum(B) <= 0; continue; end
        local H
        if ts < 1e-3          # scattering depth negligible: S = ε S_th (the scattered light, a fraction ≲ ts of the emission, is dropped)
            S = max.(κ, 1e-300) .* B ./ χ; Ib = B[1]
            Io = [SP.emergent(ray.z, χ[ray.idx], S[ray.idx], ray.core, Ib; jm = ray.jm) for ray in R]
            L = 0.0; for j in 1:length(R)-1; L += 0.5 * (Io[j] * ps[j] + Io[j+1] * ps[j+1]) * (ps[j+1] - ps[j]); end
            Le[i] = 8π^2 * L
            ip > 0 && ((_, H) = SP.sweep_moments(atm, R, χ, S, Ib))
        else
            # (scattering-dominated wavelengths, Rayleigh in the blue: the Λ-iteration stalls and the full Λ matrix is solved, its = −1)
            Le[i], _, it, m = SP.monochromatic(atm, op, R, ν; maxit = 60, moments = true, maxit_direct = true)
            its[i] = it; H = m.H
        end
        if ip > 0
            w = (rph - atm.r[ip]) / (atm.r[ip+1] - atm.r[ip])
            Lp[i] = 16π^2 * ((1 - w) * atm.r[ip]^2 * H[ip] + w * atm.r[ip+1]^2 * H[ip+1])
        end
    end
    (; Le, Lp, rτ, τin, its, kc, np)
end

const HLL = Ref{Any}(nothing); const HLV = Ref{Any}(nothing)
"""
H₂ line luminosities [erg s⁻¹] per line at the outer boundary and at the photosphere interface (the run's ALI line transfer with
the run's core merging: `linecool` of popiii_m2c.jl), and the transfer's diagnostics.
"""
function h2lines(x, τ, kph)
    HLL[] === nothing && (HLL[] = H2Lines.load(); HLV[] = H2Levels.load())
    n = length(x.ρ); d = densities(x); rr = x.r; uu = x.u; dd = x.ρ; tt = x.T
    jc = findlast(>(LTAUCUT), τ[1:n]); jc = jc === nothing ? 0 : jc; jc = jc >= n - 4 ? 0 : jc     # the driver's: outermost zone whose inner interface is deeper than LTAUCUT
    if jc < 2
        jc = 0; nc = 0; r2 = rr; u2 = uu; t2 = tt; h2 = d.nH2; nH = d.nH; nHe = d.nHe; nHp = d.nHp; d2 = dd
    else
        nc = 1; rin = rr[jc+1]; r2 = [0.0; rr[jc+1:end]]; u2 = [0.0; uu[jc+1:end]]
        mc = sum(dd[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:jc); d2 = [mc / (4π / 3 * rin^3); dd[jc+1:end]]
        t2 = [tt[jc]; tt[jc+1:end]]; h2 = [0.0; d.nH2[jc+1:end]]
        nH = [d.nH[jc]; d.nH[jc+1:end]]; nHe = [d.nHe[jc]; d.nHe[jc+1:end]]; nHp = [d.nHp[jc]; d.nHp[jc+1:end]]
    end
    _, Lr, Ll = H2Lines.cooling(HLL[], r2, u2, d2, t2, h2; nbin = NBIN, nlte = (LV = HLV[], nH, nHe, nHp), nlte_mode = :ali)
    la = H2Lines.LAST[]
    kp = kph === nothing ? 0 : kph - jc + nc                        # the photosphere interface on the reduced grid
    # the gas's net collisional loss to the levels at the ALI populations, summed outwards: what the run books as line cooling
    # (LINECOLL) and as L_lines. Equal to the line flux in statistical equilibrium; the formal solution's flux differs by its
    # discretisation error where the lines are thick (up to a factor 2), so both are kept.
    Lc = zeros(length(d2) + 1)
    for j in eachindex(d2)
        Lc[j+1] = Lc[j] + H2Levels.collcool(HLV[], t2[j], view(la.pops, :, j), nH[j], h2[j], nHe[j], nHp[j]) * 4π / 3 * (r2[j+1]^3 - r2[j]^3)
    end
    (; ν = HLL[].ν, Le = Ll[:, end], Lp = kp >= 1 ? Ll[:, kp] : fill(NaN, size(Ll, 1)), Ltot = Lr[end], Lcoll = Lc[end], Lcoll_ph = kp >= 1 ? Lc[kp] : NaN,
       jc, iters = la.iters, converged = la.converged, nthick = la.nthick)
end

"Everything stored per snapshot."
function spectra(x; lines = true)
    τ = greytau(x); kph = findlast(>(2 / 3), τ)
    t0 = time(); c = continuum(x, τ, kph); tc = time() - t0
    t0 = time()
    l = try lines ? h2lines(x, τ, kph) : nothing catch e; @warn "line transfer failed" step = x.step e; nothing end
    tl = time() - t0
    Lrun = 4π .* x.r .^ 2 .* x.fr                                    # the run's own (grey) radiative luminosity at the interfaces
    # (edge: two interfaces inside the boundary; the boundary interface itself carries the free-streaming flux of the background field)
    (; step = x.step, t = x.t, ν = NU, cont = c, lines = l, kph = kph === nothing ? 0 : kph, rph = kph === nothing ? NaN : x.r[kph],
       Lrun_edge = Lrun[end-2], Lrun_ph = kph === nothing ? NaN : Lrun[kph], Teff = x.Teff, tc, tl)
end
key(x) = @sprintf("%09d_%.10e", x.step, x.t)
end # module

"""
Write <out>/spec.json and <out>/spec.f32 from the cached spectra of the snapshots L (surviving branch only; a restart leaves
spectra of abandoned branches in the cache). spec.f32 is Float32, spectrum after spectrum:
  ν L_ν [L_sun] of the continuum at the edge (nν values), at the photosphere (nν), r(τ_ν = 2/3) [R_sun] (nν),
  H₂ line luminosities [L_sun] at the edge (nl) and through the photosphere (nl).
The line values are the formal solution's, scaled by one factor per snapshot (`lscale`) to the collisional total (see h2lines).
spec.json: wavelengths [μm] of the continuum grid and of the lines, and per spectrum its step, time and the scalars SPSCAL
(luminosities in L_sun; `Lc_*` = ∫ L_ν dν of the continuum spectrum, `Lrun_*` = the run's own grey radiative luminosity there).
"""
const SPSCAL = ["step", "t_yr", "kph", "r_ph", "T_eff", "Lrun_edge", "Lrun_ph", "Lc_edge", "Lc_ph", "Ll_edge", "Ll_ph", "Ll_formal", "Ll_run", "lscale",
                "ndirect", "ali_iters"]
function export_spec(L, cache, out)
    mkpath(out); nν = length(SnapSpec.NU); q = SnapSpec.QB - 1 / SnapSpec.QB; rows = Any[]
    λl = 1e4 .* SnapSpec.CO.c ./ H2Lines.load().ν; nl = length(λl)
    open(joinpath(out, "spec.f32.tmp"), "w") do io
        for x in L
            f = joinpath(cache, SnapSpec.key(x) * ".jls"); isfile(f) || continue
            s = try deserialize(f) catch; continue end
            c = s.cont; l = s.lines
            band(Lν) = sum(Lν[j] * SnapSpec.NU[j] * q for j in 1:nν)                    # ∫ L_ν dν, band sums
            sc = (l === nothing || l.Ltot <= 0) ? NaN : l.Lcoll / l.Ltot
            Teff = (s.kph > 0 && s.Lrun_ph > 0) ? (s.Lrun_ph / (4π * s.rph^2 * sb))^0.25 : NaN
            write(io, Float32.(c.Le .* SnapSpec.NU ./ Ls)); write(io, Float32.(c.Lp .* SnapSpec.NU ./ Ls)); write(io, Float32.(c.rτ ./ Rs))
            if l === nothing
                write(io, fill(NaN32, 2nl))
            else
                write(io, Float32.(l.Le .* sc ./ Ls)); write(io, Float32.(l.Lp .* (isfinite(l.Lcoll_ph) && sum(l.Lp) > 0 ? l.Lcoll_ph / sum(l.Lp) : 1.0) ./ Ls))
            end
            push!(rows, [s.step, s.t / yr, s.kph, s.rph / Rs, Teff, s.Lrun_edge / Ls, s.Lrun_ph / Ls, band(c.Le) / Ls, s.kph > 0 ? band(c.Lp) / Ls : NaN,
                         l === nothing ? NaN : l.Lcoll / Ls, l === nothing ? NaN : l.Lcoll_ph / Ls, l === nothing ? NaN : l.Ltot / Ls, x.Llines / Ls, sc,
                         count(==(-1), c.its), l === nothing ? NaN : l.iters])
        end
    end
    js(v) = isfinite(v) ? @sprintf("%.6g", v) : "null"
    open(joinpath(out, "spec.json.tmp"), "w") do io
        print(io, "{\"nnu\":$nν,\"nl\":$nl,\"lam\":[", join((@sprintf("%.6g", 1e4 * SnapSpec.CO.c / ν) for ν in SnapSpec.NU), ","), "],\"lamline\":[",
              join((@sprintf("%.6g", v) for v in λl), ","), "],\"scalars\":[", join(("\"$s\"" for s in SPSCAL), ","), "],\"spectra\":[",
              join(("[" * join(js.(r), ",") * "]" for r in rows), ","), "]}")
    end
    for f in ("spec.f32", "spec.json"); mv(joinpath(out, f * ".tmp"), joinpath(out, f); force = true); end
    println("exported $(length(rows)) spectra to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # snapspec.jl <run dir> [stepmod]   or   snapspec.jl --old <snaps file> <out dir> [stepmod]
    # Spectra of the snapshots whose step number is a multiple of stepmod (default 160), and of the last one; cached in
    # <dir>/spec/<step>_<t>.jls; SPECMAX=<n> stops after n new ones (0: export only). Then the export (<dir>/export/spec.*),
    # also after every 20 new spectra.
    old = ARGS[1] == "--old"
    L = old ? lineage_old(ARGS[2]) : lineage(ARGS[1]); dir = old ? ARGS[3] : ARGS[1]
    smod = length(ARGS) > (old ? 3 : 1) ? parse(Int, ARGS[old ? 4 : 2]) : 160
    cache = joinpath(dir, "spec"); mkpath(cache)
    sel = [i for i in eachindex(L) if L[i].step % smod == 0]; (isempty(sel) || sel[end] != length(L)) && !isempty(L) && push!(sel, length(L))
    # coarse to fine, so that a partly filled cache covers the whole run evenly
    order = Int[]; for s in (64, 32, 16, 8, 4, 2, 1), j in 1:s:length(sel); sel[j] in order || push!(order, sel[j]); end
    nmax = parse(Int, get(ENV, "SPECMAX", "100000")); ndone = 0
    for i in order
        ndone >= nmax && break
        f = joinpath(cache, SnapSpec.key(L[i]) * ".jls"); isfile(f) && continue
        s = SnapSpec.spectra(L[i]); tmp = f * ".tmp"; serialize(tmp, s); mv(tmp, f; force = true)
        @printf("step %d: continuum %.1f s (%d points, cut zone %d, direct solutions %d), lines %.1f s%s\n", s.step, s.tc, s.cont.np, s.cont.kc,
                count(==(-1), s.cont.its), s.tl, s.lines === nothing ? " (failed)" : @sprintf(" (ALI %d, converged %s)", s.lines.iters, s.lines.converged))
        flush(stdout); global ndone += 1
        ndone % 20 == 0 && export_spec(L[sort(sel)], cache, joinpath(dir, "export"))      # a growing export while the cache fills
    end
    println("spectra: $ndone computed, $(count(i -> isfile(joinpath(cache, SnapSpec.key(L[i]) * ".jls")), sel)) of $(length(sel)) selected snapshots cached")
    export_spec(L[sel], cache, joinpath(dir, "export"))
end
