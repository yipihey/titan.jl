# S8a: coupling of the Monte Carlo Lyα transfer (LyaMC.jl) to the S7b feedback path (RT_LYA=mc). Loaded by s7b_env.jl (module S8a in Main).
#
#   inner column   Lyα of the H II region (0.68 per absorbed ionising photon) and the stellar continuum in the window |ν/ν_α − 1| < RT.LYA_MC_WIN pass the neutral free-fall column
#                  between the H II region and the first grid interface (τ_0 ~ 1e13, aτ_0 ~ 1e10: out of reach of a direct Monte Carlo). The column is transported in a *scaled* medium:
#                  in the diffusive regime aτ_0 ≫ 1 the transfer equation (Harrington 1973) is invariant under x → x/λ, τ_0 → τ_0/λ³ (the escape frequency is x_esc ∝ (aτ_0)^{1/3}),
#                  so the column density is divided by s = λ³, the destruction probability per scattering multiplied by s (the number of scatterings of a photon is ∝ τ_0 in the thick limit,
#                  ε N_sc invariant: test s8a_test.jl), the recoil parameter by λ (g x_esc invariant), the flow velocities in Doppler units divided by λ, and the spectra rescaled by λ.
#                  The Monte Carlo yields response functions (emergent spectrum per injected photon) of the line and of NWIN sub-windows of the continuum, cached while the column τ_0 stays
#                  within 25 %. With ε τ_0 ≫ 1 (two-photon channel, ε N_sc, N_sc ≈ 1.6 τ_0) all photons that get into the core are destroyed; only those injected beyond ≈ x_esc escape.
#   grid           LyaMC on the hydro zones (moving medium, composition of the sweep): sources = emergent spectrum of the inner column at the first interface, recombination and collisional
#                  Lyα of the zones; inner boundary = reflecting with the albedo of the thick column; outputs: force per zone (-> momentum row), J_α, T_c, escape, destruction.
module S8a
using Printf
const L = Main.LyaMC; const RT = Main.RadTransfer; const CN = Main.ChemNet
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
const TAUMC = envf("S8A_TAUMC", 2e6); const NPK_LINE = Int(envf("S8A_NPK_LINE", 320)); const NPK_WIN = Int(envf("S8A_NPK_WIN", 24)); const NPK_GRID = Int(envf("S8A_NPK_GRID", 300)); const NPK_IN = Int(envf("S8A_NPK_IN", 300))
const NBATCH = Int(envf("S8A_NBATCH", 8)); const NEVERY = Int(envf("S8A_NEVERY", 25)); const XE_IN = envf("S8A_XE_IN", 1e-4); const NZIN = Int(envf("S8A_NZIN", 6))
const EPSTAU_MAX = envf("S8A_EPSTAU_MAX", 40.0); const GRIDSRC = get(ENV, "S8A_GRIDSRC", "all"); const NWIN = 24; const SAFE = envf("S8A_SAFE", 4.0)
const mHyd = CN.mH

mutable struct Inner
    valid::Bool; tau0::Float64; resp_line::Vector{Float64}; resp_win::Matrix{Float64}; wedges::Vector{Float64}; lam::Float64
    flines::Float64; fwin::Vector{Float64}; nsc::Float64; calls::Int; t_cpu::Float64; Ncol::Float64; T::Float64; stats::NamedTuple
end
const INNER = Inner(false, 0.0, Float64[], zeros(0, 0), Float64[], 1.0, 0.0, Float64[], 0.0, 0, 0.0, 0.0, 0.0, (;))
const LASTGRID = Ref{Any}(nothing)          # (rc, accel) of the last grid Monte Carlo (held between calls)
const LOG = Ref{Any}(nothing)

"Column density N_H of the free-fall column ρ ∝ r^{-3/2} from `r1` to `r2` normalised to n_H(r2) = nH2."
column(nH2, r1, r2) = 2 * nH2 * r2^1.5 * (r1^-0.5 - r2^-0.5)

"""
    inner_response!(Rstar, RII, redge, nHedge, uedge, T) -> Inner

Response functions of the free-fall column (R_II → r_e, ρ ∝ r^{-3/2}, inflow ∝ r^{-1/2}) by the scaled Monte Carlo, cached in `INNER`.
"""
function inner_response!(Rstar, RII, redge, nHedge, uedge, T)
    Ncol = column(nHedge * (1 - XE_IN), RII, redge); tau0 = Ncol * L.sig0(T)
    epsin = L.eps_collisional(T, XE_IN * nHedge, XE_IN * nHedge)
    if INNER.valid && abs(log(tau0 / INNER.tau0)) < log(1.25) && abs(log(T / INNER.T)) < 0.2
        return INNER
    end
    t0 = time()
    s = max(tau0 / TAUMC, 1.0); lam = cbrt(s)
    nz = NZIN
    r = vcat([Rstar], RII .* (redge / RII) .^ range(0, 1, length = nz + 1))
    rc = sqrt.(r[2:end-1] .* r[3:end])
    nH = vcat(0.0, nHedge .* (rc ./ redge) .^ -1.5)                       # H II zone [R*, R_II]: no neutral H
    nHI = nH .* (1 - XE_IN) ./ s
    # flow: inward, u = -u_e (r/r_e)^(-1/2) (Doppler units / λ in the scaled medium)
    u = -uedge .* (r ./ redge) .^ -0.5 ./ lam
    Z = L.build_zones(r, nHI, fill(T, nz + 1), u; ne = nH .* XE_IN, np = nH .* XE_IN, eps_extra = s, g_extra = lam, wscale = lam, safe = SAFE)
    a = L.damp(T); b = L.bth(T)
    xesc = max((a * tau0 / s)^(1 / 3), 3.0)
    wmax = 14 * xesc * b
    srcs = [L.src_sphere(RII, 1.0; npk = NPK_LINE, name = :line)]
    win = RT.LYA_MC_WIN[]; edges = range(-win * L.c, win * L.c, length = NWIN + 1)             # sub-windows in w = c Δν/ν
    for i in 1:NWIN
        push!(srcs, L.src_shell(Rstar, 1.0, [edges[i] / lam, edges[i+1] / lam], [1.0]; npk = NPK_WIN, name = Symbol("win$i")))
    end
    res = L.run(Z, srcs; nbatch = NBATCH, seed = 11, opts = L.Opts(inner = :absorb, wmax = wmax, nbw = 700), store = false)
    INNER.resp_line = res.spec_out_src[1, :]; INNER.resp_win = res.spec_out_src[2:end, :]
    INNER.wedges = collect(range(-wmax, wmax, length = 701)) .* lam
    INNER.lam = lam; INNER.tau0 = tau0; INNER.valid = true; INNER.calls += 1; INNER.T = T; INNER.Ncol = Ncol
    INNER.flines = sum(INNER.resp_line); INNER.fwin = vec(sum(INNER.resp_win, dims = 2)); INNER.nsc = res.nsc_per_photon
    INNER.t_cpu += time() - t0
    INNER.stats = (; s, lam, tau0, tau_mc = tau0 / s, xesc, f_out = res.f_out, f_in = res.f_in, f_2g = res.f_2g, nsc = res.nsc_per_photon, wall = time() - t0, nz = nz,
                   nskip = count(>(0), Z.xc), eps = Z.eps[2] / 1.0)
    @printf("S8A inner column: N_H = %.2e, tau0 = %.2e, scale s = %.2e (lambda %.1f), x_esc(scaled) = %.1f, line escape %.4f (star-absorbed %.4f), N_sc/photon %.2e, %.1f s\n",
            Ncol, tau0, s, lam, xesc, INNER.flines, res.f_in, res.nsc_per_photon, time() - t0)
    INNER
end

"Emergent spectrum of the column: (edges [cm/s], photons/s per bin) from the line rate and the continuum rates of the NWIN sub-windows. The cached response is stretched by (τ_0/τ_0,cache)^{1/3}."
function emergent(inn::Inner, rate_line, rate_win, tau0_now)
    st = cbrt(tau0_now / inn.tau0)
    wts = rate_line .* inn.resp_line
    for i in 1:NWIN; wts .+= rate_win[i] .* view(inn.resp_win, i, :); end
    (inn.wedges .* st, wts)
end

"Albedo of the thick column seen from outside: photons that enter it are returned with probability 1 − 1/(1 + 3τ(w)/4) (two-stream estimate), τ(w) = N_H σ(w)."
function column_albedo(Ncol, T)
    a = L.damp(T); b = L.bth(T); s0 = L.sig0(T)
    w -> (x = w / b; τ = Ncol * s0 * L.voigt(a, x) * L.Cwing(w / L.c); 1 - 1 / (1 + 0.75 * τ))
end

"""
    window_rates(g, Lout, win) -> Vector{Float64}

Photon rates [s⁻¹] of the stellar spectrum `Lout` (power per frequency bin) in NWIN equal sub-windows of |ν/ν_α − 1| < win.
"""
function window_rates(g, Lout, win)
    ν0 = L.ν0; rates = zeros(NWIN)
    for f in 1:RT.nfreq(g)
        Δ = g.ν[f] / ν0 - 1; abs(Δ) < win || continue
        i = clamp(floor(Int, (Δ + win) / (2win) * NWIN) + 1, 1, NWIN)
        rates[i] += Lout[f] / (L.hP * g.ν[f])
    end
    rates
end

"""
    grid_mc!(...) -> (accel [cm s⁻²] per zone, result)

Monte Carlo on the hydro zones: `rn` interfaces, `tn`, `dn` (g cm⁻³), `un` (cm/s, > 0 outwards), `X` composition (per H nucleus), `emerge` = (edges, rates) entering at the first
interface, `Ncol` = H I column of the inner column (albedo), `mz` zone masses [g].
"""
function grid_mc!(rn, tn, dn, un, X, emerge, Ncol, Tin, mz; seed = 1)
    n = length(tn); nH = dn ./ CN.μH()
    g(s) = haskey(X, s) ? X[s] : zeros(n)
    nHI = nH .* g(:H); ne = nH .* g(:e); np = nH .* g(:Hp); nH2 = nH .* g(:H2); nHm = nH .* g(:Hm)
    λHm = 1.2157e-1
    kabs = [nHm[k] * L.sigma_Hm(λHm) for k in 1:n]
    Z = L.build_zones(rn, nHI, tn, un; ne = ne, np = np, kabs = kabs, nH2 = nH2, safe = SAFE)
    V = Z.vol
    rates = grid_rates(tn, dn, X, V)
    srcs = L.Source[]
    wedges, wts = emerge
    Rin = sum(wts)
    Rin > 0 && push!(srcs, L.src_shell(rn[1], Rin, wedges, wts; npk = NPK_IN, name = :inner))
    sg = sum(rates)
    sg > 0 && push!(srcs, L.src_zones(rates; npk = NPK_GRID, name = :grid))
    isempty(srcs) && return (zeros(n), nothing)
    wmax = max(abs(wedges[1]), abs(wedges[end]), 3e8)
    res = L.run(Z, srcs; nbatch = NBATCH, seed = seed, opts = L.Opts(inner = :reflect, albedo = column_albedo(Ncol, Tin), wmax = wmax, nbw = 1200))
    @printf("S8A grid MC: %d zones (%d with core skipping), sources %s, N_sc/photon %.2e, out %.3f in %.3f 2γ %.3f H⁻ %.2e H2 %.2e cut %.3f, F_tot c/L_src %.3f, %d packets\n", n, count(>(0), Z.xc),
            string([(s.name, round(s.rate, sigdigits = 3)) for s in srcs]), res.nsc_per_photon, res.f_out, res.f_in, res.f_2g, res.f_Hm, res.f_H2, res.f_lost, sum(res.force) * L.c / (res.emit * L.E0), sum(s.npk for s in srcs))
    (res.force ./ mz, res)
end

"Lyα photon rates [s⁻¹] of the grid zones: case-B recombination (Cantalupo et al. 2008 fraction) + collisional 1s → 2p excitation by electrons (`S8A_GRIDSRC` = all | rec | none)."
function grid_rates(tn, dn, X, V)
    n = length(tn); nH = dn ./ CN.μH(); g(s) = haskey(X, s) ? X[s] : zeros(n)
    ne = nH .* g(:e); np = nH .* g(:Hp); nHI = nH .* g(:H); rates = zeros(n)
    GRIDSRC == "none" && return rates
    for k in 1:n
        rec = L.f_lya_caseB(tn[k]) * L.alphaB(tn[k]) * ne[k] * np[k]
        col = GRIDSRC == "all" ? ne[k] * nHI[k] * L.gamma_1s2p(tn[k]) : 0.0
        rates[k] = (rec + col) * V[k]
    end
    rates
end

"Interpolate the held acceleration of the last Monte Carlo (zone-centre radii `rc`) onto the zones of the current grid."
function held_accel(rc_now)
    LASTGRID[] === nothing && return zeros(length(rc_now))
    rc0, a0 = LASTGRID[][1], LASTGRID[][2]
    out = zeros(length(rc_now))
    for (i, r) in enumerate(rc_now)
        if r <= rc0[1]; out[i] = a0[1]
        elseif r >= rc0[end]; out[i] = a0[end]
        else
            j = searchsortedlast(rc0, r); w = log(r / rc0[j]) / log(rc0[j+1] / rc0[j]); out[i] = (1 - w) * a0[j] + w * a0[j+1]
        end
    end
    out
end

end # module
