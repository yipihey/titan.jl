# S7b acceptance tests: coupled protostellar feedback (RadSweep.jl in RadTransfer, s7b_feedback.jl, the momentum-row hook TitanResidFast.AEXT).
#   TITAN_NOLIB=1 julia --project=env -t 4 s7b_test.jl [a|b|c|d|all]
# a  Omukai & Inutsuka (2002): H II region in free-fall accretion, R_II = R_in exp[(Q − Ṁ/m_p)/Q_crit]; critical Ṁ for confinement
# b  Omukai & Inutsuka (2002): uniform-velocity solution, ionising radiation force = gravity at Q_crit (eq. 28), R_II,crit (eq. 25)
# c  radiation-pressure shell in TITAN (AEXT hook): thin and thick shell vs analytic momentum balance
# d  multi-species sweep: photon and energy conservation, He ionisation structure
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl")); include(joinpath(J, "RadTransfer.jl")); include(joinpath(J, "Perf.jl"))
const CN = ChemNet; const RT = RadTransfer; const A = AtomicHHe
const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-30s %s\n", ok ? "PASS" : "FAIL", name, detail); flush(stdout)
end
const G = 6.67430e-8; const Msun = 1.98847e33; const Rsun = 6.957e10; const YR = 3.15576e7; const mp = CN.mp; const AU = 1.495978707e13; const c = CN.c; const hP = CN.hP; const eV = CN.eV
A.ONSPOT[] = false
which = isempty(ARGS) ? "all" : ARGS[1]

net, photos = A.atomic()
gr = RT.default_grid(dUV = 0.05)
const IA = RT.ion_atoms(photos, gr; yHe = 0.0, dr = A.k_di_HeII)
const IAHe = RT.ion_atoms(photos, gr; dr = A.k_di_HeII)

"R of x = 1/2 (interpolated linearly in ln r between zone centres), the outermost crossing from inside."
function rfront(x, r)
    n = length(x); rc = [sqrt(r[k] * r[k+1]) for k in 1:n]; kf = findlast(>=(0.5), x)
    (kf === nothing) && return r[1]; kf == n && return r[end]
    w = (x[kf] - 0.5) / (x[kf] - x[kf+1]); exp(log(rc[kf]) + w * (log(rc[kf+1]) - log(rc[kf])))
end

# ====================================================================== a. free-fall flow (OI02 §2)
function test_a()
    M = 30Msun; Rin = 10Rsun; T = 1e4; αB = RT.ion_alpha(IA, T)[1]
    function run_a(Mdot, Q; nz = 600, rmaxfac = 3e4, Eph = 20.0)
        r = Rin .* exp.(range(0, log(rmaxfac), length = nz + 1))
        rc = sqrt.(r[1:end-1] .* r[2:end]); v = sqrt.(2G * M ./ rc); ρ = Mdot ./ (4π .* rc .^ 2 .* v); nH = ρ ./ mp
        src = RT.mono_source(gr, Eph, Q)
        # the flux of H nuclei that crosses the zone (neutral gas arrives from outside), pure H: Ṁ/m_p
        s = RT.ion_sweep(IA, src.L, r, fill(T, nz), nH; mdotH = fill(Mdot / mp, nz), passes = 8)
        (; s, r, nH, RII = rfront(s.x, r))
    end
    Qcrit(Mdot) = αB * Mdot^2 / (8π * G * mp^2 * M)
    @printf("a: M = %.0f Msun, R_in = %.0f Rsun, α_B(1e4 K) = %.3e, Q_crit(1e-3 Msun/yr) = %.3e s^-1 (OI02 eq. 9: %.2e)\n", M / Msun, Rin / Rsun, αB, Qcrit(1e-3Msun / YR), 1e51 * (100Msun / M) * 1.0)
    worst = 0.0; worstphot = 0.0; lines = String[]
    for Mdot in (1e-3Msun / YR, 1e-2Msun / YR), q in (0.5, 1.0, 2.0, 3.0, 5.0)
        Qc = Qcrit(Mdot); Q = q * Qc + Mdot / mp                     # Q − Ṁ/m_p = q Q_crit
        o = run_a(Mdot, Q; nz = 800)
        RIIan = Rin * exp(q)
        led = o.s.ledger
        pb = abs(led.Nsrc - led.Nesc - led.Nion) / led.Nsrc
        # photon budget: ionisations = recombinations + advected neutral flux
        bal = abs(led.Nion - led.Nrec - led.Nadv) / led.Nion
        push!(lines, @sprintf("   Mdot = %.0e: Q/Q_crit = %.1f (Q = %.2e): R_II/R_in = %.3f (analytic e^q = %.3f), ratio %.3f; photons %.1e, ionisations vs recombinations + advection %.1e, escaping %.1e", Mdot * YR / Msun, q + Mdot / mp / Qc, Q, o.RII / Rin, exp(q), o.RII / RIIan, pb, bal, led.Nesc / Q))
        worst = max(worst, abs(o.RII / RIIan - 1)); worstphot = max(worstphot, pb, bal)
    end
    foreach(println, lines)
    check("a1 R_II vs OI02 eq. 7", worst < 0.03, @sprintf("max |R_II/(R_in e^{(Q−Ṁ/m_p)/Q_crit}) − 1| = %.3f", worst))
    check("a2 photon budget", worstphot < 1e-6, @sprintf("max %.1e", worstphot))
    # critical Mdot: for fixed Q and M, R_II grows exponentially as Mdot falls; the HII region reaches R_out = 1e4 R_in when Mdot ≤ Mdot_crit
    Q = 1e50; Rout = 1e4 * Rin
    Mcrit_an = sqrt(8π * G * mp^2 * M * (Q - 0.0) / (αB * log(Rout / Rin))) * (1.0)           # neglecting Ṁ/m_p
    f(Mdot) = (o = run_a(Mdot, Q; nz = 600, rmaxfac = 3e4); o.RII - Rout)
    lo, hi = 0.2 * Mcrit_an, 5 * Mcrit_an
    for _ in 1:40; mid = sqrt(lo * hi); f(mid) > 0 ? (lo = mid) : (hi = mid); end
    Mc = sqrt(lo * hi)
    # analytic with the advection term included: Q − Ṁ/m_p = Q_crit(Ṁ) ln(R_out/R_in)  (solve)
    g(Md) = (Q - Md / mp) - Qcrit(Md) * log(Rout / Rin)
    l2, h2 = 1e-6Msun / YR, 1e-1Msun / YR
    for _ in 1:200; m2 = sqrt(l2 * h2); g(m2) > 0 ? (l2 = m2) : (h2 = m2); end
    Man = sqrt(l2 * h2)
    check("a3 critical Mdot", abs(Mc / Man - 1) < 0.03, @sprintf("Q = 1e50, M = 30 Msun: Ṁ_crit (R_II = 1e4 R_in) = %.4e Msun/yr, analytic (eq. 5, 6) %.4e, ratio %.4f", Mc * YR / Msun, Man * YR / Msun, Mc / Man))
end

# ====================================================================== b. uniform-velocity solution (OI02 §3): ionising radiation force = gravity
function test_b()
    M = 100Msun; Rin = 10Rsun; T = 1e4; αB = RT.ion_alpha(IA, T)[1]; Eph = 20.0; hν = hP * gr.ν[RT.ibin(gr, Eph)]      # the source sits in one bin: photon energy = bin centre
    RH = αB * hν / (4π * c * G * mp^2)
    @printf("b: α_B = %.3e, hν = %.1f eV, Haehnelt radius R_H = %.2f pc (OI02: ~26 pc for 13.6 eV)\n", αB, Eph, RH / 3.0857e18)
    nz = 800; worst_a = 0.0; worst_R = 0.0
    for RIIfac in (10.0, 30.0, 100.0)
        r = Rin .* exp.(range(0, log(3RIIfac), length = nz + 1)); rc = sqrt.(r[1:end-1] .* r[2:end])
        ρ = M ./ (4π * RH .* rc .^ 2); nH = ρ ./ mp                                       # eq. 21
        Q = αB * M^2 / (4π * RH^2 * mp^2) * (1 / Rin - 1 / (RIIfac * Rin))                  # ∫ α n² dV to R_II
        src = RT.mono_source(gr, Eph, Q)
        s = RT.ion_sweep(IA, src.L, r, fill(T, nz), nH)
        RII = rfront(s.x, r)
        V = [RT.volume(r, k) for k in 1:nz]
        arad = [sum((s.Ld[k, f] - s.Ld[k+1, f]) for f in 1:RT.nfreq(gr)) / (c * ρ[k] * V[k]) for k in 1:nz]
        gr_ = G * M ./ rc .^ 2
        inside = findall(k -> rc[k] < 0.9 * RII && rc[k] > 1.2 * Rin, 1:nz)
        dev = maximum(abs.(arad[inside] ./ gr_[inside] .- 1))
        worst_a = max(worst_a, dev); worst_R = max(worst_R, abs(RII / (RIIfac * Rin) - 1))
        @printf("   R_II,crit = %.0f R_in: Q = %.3e, sweep R_II/R_II,crit = %.4f, a_rad/g in (1.2 R_in, 0.9 R_II): %.4f … %.4f\n", RIIfac, Q, RII / (RIIfac * Rin), minimum(arad[inside] ./ gr_[inside]), maximum(arad[inside] ./ gr_[inside]))
    end
    check("b1 uniform-velocity solution", worst_a < 0.01 && worst_R < 0.03, @sprintf("a_rad = G M/r² everywhere inside the H II region to %.1e; R_II = R_II,crit to %.1e", worst_a, worst_R))
    # eq. 30: Q_crit = (hν/c)^{-1} G M²/(R_in R_H) for R_in ≪ R_II
    Qeq30 = G * M^2 / (Rin * RH) * c / hν
    r = Rin .* exp.(range(0, log(1e4), length = 1200 + 1)); rc = sqrt.(r[1:end-1] .* r[2:end]); ρ = M ./ (4π * RH .* rc .^ 2)
    Qint = αB * M^2 / (4π * RH^2 * mp^2) * (1 / Rin)
    check("b2 OI02 eq. 30", abs(Qint / Qeq30 - 1) < 1e-9, @sprintf("Q_crit(R_II → ∞) = α M²/(4π R_H² m_p² R_in) = %.4e vs (hν/c)^{-1} G M²/(R_in R_H) = %.4e", Qint, Qeq30))
end

# ====================================================================== d. multi-species sweep
function test_d()
    kpc = 3.0857e21
    # d1: Strömgren sphere (pure H, monochromatic 20 eV, uniform), steady state: R_S = (3Q/4π α n²)^{1/3}
    T = 1e4; α = RT.ion_alpha(IA, T)[1]; nH0 = 1e-3; Q = 5e48
    nz = 400; r = collect(range(0, 6.6kpc, length = nz + 1)); r[1] = 1e-3 * r[2]
    src = RT.mono_source(gr, 20.0, Q)
    s = RT.ion_sweep(IA, src.L, r, fill(T, nz), fill(nH0, nz))
    RS = cbrt(3Q / (4π * α * nH0^2)); rc = [0.5 * (r[k] + r[k+1]) for k in 1:nz]
    kf = findlast(>=(0.5), s.x); RI = rc[kf] + (s.x[kf] - 0.5) / (s.x[kf] - s.x[kf+1]) * (rc[kf+1] - rc[kf])
    # steady on-the-spot profile of the same monochromatic problem: (1 − x) Γ = α n x², Γ = Q e^{−τ} σ/(4π r²), dτ/dr = n σ (1 − x) (RK4 in ln r)
    σ1 = IA.sig[1][RT.ibin(gr, 20.0)]
    xhii(r_, τ) = (Γ = Q * exp(-τ) * σ1 / (4π * r_^2); a_ = α * nH0; y = (2a_ + Γ - sqrt((2a_ + Γ)^2 - 4a_^2)) / (2a_) ; 1 - y)       # ionised fraction
    lr = collect(range(log(1e-3 * RS), log(1.3RS), length = 20001)); τ = 0.0; xr = zeros(length(lr)); xr[1] = xhii(exp(lr[1]), 0.0)
    f_(l, τ_) = (r_ = exp(l); r_ * nH0 * σ1 * (1 - xhii(r_, τ_)))
    for i in 1:length(lr)-1
        h = lr[i+1] - lr[i]; k1 = f_(lr[i], τ); k2 = f_(lr[i] + h / 2, τ + h / 2 * k1); k3 = f_(lr[i] + h / 2, τ + h / 2 * k2); k4 = f_(lr[i+1], τ + h * k3)
        τ += h / 6 * (k1 + 2k2 + 2k3 + k4); xr[i+1] = xhii(exp(lr[i+1]), τ)
    end
    jf = findlast(>=(0.5), xr); RIode = exp(lr[jf] + (xr[jf] - 0.5) / (xr[jf] - xr[jf+1]) * (lr[jf+1] - lr[jf]))
    dev = 0.0
    for k in 1:nz
        rc[k] > 0.05 * RS && rc[k] < 1.2 * RS && (xo = xr[searchsortedlast(lr, log(rc[k]))]; dev = max(dev, abs((1 - s.x[k]) - (1 - xo)) / (1 - xo + 1e-30) * (1 - xo > 1e-3)))
    end
    check("d1 Strömgren sphere", abs(RI / RIode - 1) < 0.01 && dev < 0.05, @sprintf("steady on-the-spot profile (ODE): r(x = ½) = %.4f kpc, sweep %.4f (ratio %.4f), max rel. deviation of x_HI(r) (x_HI > 1e-3) %.1e; (R_S = (3Q/4πα n²)^{1/3} = %.4f kpc, r_I/R_S = %.3f: the 20 eV photons ionise a thicker shell)", RIode / kpc, RI / kpc, RI / RIode, dev, RS / kpc, RIode / RS))
    # d2: H + He, blackbody 1e5 K, Q(H) = 1e49, n_H = 100: nested He²⁺ ⊂ He⁺ ⊂ H⁺, photon budget, energy
    nH = 100.0; Q = 1e49; nz = 600
    src2 = RT.blackbody_photons(gr, 1e5, Q)
    RSH = cbrt(3Q / (4π * α * nH^2)); r = collect(range(0, 1.5RSH, length = nz + 1)); r[1] = 1e-4 * r[2]; rc = [0.5 * (r[k] + r[k+1]) for k in 1:nz]
    Tf = (G, x, nH_, T_) -> 1e4
    s2 = RT.ion_sweep(IAHe, src2.L, r, fill(1e4, nz), fill(nH, nz); Tfun = Tf)
    function rpos(v); kf = findlast(>=(0.5), v); kf === nothing && return 0.0; kf == nz && return r[end]; rc[kf] + (v[kf] - 0.5) / (v[kf] - v[kf+1]) * (rc[kf+1] - rc[kf]); end
    RSH = RSH; RH = rpos(s2.x); RHe1 = rpos(s2.a1 .+ s2.a2); RHe2 = rpos(s2.a2)
    led = s2.ledger
    bal = abs(led.Nion - led.Nrec) / led.Nion; pb = abs(led.Nsrc - led.Nesc - led.Nion - led.Nother) / led.Nsrc
    check("d2 H+He structure", RHe2 < RHe1 && 0.8 < RHe1 / RH < 1.1 && bal < 1e-6, @sprintf("R(H⁺) = %.3f R_S(H-only), R(He⁺) = %.3f, R(He²⁺) = %.3f R_S; ionisations vs recombinations %.1e; photon closure %.1e; escaping %.1e", RH / RSH, RHe1 / RSH, RHe2 / RSH, bal, pb, led.Nesc / Q))
    # energy: heating + ionisation energy = absorbed energy (all absorbers; no other opacity)
    Eion = IAHe.Eion; hν = hP .* gr.ν
    Ein = sum(src2.L[IAHe.f0:end]); Eout = sum(s2.Ld[end, IAHe.f0:end])
    Eion_sink = sum(s2.A[:, 1]) * Eion[1] + sum(s2.A[:, 2]) * Eion[2] + sum(s2.A[:, 3]) * Eion[3]
    Eheat = sum(s2.G[k] * RT.volume(r, k) for k in 1:nz)
    check("d3 energy", abs(Ein - Eout - Eion_sink - Eheat) / Ein < 1e-9, @sprintf("(E_in − E_out − ionisation − heating)/E_in = %.1e (mean photoelectron energy %.2f eV per ionisation)", (Ein - Eout - Eion_sink - Eheat) / Ein, Eheat / (sum(s2.A) * eV)))
    # d4: front crossing many zones: the sweep result is independent of the zone size (equilibrium, no lag)
    Rs = Float64[]
    for nz_ in (100, 200, 800)
        r_ = collect(range(0, 6.6kpc, length = nz_ + 1)); r_[1] = 1e-3 * r_[2]
        s_ = RT.ion_sweep(IA, src.L, r_, fill(T, nz_), fill(1e-3, nz_)); Vz = [RT.volume(r_, k) for k in 1:nz_]
        push!(Rs, cbrt(3 * sum(s_.x .* Vz) / (4π)))
    end
    check("d4 resolution", maximum(abs.(Rs ./ Rs[end] .- 1)) < 5e-3, @sprintf("volume-equivalent ionised radius at 100/200/800 zones: %s kpc (a front of %.0f zones per call at 800)", join([@sprintf("%.4f", x / kpc) for x in Rs], " / "), 0.5 * 800 * RS / (6.6kpc) ))
end

# ====================================================================== e. star record -> inner column -> envelope field (S7b module)
function test_e()
    include(joinpath(J, "s7b_exchange.jl")); include(joinpath(J, "s7b_feedback.jl"))
    Base.invokelatest(test_e_body)
end
function test_e_body()
    S = Main.S7b
    # e1: exchange record round trip
    d = Dict("a" => 1.5, "b" => -2e-30, "c" => 3.0e40); path = joinpath(mktempdir(), "rec.txt"); S.write_rec(path, d); r2 = S.wait_rec(path)
    check("e1 record round trip", all(r2[k] == d[k] for k in keys(d)), "write_rec/wait_rec exact")
    # e2: inner column: bolometric conservation (absorbed power re-emitted), photon balance, confinement
    M = 60Msun; R = 10Rsun; rs = 300AU; Mdot = 5e-3Msun / YR; u_s = 15e5; rho_s = Mdot / (4π * rs^2 * u_s)
    rec = Dict("Lstar" => 1e6 * S.Lsun, "R" => R, "Linf" => 1.2e6 * S.Lsun, "Tph" => 8000.0)
    src = S.star_source(gr, rec; mode = :star)
    inn = S.inner_column(IAHe, src, M, R, rs, rho_s, u_s; reemit = (gr, 8000.0))
    Lin = sum(src.L); Lout = sum(inn.Lout)
    check("e2 inner column", abs(Lout - Lin) / Lin < 1e-3 && inn.Qesc <= inn.Qstar, @sprintf("L_out/L_in = %.12f (absorbed %.2e of %.2e erg/s re-emitted); Q* = %.2e, escaping %.2e; R_II/R* = %.3f; Q/Q_crit = %.2e; ionised flux Ṁ/m_p = %.2e", Lout / Lin, inn.Labs, Lin, inn.Qstar, inn.Qesc, inn.RII / R, inn.Qstar / inn.Qcrit, Mdot / mp))
    # e3: breakout: the same star in a thin envelope (Mdot ÷ 300): the H II region reaches r_s and ionising photons enter the grid
    inn2 = S.inner_column(IAHe, src, M, R, rs, rho_s / 300, u_s; reemit = (gr, 8000.0))
    check("e3 breakout", inn2.Qesc > 0.5 * inn2.Qstar && inn2.Qstar / inn2.Qcrit > 1, @sprintf("Mdot/300: Q/Q_crit = %.1f, escaping fraction %.3f, R_II/R* = %.0f", inn2.Qstar / inn2.Qcrit, inn2.Qesc / inn2.Qstar, inn2.RII / R))
    # e4: envelope field: the sweep's beam and RadTransfer.transfer's direct beam agree (same composition); photons and momentum
    n = 60; r = collect(range(rs, 3rs, length = n + 1)); rc = 0.5 .* (r[1:end-1] .+ r[2:end]); nHz = rho_s / 300 ./ CN.μH() .* (rc ./ rs) .^ (-1.5)
    net2, ph2 = A.atomic(); xn = zeros(length(net2.free)); for (i, s) in enumerate(net2.free); s === :H2 && (xn[i] = 0.01); s === :Hp && (xn[i] = 1e-6); end
    xnet = [copy(xn) for _ in 1:n]; T = fill(1500.0, n); rho = nHz .* CN.μH()
    ab = RT.absorbers(ph2, gr)
    f = S.envfields(gr, IAHe, ab, net2, ph2, r, T, rho, nothing, xnet, inn2.Lout; opts = RT.Opts(diffuse = false), passes = 1)
    Ldir = f.rf.Ldir; Lsw = f.sw.Ld
    dev = maximum(abs.(Ldir[end, IAHe.f0:end] .- Lsw[end, IAHe.f0:end]) ./ (Lsw[1, IAHe.f0:end] .+ 1e-300))
    check("e4 sweep vs transfer", dev < 1e-3, @sprintf("transmitted ionising power: max relative difference of the two direct beams %.1e; lit zones (x > 0.5) %d of %d; |ΣaM|c = %.3e of L_abs = %.3e erg/s", dev, count(>=(0.5), f.sw.x), n, abs(sum(f.arad .* (4π / 3 .* (r[2:end] .^ 3 .- r[1:end-1] .^ 3) .* rho))) * CN.c, sum(inn2.Lout) - sum(f.rf.Ldir[end, :])))
end
which in ("a", "all") && test_a()
which in ("b", "all") && test_b()
which in ("e", "all") && test_e()
which in ("d", "all") && test_d()
if which in ("c", "all")
    include(joinpath(J, "s7b_shell.jl"))
    Base.invokelatest(test_c)
end
println("S7B: ", all(last, RESULTS) ? "ALL PASS" : "$(count(!last, RESULTS)) FAIL")
