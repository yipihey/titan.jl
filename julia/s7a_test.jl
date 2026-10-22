# S7a acceptance tests of the frequency-dependent transfer and photochemistry (RadTransfer, s7a_hii.jl, s7a_hydro.jl). Plain Julia:
#   TITAN_NOLIB=1 julia --project=env -t 4 s7a_test.jl            (one process; ≈ 10 min incl. the hydro run; S7A_HYDRO=skip reads an existing log)
# Output: PASS/FAIL per check with the numbers; last line "S7A: ALL PASS" / "S7A: n FAIL".
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl")); include(joinpath(J, "RadTransfer.jl")); include(joinpath(J, "s7a_hii.jl"))
const CN = ChemNet; const RT = RadTransfer; const A = AtomicHHe; const yHe = CN.yHe
const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-30s %s\n", ok ? "PASS" : "FAIL", name, detail); flush(stdout)
end
const Myr = 3.15576e13; const kpc = 3.0857e21
A.ONSPOT[] = false
Bν(ν, T) = 2CN.hP * ν^3 / CN.c^2 / expm1(CN.hP * ν / (CN.kB * T))

# ====================================================================== formal solver sanity: uniform emitting sphere
let worst_c = 0.0, worst_J = 0.0, worst_Jc = 0.0
    g1 = RT.edge_grid([13.5985, 13.6015])
    for r1 in (0.0, 0.05), τR in (0.01, 1.0, 10.0)
        n = 60; R = 1.0; r = collect(range(0, R, length = n + 1)); r[1] = r1; χ = τR / R; Sf = 3.0
        op = RT.new_opacity(n, 1, 0); op.kabs .= χ; op.eta .= χ * Sf
        V = [RT.volume(r, k) for k in 1:n]; Jd = zeros(n, 1); F = zeros(n + 1, 1)
        med = RT.Medium(r, fill(1e4, n), fill(1.0, n), Dict{Symbol,Vector{Float64}}())
        RT.solve_diffuse!(Jd, F, g1, op, zeros(n, 1), V, r, med, RT.Opts(), nothing)
        Lout = 4π * R^2 * F[end, 1]; Em = sum(4π * op.eta[k, 1] * V[k] for k in 1:n); Ab = sum(4π * op.kabs[k, 1] * Jd[k, 1] * V[k] for k in 1:n)
        worst_c = max(worst_c, abs(Em - Ab - Lout) / Em)
        r1 == 0 && (worst_J = max(worst_J, abs(Jd[1, 1] - Sf * (1 - exp(-τR))) / (Sf * (1 - exp(-τR)))))
        r1 == 0 && τR == 10.0 && (worst_Jc = abs(Jd[end, 1] / Sf - 0.5) / 0.5)        # surface of a thick sphere: J = S/2
    end
    check("solver: uniform sphere", worst_c < 1e-6 && worst_J < 5e-3 && worst_Jc < 0.2, @sprintf("budget closure %.1e (conserved), J(0) vs S(1−e^−τ) %.1e, thick-surface J/S−1/2 %.2f (zone mean)", worst_c, worst_J, worst_Jc))
end

# ====================================================================== a. Iliev et al. 2006 Test 1 (monochromatic 13.6 eV, T = 1e4 K, case B, x0 = 1.2e-3)
function iliev1(; nz = 200, dt0 = 0.25)
    g = RT.edge_grid([13.5985, 13.6015]); src = RT.mono_source(g, 13.6, 5e48)
    σ = [RT.binavg(ν -> ChemRates.σ_verner(ν, 13.598, ChemRates.VERNER_HI), g, 1)]
    HII.AB_FIXED[] = 2.59e-13
    r = collect(range(0, 6.6kpc, length = nz + 1)); nH = fill(1e-3, nz)
    st = HII.State(fill(1.2e-3, nz), fill(1e4, nz))
    trec = 1 / (2.59e-13 * 1e-3); rS = cbrt(3 * 5e48 / (4π * 2.59e-13 * 1e-6))
    t = 0.0; dt = dt0 * Myr; outs = [10.0, 30.0, 100.0, 200.0, 500.0]; io = 1; res = Tuple{Float64,Float64,Float64,Float64}[]
    ph = (0.0, 0.0)
    while t < 500Myr - 1
        led = HII.sweep!(st, g, σ, src, r, nH, dt; energy = false, collisions = false); t += dt
        ph = (led.Nsrc, led.Nabs)
        if io <= length(outs) && t >= outs[io] * Myr - 1e-3 * dt
            an = rS * (1 - exp(-t / trec))^(1 / 3); push!(res, (t / Myr, HII.rfront(st, r) / an, HII.rfront_vol(st, r, nH) / an, an / kpc)); io += 1
        end
    end
    (; res, st, r, σ = σ[1], rS, trec)
end
"Steady-state Strömgren profile of the monochromatic, isothermal case (photoionisation equilibrium, on the spot): Γ(r) (1 − x) = α n x², dτ/dr = n σ (1 − x)."
function stromgren_profile(σ, rr; Ndot = 5e48, α = 2.59e-13, n = 1e-3, rmin = 1e-3 * 5.4kpc)
    # integrate in r from rmin to rr[end] (RK4 on τ), returning x_HI at the radii rr
    function y_of(r, τ)             # neutral fraction from (1 − y)² α n = Γ y,  Γ = Ndot e^{-τ} σ /(4π r²)
        Γ = Ndot * exp(-τ) * σ / (4π * r^2); a = α * n; (2a + Γ - sqrt((2a + Γ)^2 - 4a^2)) / (2a)
    end
    out = zeros(length(rr)); τ = 0.0; r = rmin; j = 1
    f(r, τ) = n * σ * y_of(r, τ)
    while j <= length(rr)
        h = min(rr[j] - r, 1e-3 * 5.4kpc); h <= 0 && (out[j] = y_of(r, τ); j += 1; continue)
        k1 = f(r, τ); k2 = f(r + h / 2, τ + h / 2 * k1); k3 = f(r + h / 2, τ + h / 2 * k2); k4 = f(r + h, τ + h * k3)
        τ += h / 6 * (k1 + 2k2 + 2k3 + k4); r += h
    end
    out
end
let
    t0 = time(); R1 = iliev1(); wall = time() - t0
    rmax = maximum(abs(x[2] - 1) for x in R1.res); rvmax = maximum(abs(x[3] - 1) for x in R1.res)
    check("a1 Iliev T1 r_I(t)", rmax < 0.05, @sprintf("t = %s Myr: r_I(x = 0.5)/analytic = %s; volume-equivalent %s; wall %.0f s", join([@sprintf("%.0f", x[1]) for x in R1.res], "/"),
          join([@sprintf("%.3f", x[2]) for x in R1.res], "/"), join([@sprintf("%.3f", x[3]) for x in R1.res], "/"), wall))
    r = R1.r; nz = length(R1.st.x); rc = [0.5 * (r[k] + r[k+1]) for k in 1:nz]
    sel = [k for k in 1:nz if 0.25 < rc[k] / R1.rS < 0.85]
    yan = stromgren_profile(R1.σ, rc[sel]); ysim = 1 .- R1.st.x[sel]
    rat = ysim ./ yan; kk = [findfirst(k -> rc[k] / R1.rS >= q, sel) for q in (0.3, 0.5, 0.7, 0.8)]
    check("a2 Iliev T1 x_HI(r) 500 Myr", maximum(abs.(log.(rat))) < 0.5, @sprintf("x_HI(sim)/x_HI(steady Strömgren) at r/rS = 0.3/0.5/0.7/0.8: %s (sim %s; analytic %s); range of ratio %.2f–%.2f",
          join([@sprintf("%.2f", rat[k]) for k in kk], "/"), join([@sprintf("%.2e", ysim[k]) for k in kk], "/"), join([@sprintf("%.2e", yan[k]) for k in kk], "/"), minimum(rat), maximum(rat)))
end

# ====================================================================== b. Iliev et al. 2006 Test 2 (10^5 K blackbody, T evolves)
function iliev2(; nz = 200, dt0 = 0.25, dl = 0.05, tend = 500.0)
    E = 13.598 .* exp.(range(0, log(100 / 13.598), step = dl)); E[1] = 13.5985
    g = RT.edge_grid(E); src = RT.blackbody_photons(g, 1e5, 5e48)
    σ = [RT.binavg(ν -> ChemRates.σ_verner(ν, 13.598, ChemRates.VERNER_HI), g, i) for i in 1:RT.nfreq(g)]
    HII.AB_FIXED[] = 0.0
    HII.init_ek([(CN.χH * (1 - 1 / n^2), 2.0n^2, CN.χH / n^2, ν -> A.σ_hyd(n, 1, ν)) for n in 2:15])
    r = collect(range(0, 6.6kpc, length = nz + 1)); nH = fill(1e-3, nz)
    st = HII.State(fill(1e-8, nz), fill(100.0, nz))
    trec = 1 / (2.59e-13 * 1e-3); rS = cbrt(3 * 5e48 / (4π * 2.59e-13 * 1e-6))
    t = 0.0; dt = dt0 * Myr; outs = [10.0, 30.0, 100.0, 200.0, 500.0]; io = 1; res = []; led = nothing
    while t < tend * Myr - 1
        led = HII.sweep!(st, g, σ, src, r, nH, dt; energy = true, collisions = true); t += dt
        if io <= length(outs) && t >= outs[io] * Myr - 1e-3 * dt
            rI = HII.rfront(st, r); an = rS * (1 - exp(-t / trec))^(1 / 3); kf = findlast(k -> st.x[k] > 0.5, 1:nz)
            push!(res, (t = t / Myr, ratio = rI / an, Tin = st.T[max(1, round(Int, 0.3kf))], Tout = st.T[round(Int, 0.9kf)], xhi = 1 - st.x[round(Int, 0.5kf)],
                        Tmax = maximum(st.T[1:kf]), Tpre = st.T[min(nz, kf + 8)], xpre = st.x[min(nz, kf + 8)])); io += 1
        end
    end
    (; res, st, r, rS, trec, led, g)
end
let
    t0 = time(); R2 = iliev2(); wall = time() - t0
    rt = [x.ratio for x in R2.res]
    # Iliev et al. 2006 Fig. 15 (top): r_num/r_analytic(1e4 K) of the participating codes: 0.93–0.97 at t → 0, rising to 1.0–1.10 at 4 t_rec (read off the plot, ±0.02)
    check("b1 Iliev T2 r_I(t)", all(0.88 .< rt .< 1.13) && rt[end] > 1.0 && rt[end] > rt[1],
          @sprintf("t = %s Myr: r_I/r_an(1e4 K) = %s (codes 0.93–1.10, paper Fig. 15); wall %.0f s", join([@sprintf("%.0f", x.t) for x in R2.res], "/"), join([@sprintf("%.3f", x) for x in rt], "/"), wall))
    Tin = [x.Tin for x in R2.res]; Tmx = [x.Tmax for x in R2.res]
    check("b2 Iliev T2 temperature", all(1.0e4 .< Tin[2:end] .< 4.0e4) && all(Tmx .< 4.5e4),
          @sprintf("T at 0.3 r_I = %s K, at 0.9 r_I = %s K; max T in the HII region %s K (paper: 1–4e4 K inside); neutral fraction at 0.5 r_I %s; pre-front gas (8 zones ahead) T %s K, x %s",
                   join([@sprintf("%.0f", x) for x in Tin], "/"), join([@sprintf("%.0f", x.Tout) for x in R2.res], "/"), join([@sprintf("%.0f", x) for x in Tmx], "/"),
                   join([@sprintf("%.1e", x.xhi) for x in R2.res], "/"), join([@sprintf("%.0f", x.Tpre) for x in R2.res], "/"), join([@sprintf("%.1e", x.xpre) for x in R2.res], "/")))
    L = R2.led
    check("b3 sweep photon budget", abs(L.Nsrc - L.Nabs - L.Nesc) / L.Nsrc < 1e-12, @sprintf("(emitted − absorbed − escaping)/emitted = %.1e at 500 Myr (last sweep)", (L.Nsrc - L.Nabs - L.Nesc) / L.Nsrc))
end

# ====================================================================== the general transfer in a static HII region: d (conservation), cost g, h (ChemNet consistency), a3 (diffuse = on the spot)
net, photos = A.atomic()
photosG = [p for p in photos if occursin("ground", p.name)]            # ground-level photo pairs only (a clean ionising-photon ledger)
function hii_medium(n; nHc = 100.0, Teff = 1e5, Ndot = 5e48)
    rS = cbrt(3 * Ndot / (4π * 2.6e-13 * nHc^2))
    r = collect(10 .^ range(log10(1e-3 * rS), log10(3rS), length = n + 1)); r[1] = 0.0
    rc = [0.5 * (r[k] + r[k+1]) for k in 1:n]; ion = rc .< rS
    X = Dict{Symbol,Vector{Float64}}()
    X[:Hp] = [i ? 1.0 - 1e-4 : 1e-4 for i in ion]; X[:H] = 1 .- X[:Hp]
    X[:He] = [i ? 0.0 : yHe for i in ion]; X[:Hep] = [i ? 0.2yHe : 0.0 for i in ion]; X[:Hepp] = [i ? 0.8yHe : 0.0 for i in ion]
    X[:e] = X[:Hp] .+ X[:Hep] .+ 2 .* X[:Hepp]; X[:H2] = zeros(n); X[:Hm] = zeros(n)
    T = [i ? 1.5e4 : 100.0 for i in ion]
    RT.Medium(r, T, fill(nHc, n), X), rS
end
let
    g = RT.default_grid(); med, rS = hii_medium(120); n = length(med)
    src = RT.blackbody_photons(g, 1e5, 5e48); rho = fill(100.0 * CN.μH(), n); V = [RT.volume(med.r, k) for k in 1:n]
    # --- d1: direct beam, ground-level absorbers only: ΣΓ n V = emitted − escaping (independent of the J ledger), energy ∫(heating + ionisation sink) = absorbed
    abG = RT.absorbers(photosG, g); oD = RT.Opts(diffuse = false, co = false, emit = false)
    rf = RT.transfer(g, abG, src, med, oD; rho)
    Nabs_G = sum(sum(rf.Gamma[ia][k] * med.nH[k] * RT.xs(med, abG[ia].R, k) * V[k] for k in 1:n) for ia in eachindex(abG))
    L = rf.ledger; Nin = L.Nsrc; Nout = L.Ndir_out
    e1 = abs(Nabs_G - (Nin - Nout)) / Nin
    Eabs_flux = sum(src.L[f] - rf.Ldir[end, f] for f in 1:RT.nfreq(g) if g.edges[f] >= CN.χH / CN.hP * (1 - 1e-9))
    Eheat = sum(rf.heat .* V) + sum(rf.ionrate .* V); e2 = abs(Eheat - Eabs_flux) / Eabs_flux
    check("d1 photons (direct)", e1 < 1e-3, @sprintf("(ΣΓ_p n_p V − (Ṅ_emitted − Ṅ_escaping))/Ṅ = %.1e; emitted %.4e, escaping %.3e, absorbed %.4e s⁻¹ (tolerance 1e-3)", e1, Nin, Nout, Nabs_G))
    check("d2 energy (direct)", e2 < 1e-6, @sprintf("(∫heating + ionisation sink − absorbed energy)/absorbed = %.1e (tolerance 1e-6); heating %.4e + ionisation %.4e erg/s", e2, sum(rf.heat .* V), sum(rf.ionrate .* V)))
    # --- d3: direct + diffuse, all absorbers and continua, scattering iterated
    ab = RT.absorbers(photos, g)
    rf2 = RT.transfer(g, ab, src, med, RT.Opts(); rho); L2 = rf2.ledger
    e3 = (L2.Nsrc + L2.Nem - L2.Ndir_out - L2.Ndif_out - L2.Nabs) / L2.Nsrc
    e4 = abs(L2.Eheat - L2.Eabs) / L2.Eabs
    check("d3 photons (+ diffuse)", abs(e3) < 1e-3 && e4 < 1e-6, @sprintf("(Ṅ_src + Ṅ_gas emission − escaping − absorbed)/Ṅ_src = %.1e; gas emission %.3e s⁻¹; scattering sweeps %d; energy closure %.1e", e3, L2.Nem, rf2.iters, e4))
    # --- i1: ChemNet.step with the per-zone env and the atomic cooling (photo-heating through env.J): the energy gained over a short step = the transfer's absorbed energy
    let k = 105, ze1 = RT.ZoneEnv(g, photos)
        RT.zonehook(ze1, rf2)(k, nothing)
        T0 = med.T[k]; nHk = med.nH[k]; ρk = nHk * CN.μH()
        x0 = [RT.xs(med, s, k) for s in net.free]
        ΓH = sum(rf2.Gamma[ia][k] for (ia, a) in enumerate(ab) if a.R === :H); Γ0 = sum(rf2.Gamma[ia][k] * RT.xs(med, ab[ia].R, k) for ia in eachindex(ab))
        dt = 1e-4 / Γ0
        cool1 = A.cooling(photos, ze1)
        u0 = CN.eos(ρk, T0, x0, net)[1]
        x1, T1 = CN.step(x0, ρk, u0, T0, dt, net; env = ze1, photos, cooling = cool1)
        u1 = CN.eos(ρk, T1, x1, net)[1]
        Eexp = (rf2.heat[k] - rf2.heat_other[k] + rf2.ionrate[k]) * dt / ρk          # absorbed energy per gram of the network's absorbers
        ip = findfirst(==(:Hp), net.free)
        check("i1 ChemNet.step with ZoneEnv", abs((u1 - u0) / Eexp - 1) < 0.05 && x1[ip] > x0[ip],
              @sprintf("zone %d (T = %.0f K, x_H+ = %.1e), dt = 1e-4/Γ = %.2e s: Δu·ρ/dt = %.3e vs transfer absorbed energy (heating + ionisation) %.3e erg cm⁻³ s⁻¹ (ratio %.3f; the difference is recombination/line cooling); Δx_H+ = %.3e (Γ_H x_H dt = %.3e)",
                      k, T0, x0[ip], dt, (u1 - u0) * ρk / dt, Eexp * ρk / dt, (u1 - u0) / Eexp, x1[ip] - x0[ip], ΓH * med.X[:H][k] * dt))
    end
    # --- radiation force: momentum conservation of the direct beam; compare (1/ρc)∫κF with the divergence form in a thin medium
    Mom = sum(rf.arad .* rho .* V); Mabs = sum(sum(rf.Ldir[k, :] .- rf.Ldir[k+1, :]) for k in 1:n) / CN.c
    medT, rST = hii_medium(60; nHc = 1e-2); rhoT = fill(1e-2 * CN.μH(), 60); VT = [RT.volume(medT.r, k) for k in 1:60]
    rfT = RT.transfer(g, ab, src, medT, RT.Opts(diffuse = false); rho = rhoT, op = nothing)
    opT = RT.new_opacity(60, RT.nfreq(g), length(ab)); RT.fill_opacity!(opT, g, ab, medT, RT.Opts(diffuse = false))
    aF = [sum((opT.kabs[k, f] + opT.ssc[k, f]) * 0.5 * (rfT.Ldir[k, f] + rfT.Ldir[k+1, f]) / (4π * ((0.5 * (medT.r[k]^3 + medT.r[k+1]^3))^(1 / 3))^2) / CN.c / rhoT[k] for f in 1:RT.nfreq(g)) for k in 1:60]
    kk = [k for k in 5:60 if medT.r[k+1] < 0.7rST]; rdev = maximum(abs.(aF[kk] ./ rfT.arad[kk] .- 1))
    check("d4 radiation force", abs(Mom / Mabs - 1) < 1e-10 && rdev < 0.2, @sprintf("Σ a ρ V = %.4e dyn vs (L_in − L_out)/c = %.4e (direct beam, ratio − 1 = %.1e); (1/ρc)Σχ L_mid/(4π r²) vs the divergence form in a thin n_H = 1e-2 medium: max dev %.1e over the %d thin ionised-interior zones (they differ where a zone is thick: the arithmetic mean of L); a(r) at 1e17 / 1e18 / 1e19 cm = %.2e / %.2e / %.2e cm s⁻²",
          Mom, Mabs, Mom / Mabs - 1, rdev, length(kk), rf2.arad[findfirst(>(1e17), med.r) - 1], rf2.arad[findfirst(>(1e18), med.r) - 1], rf2.arad[findfirst(>(1e19), med.r) - 1]))
    # --- g: cost
    nz = 200; med2, _ = hii_medium(nz); rho2 = fill(100.0 * CN.μH(), nz)
    RT.transfer(g, ab, src, med2, RT.Opts(); rho = rho2)
    tf = @elapsed RT.transfer(g, ab, src, med2, RT.Opts(); rho = rho2)
    td = @elapsed RT.transfer(g, ab, src, med2, RT.Opts(diffuse = false); rho = rho2)
    # per-zone ChemNet env: set_zone! (J at the quadrature nodes, bin integrals) + the Milne rates and energies of all photo pairs at 6 temperatures (one Newton step's worth), per zone
    zee = RT.ZoneEnv(g, photos); Jz = rf2.J[40, :]
    RT.set_zone!(zee, Jz); for p in photos; CN.photorates_t(p, 1.5e4, zee); CN.photoenergy_t(p, 1.5e4, zee); end
    tz = @elapsed for rep in 1:20
        RT.set_zone!(zee, Jz)
        for Tt in (1.5e4, 1.52e4, 1.55e4, 1.5e4 * 1.03, 1.5e4 * 1.01, 1.5e4 * 0.99), p in photos; CN.photorates_t(p, Tt, zee); CN.photoenergy_t(p, Tt, zee); end
    end
    check("g cost", true, @sprintf("one transfer call, %d zones × %d bins × %d absorbers, %d threads: full (opacity + direct + diffuse + rates) %.2f s; direct only (with opacity, co) %.2f s; per-zone ChemNet env (set_zone! + rates/energies of %d photo pairs at 6 temperatures) %.1f ms per zone = %.2f s per 200 zones", nz, RT.nfreq(g), length(ab), Threads.nthreads(), tf, td, length(photos), 1e3 * tz / 20, 200 * tz / 20))
    # --- a3: with the diffuse field the recombination continuum is absorbed on the spot (case B emerges): Γ_diff n_HI ≈ α_1 n_p n_e in the ionised interior
    gb = RT.default_grid(); abH = RT.absorbers([photosG[findfirst(p -> p.R === :H, photosG)]], gb)
    n3 = 150; med3, rS3 = hii_medium(n3); Xs = med3.X
    for k in 1:n3; Xs[:Hp][k] = Xs[:Hp][k] > 0.5 ? 1 - 1e-4 : Xs[:Hp][k]; end
    srcH = RT.blackbody_photons(gb, 1e5, 5e48)
    rfD = RT.transfer(gb, abH, srcH, med3, RT.Opts(co = false, emit = true, stim = true))
    rf0 = RT.transfer(gb, abH, srcH, med3, RT.Opts(co = false, emit = true, stim = true, diffuse = false))
    ph = photosG[findfirst(p -> p.R === :H, photosG)]
    rat = Float64[]; Sdiff = 0.0; Semit = 0.0; V3 = [RT.volume(med3.r, k) for k in 1:n3]
    ph = photosG[findfirst(p -> p.R === :H, photosG)]
    for k in 1:n3
        T = med3.T[k]; _, α1 = CN.photorates(ph, T, (Trad = 0.0, J = ν -> 0.0, static = false))
        nH = med3.nH[k]; nHI = nH * Xs[:H][k]; em = α1 * nH^2 * Xs[:Hp][k] * Xs[:e][k]
        ab_ = (rfD.Gamma[1][k] - rf0.Gamma[1][k]) * nHI
        Sdiff += ab_ * V3[k]; Semit += em * V3[k]
        rc = 0.5 * (med3.r[k] + med3.r[k+1]); (0.3 < rc / rS3 < 0.8) && push!(rat, ab_ / em)
    end
    check("a3 diffuse ⇒ on the spot", abs(Sdiff / Semit - 1) < 0.05, @sprintf("Σ(diffuse photons absorbed by H) / Σ(ground-level recombinations α_1 n_p n_e V) = %.4f (case B emerges when the recombination continuum is transported); locally at 0.3–0.8 r_S the ratio is %.2f–%.2f (the interior is thin: mean free path 1/(n_HI σ) = %.1e cm > r_S = %.1e cm, so those photons are absorbed at the front)",
          Sdiff / Semit, minimum(rat), maximum(rat), 1 / (100.0 * 1e-4 * 6.3e-18), rS3))
    # --- h: ChemNet consistency of the per-zone env
    ze = RT.ZoneEnv(g, photos); eR = 0.0; eA = 0.0; emaxI = 0.0; emaxC = 0.0; nzk = 0
    zex = RT.ZoneEnv(g, photos; exact = true)
    for k in (30, 60, 90)
        T = med.T[k]; Jb = rf2.J[k, :]; RT.set_zone!(ze, Jb); RT.set_zone!(zex, Jb)
        envh = (Trad = 0.0, J = ze.J, static = false)
        for (ia, p) in enumerate(photos)
            R1, A1 = CN.photorates(p, T, envh); R2, A2 = CN.photorates_t(p, T, ze); R3, A3 = CN.photorates_t(p, T, zex)
            Rg = rf2.Gamma[ia][k]
            (R1 > 1e-30 && R3 > 1e-30) || continue
            eR = max(eR, abs(R1 / R3 - 1)); emaxI = max(emaxI, abs(R2 / R3 - 1)); emaxC = max(emaxC, abs(Rg / R3 - 1))
            A1 > 1e-40 && (eA = max(eA, abs(A3 / A1 - 1)); emaxI = max(emaxI, abs(A2 / A3 - 1)))
            nzk += 1
        end
    end
    check("h1 ChemNet rates vs transfer", emaxC < 1e-10 && emaxI < 1e-3 && eA < 3e-2 && eR < 0.15,
          @sprintf("env R (bins of the transfer) vs the transfer's Γ: %.1e (photon conserving by construction); ChemNet.photorates (its own node quadrature of the binned J) vs env: R %.1e (aliasing of ChemNet's nodes against the J bins, 1 node per bin), Milne α %.1e (T ≥ 1.5e4 K; at 100 K both quadratures under-resolve the threshold); T-tabulated vs exact %.1e; %d photo pairs × zones", emaxC, eR, eA, emaxI, nzk))
    # heating consistency: AtomicHHe cooling (rad term, env.J = J_tot) vs the transfer's photo-heating: no double counting (the Photo pairs carry −n_X RE in `cooling`, the transfer reports the same absorbed energy)
    worst = 0.0
    for k in (30, 60, 90)
        T = med.T[k]; RT.set_zone!(ze, rf2.J[k, :]); nH = med.nH[k]
        d = Dict(s => v[k] for (s, v) in med.X)
        envJ = ze; env0 = (Trad = 0.0, J = ν -> 0.0, static = false)
        pJ = A.parts(med.nH[k] * CN.μH(), T, d; photos, env = envJ); p0 = A.parts(med.nH[k] * CN.μH(), T, d; photos, env = env0)
        # rad(J) − rad(0) = −Σ n_X RE + stimulated-recombination energy changes; compare −Σ n_X RE with the transfer's absorbed energy of the same absorbers
        ΔRE = 0.0
        for p in photos; RE, _, _ = CN.photoenergy_t(p, T, envJ); ΔRE += nH * get(d, p.R, 0.0) * RE; end
        Eabs_k = rf2.heat[k] - rf2.heat_other[k] + rf2.ionrate[k]
        worst = max(worst, abs(ΔRE / Eabs_k - 1))
    end
    check("h2 heating, no double count", worst < 3e-2, @sprintf("Σ_p n_p RE_p (what `AtomicHHe.cooling` subtracts via env.J) / (transfer: heating + ionisation sink of the same absorbers) − 1: max %.1e at 3 zones; the transfer's own heating is NOT added when the atomic cooling is used with this env", worst))
end

# ====================================================================== e. Lyman–Werner self-shielding
let
    g = RT.default_grid(); iLW = RT.ibin(g, 12.87); Jν = 1e-21
    # independent implementations of DB96 eq. 37 and WGHB11 eq. 10
    fdb(N, T) = (x = N / 5e14; b5 = sqrt(2 * 1.380649e-16 * T / (2 * 1.6735575e-24)) / 1e5; 0.965 / (1 + x / b5)^2 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x)))
    fwg(N, T) = (x = N / 5e14; b5 = sqrt(2 * 1.380649e-16 * T / (2 * 1.6735575e-24)) / 1e5; 0.965 / (1 + x / b5)^1.1 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x)))
    worst = 0.0; line = String[]; Ns = 10.0 .^ (14:21)
    for N in Ns, T in (300.0, 1e3, 5e3)
        # one zone with column N from the source to its centre (direct field), J(12.87 eV) = 1e-21 in the LW bin
        dr = 1e15; nH2 = 2N / dr; r = [1e17, 1e17 + dr]
        med = RT.Medium(r, [T], [nH2 * 1e3], Dict(:H2 => [1e-3]))
        Jdir = zeros(1, RT.nfreq(g)); Jdir[1, iLW] = Jν
        k, JLW, fs = RT.lw_rates(g, med, Jdir, zeros(1, RT.nfreq(g)), nothing)
        kex = 1.38e9 * Jν * fwg(N, T)
        worst = max(worst, abs(k[1] / kex - 1))
        T == 1e3 && push!(line, @sprintf("%.0e: %.2e (f %.3f; DB96 %.3f)", N, k[1], fs[1], fdb(N, T)))
    end
    mono = all(diff([RT.fsh_lw(N, 1e3) for N in Ns]) .< 0) && RT.fsh_lw(1e14, 1e3) > 0.9 && RT.fsh_lw(1e21, 1e3) < 1e-4
    check("e1 k_LW(N_H2)", worst < 1e-12 && mono, @sprintf("k_LW [s⁻¹] = 1.38e9 J f_sh(N, T) for J = 1e-21 (T = 1e3 K, WGHB11 α = 1.1): %s; max |k/formula − 1| over N = 1e14…1e21 × T = 300, 1e3, 5e3 K: %.1e", join(line, "  "), worst))
    rdb = [fwg(N, 1e3) / fdb(N, 1e3) for N in Ns]
    check("e2 WGHB11 vs DB96", true, @sprintf("f_sh(WGHB11)/f_sh(DB96) at T = 1e3 K, N = 1e14…1e21: %s (WGHB11: α = 1.1 instead of 2, ~15%% agreement with 3D halo simulations)", join([@sprintf("%.2f", x) for x in rdb], "/")))
    # the reaction in the network (opt-in): H₂ + γ → 2H with rate env.kLW
    netm, phm = ChemRates.modern(); netL = RT.with_lw(netm)
    x = [0.0, 0.0, 1e-3, 0.0]; nH = 1e3; T = 100.0
    envL = (Trad = 0.0, J = ν -> 0.0, static = true, kLW = 3e-10); env0 = (Trad = 0.0, J = ν -> 0.0, static = true)
    f1 = CN.dxdt(x, nH, T, netL, envL, []); f0 = CN.dxdt(x, nH, T, netm, env0, [])
    i2 = findfirst(==(:H2), netL.free)
    check("e3 network LW reaction", abs((f1[i2] - f0[i2]) / (-3e-10 * 1e-3) - 1) < 1e-10 && length(netL.reactions) == length(netm.reactions) + 1,
          @sprintf("dx_H2/dt with the LW reaction − without = %.4e vs −k x_H2 = %.4e (k = 3e-10 s⁻¹); the reaction is added only by RT.with_lw (default networks unchanged)", f1[i2] - f0[i2], -3e-10 * 1e-3))
end

# ====================================================================== f. detailed balance with the env of this module (J = B(T))
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
netB, phB = A.atomic(li = true, jb_inverse = true)
function fluxcheck(net, ph, d, nH, T, env)
    rx = net.reactions; flux = [(rate = r.k(T, env) * nH^(length(r.R) - 1); for s in r.R; rate *= d[s]; end; rate) for r in rx]
    worst = 0.0; idx = Dict((sort(r.R), sort(r.P)) => i for (i, r) in enumerate(rx))
    for (i, r) in enumerate(rx)
        j = get(idx, (sort(r.P), sort(r.R)), 0); (j == 0 || j < i) && continue
        f, b = flux[i], flux[j]; m = max(f, b); m < 1e-250 && continue
        worst = max(worst, abs(f - b) / m)
    end
    for p in ph
        R, α = CN.photorates_t(p, T, env); α *= p.esc(T, env)
        fw = R * d[p.R]; bw = α * nH^(length(p.P) - 1); for s in p.P; bw *= d[s]; end
        m = max(fw, bw); m < 1e-250 && continue
        worst = max(worst, abs(fw - bw) / m)
    end
    worst
end
let
    gT = RT.edge_grid(exp.(range(log(0.05), log(3e4), step = 0.05))); zeB = RT.ZoneEnv(gT, phB; exact = true); zeI = RT.ZoneEnv(gT, phB)
    w1 = 0.0; w2 = 0.0; ns = 0
    for T in (1e3, 1e4, 1e5, 1e6), nH in (1e0, 1e8, 1e16)
        d = saha(nH, T); RT.set_thermal!(zeB, T); w1 = max(w1, fluxcheck(netB, phB, d, nH, T, zeB))
        RT.set_thermal!(zeI, T); w2 = max(w2, fluxcheck(netB, phB, d, nH, T, zeI)); ns += 1
    end
    check("f1 detailed balance (env)", w1 < 1e-9, @sprintf("J = B(T) through ZoneEnv (exact T): max |forward − reverse|/flux over all reactions and photo pairs at the Saha composition, %d states T = 1e3–1e6 K × n_H = 1e0–1e16: %.1e (S6c a1: 6.5e-14 with the direct quadrature)", ns, w1))
    check("f2 detailed balance (tab)", w2 < 5e-3, @sprintf("same with the T-tabulated rates (0.02 dex, Catmull–Rom; the production mode): %.1e", w2))
    # relaxation to Saha with the env (pseudo-time stepping of ChemNet.steady; as S6c a2: every free species × 0.97 / 1.03, caps for H, He, Li)
    wr = 0.0; nfail = 0; nrel = 0
    for (T, nH) in ((3e3, 1e4), (1e4, 1e12), (3e4, 1e4), (1e5, 1e12), (3e5, 1e4))
        ex = saha(nH, T); RT.set_thermal!(zeB, T)
        x0 = [ex[s] * (iseven(i) ? 0.97 : 1.03) for (i, s) in enumerate(netB.free)]
        for (i, s) in enumerate(netB.free)
            s === :Hp && (x0[i] = min(x0[i], 1 - 1e-9 - 2 * ex[:H2] - ex[:H2p] - ex[:Hm]))
            s === :Hepp && (x0[i] = min(x0[i], yHe * 0.999999))
            s === :Hep && (x0[i] = min(x0[i], yHe * 0.999999 - x0[findfirst(==(:Hepp), netB.free)]))
            s === :Lip && (x0[i] = min(x0[i], CN.yLi * 0.999999))
        end
        x = try CN.steady(x0, nH * CN.μH(), T, netB; env = zeB, photos = phB, xfloor = 1e-14, tol = 1e-11, xclamp = 1e-14) catch e; (nfail += 1; continue) end
        d = CN.full(x, netB); nrel += 1
        for s in (:H, :Hp, :Hm, :H2, :H2p, :He, :Hep, :Hepp, :Li, :Lip)
            ex[s] < 1e-6 && continue; wr = max(wr, abs(d[s] / ex[s] - 1))
        end
    end
    check("f3 relaxation to Saha", wr < 5e-6 && nfail == 0, @sprintf("pseudo-transient relaxation from ±3%% perturbed states with the module's env (J = B analytic): max |x/x_Saha − 1| = %.1e over %d states (%d failed; S6c a2: 9.8e-7)", wr, nrel, nfail))
end

# ====================================================================== c. Iliev et al. 2009 Test 5 (hydrodynamic H II region expansion), TITAN hydro + this transfer/chemistry
let
    log_ = joinpath(J, "..", "lib", "popiii", "s7a_hydro_n300.csv")
    if get(ENV, "S7A_HYDRO", "run") == "run" || !isfile(log_)
        println("running s7a_hydro.jl (NZ = 300, 500 Myr) …"); flush(stdout)
        withenv("TITAN_NOLIB" => "1", "NZ" => "300", "TEND" => "500", "TAG" => "n300") do
            run(pipeline(`$(Base.julia_cmd()) --project=$(joinpath(J, "env")) -t 1 $(joinpath(J, "s7a_hydro.jl"))`; stdout = joinpath(J, "..", "lib", "popiii", "s7a_hydro_n300.log"), stderr = stdout))
        end
    end
    rows = [parse.(Float64, split(l, ",")) for l in readlines(log_)[2:end]]
    t = [r[1] for r in rows]; rI = [r[2] for r in rows]; rr = [r[3] for r in rows]; Th = [r[4] for r in rows]
    rS0 = 5.3931; kpc_ = 1.0
    spitzer(t_, Tg) = (cs = sqrt(2 * CN.kB * Tg / CN.mH); rS0 * (1 + 7 * cs * t_ * Myr / (4 * rS0 * kpc))^(4 / 7))
    sp = [spitzer(t[i], Th[i]) for i in eachindex(t)]
    i100 = argmin(abs.(t .- 100)); i200 = argmin(abs.(t .- 200)); i300 = argmin(abs.(t .- 300)); i500 = argmin(abs.(t .- 500))
    # Iliev et al. 2009 Fig. 18 (top), codes: r_I/r_S0 ≈ 0.87 (100 Myr), 1.12 (200), 1.32 (300), 1.57 (500, the codes span 1.5–1.65), read off the plot (±0.05)
    ref = Dict(100 => 0.87, 200 => 1.12, 300 => 1.32, 500 => 1.57)
    dev = [rr[i] / ref[Int(round(t[i]))] - 1 for i in (i100, i200, i300, i500)]
    check("c1 Iliev T5 r_I(t)", all(abs.(dev) .< 0.12), @sprintf("r_I/r_S0 at 100/200/300/500 Myr = %s (Iliev 2009 codes: 0.87/1.12/1.32/1.57): deviation %s; r_I = %s kpc", join([@sprintf("%.3f", rr[i]) for i in (i100, i200, i300, i500)], "/"),
          join([@sprintf("%+.1f%%", 100d) for d in dev], "/"), join([@sprintf("%.2f", rI[i]) for i in (i100, i200, i300, i500)], "/")))
    spr = [rI[i] / sp[i] for i in (i100, i200, i300, i500)]
    check("c2 Spitzer law", all(0.7 .< spr[2:end] .< 1.15) && issorted(spr) && issorted([rI[i100], rI[i200], rI[i300], rI[i500]]),
          @sprintf("r_I / [r_S0 (1 + 7 c_s t/(4 r_S0))^{4/7}] at 100/200/300/500 Myr = %s with c_s = (2kT/m_H)^½ from the mass-weighted T of the ionised gas (%s K): the Spitzer law is the r_I ≫ r_S0 asymptote and over-predicts here (r_I/r_S0 ≤ 1.7); the ratio rises towards 1 with time (criterion for t ≥ 200 Myr: 0.7–1.15, increasing)",
                   join([@sprintf("%.2f", x) for x in spr], "/"), join([@sprintf("%.0f", Th[i]) for i in (i100, i200, i300, i500)], "/")))
    # local slope at late times d ln r_I/d ln t between 300 and 500 Myr vs 4/7 (asymptote of the Spitzer solution)
    sl = log(rI[i500] / rI[i300]) / log(t[i500] / t[i300])
    check("c3 late-time slope", 0.2 < sl < 0.7, @sprintf("d ln r_I/d ln t (300 → 500 Myr) = %.2f (Spitzer asymptote 4/7 = 0.57; the Iliev et al. 2009 codes also lie below it in this box)", sl))
end

nf = count(!last, RESULTS)
println(nf == 0 ? "S7A: ALL PASS ($(length(RESULTS)) checks)" : "S7A: $nf FAIL of $(length(RESULTS))")
exit(nf == 0 ? 0 : 1)
