# Population III project, M2c: 1D collapse of the (provisional) RHFC02 P100 cloud with the RHFC02 chemistry and local
# cooling (RA04's simplified transfer: thick H₂ lines × min(1, (n/8e9)^−0.45) + CIE with escape factor) — the
# configuration of RA04's own 1D test (RA04 §2.3.3, Figs. 4–5), here in TITAN.
#   TITAN_NOLIB=1 julia --project=env popiii_m2c.jl [cooling] [nshell]     cooling ∈ ra04 | ra04lines | thin
# Initial state (provisional, see POPIII_NOTES M2b): n = 1.5 polytrope shape, ρ_c/m_H = 2.6e6, cut at 1e-3 ρ_c, scaled
# to 100 M☉; T = 270 K (ρ/ρ_c)^{2/3}; at rest; f_H2 = 5e-4, x_e = 1e-10. Geometric shell masses with the central shell
# 3e-7 M☉ (RHFC02 §3.1). Regular centre, Lagrangian, self-gravity, free surface at the edge pressure.
# Log: lib/popiii/m2c_<cooling>_n<nshell>.log (central track every 0.25 dex of n_c, ON98 Fig. 3a for comparison).
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const PC = PrimordialChem; const PM = PopIIIMaterial; const Z0 = Z0Physics
const Msun = 1.98847e33

cool = isempty(ARGS) ? "ra04" : ARGS[1]; nsh = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
io = open(joinpath(TF.ROOT, "lib", "popiii", "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), "w")
pr(f, a...) = (s = Printf.format(Printf.Format(f), a...); print(s); print(io, s); flush(io))

# ---- n = 1.5 Lane–Emden θ(ξ), integrated to θ^1.5 = 1e-3
function lane_emden(n; cut = 1e-3)
    f(ξ, y) = (y[2], -y[1]^n - 2y[2] / ξ)
    ξ = 1e-4; y = (1 - ξ^2 / 6, -ξ / 3); h = 1e-4; out = [(0.0, 1.0, 0.0)]
    while max(y[1], 0.0)^n > cut
        k1 = f(ξ, y); k2 = f(ξ + h / 2, y .+ h / 2 .* k1); k3 = f(ξ + h / 2, y .+ h / 2 .* k2); k4 = f(ξ + h, y .+ h .* k3)
        y = y .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); ξ += h
        push!(out, (ξ, y[1], -ξ^2 * y[2]))          # (ξ, θ, ξ²|θ'| ∝ m)
    end
    out
end
ic = get(ENV, "IC", "p100")
dm1 = 3e-7Msun
q_of(M, n) = let lo = 1.0001, hi = 2.0
    f(q) = dm1 * (q^n - 1) / (q - 1) - M
    for _ in 1:200; m_ = (lo + hi) / 2; f(m_) > 0 ? (hi = m_) : (lo = m_); end
    (lo + hi) / 2
end
if ic == "abn02"
    # ABN02 (Abel, Bryan & Norman 2002, Science 295, 93; arXiv astro-ph/0112088) Fig. 2, profile "9 Myr later"
    # (red dotted with circles), read off the 400-dpi rendering [read-src, approximate ±0.05 dex]: particle number
    # density n(r) (panel A), T(M), H₂ mass fraction (panel C, lower curve), v_r(M) (panel E). ρ = n μ m_H, μ = 1.22.
    lr_tab = [-3.0, -1.57, -1.28, -1.08, -0.84, -0.59, -0.35, 0.06, 0.87]          # log r [pc]
    ln_tab = [5.83, 5.81, 5.40, 5.25, 4.73, 3.99, 3.47, 2.87, 0.97]                # log n [cm⁻³]
    lM_tab = [-2.0, 1.5, 2.0, 2.5, 3.0, 3.4, 3.7, 4.0]                             # log M [M☉]
    T_tab  = [252.0, 252.0, 229.0, 205.0, 164.0, 355.0, 1120.0, 1000.0]
    lf_tab = [-2.90, -2.90, -2.98, -3.11, -3.18, -3.37, -3.65, -4.0]
    v_tab  = [-0.03, -0.32, -0.49, -0.19, -0.81, -0.70, -1.73, -0.95]               # km/s
    interp(x, xs, ys) = (x <= xs[1] ? ys[1] : x >= xs[end] ? ys[end] :
        (j = searchsortedfirst(xs, x); ys[j-1] + (x - xs[j-1]) / (xs[j] - xs[j-1]) * (ys[j] - ys[j-1])))
    pc = 3.0857e18; μ_abn = 1.22
    nr(r) = 10^interp(log10(r / pc), lr_tab, ln_tab)
    # M(r) by integrating ρ(r) on a fine log grid
    rg = [pc * 10^x for x in range(-4.0, 1.5, length = 20001)]; Mg = zeros(length(rg))
    for j in 2:length(rg)
        rm = 0.5 * (rg[j] + rg[j-1]); Mg[j] = Mg[j-1] + 4π * rm^2 * nr(rm) * μ_abn * PC.mH * (rg[j] - rg[j-1])
    end
    Mg .+= 4π / 3 * rg[1]^3 * nr(rg[1]) * μ_abn * PC.mH
    M = parse(Float64, get(ENV, "MCLOUD", "1000")) * Msun
    q = q_of(M, nsh)
    mI = [0.0; cumsum([dm1 * q^(j - 1) for j in 1:nsh])]; mI[end] = M
    r_of_m(m) = m <= 0 ? 0.0 : (j = searchsortedfirst(Mg, m); j = clamp(j, 2, length(Mg));
        rg[j-1] + (m - Mg[j-1]) / (Mg[j] - Mg[j-1]) * (rg[j] - rg[j-1]))
    rI = [r_of_m(m) for m in mI]; rI[1] = 0.0
    rho = [(mI[j+1] - mI[j]) / (4π / 3 * (rI[j+1]^3 - rI[j]^3)) for j in 1:nsh]
    mz = [0.5 * (mI[j] + mI[j+1]) / Msun for j in 1:nsh]
    Tz = [interp(log10(m), lM_tab, T_tab) for m in mz]
    fz = [10^interp(log10(m), lM_tab, lf_tab) for m in mz]
    uI = [1e5 * interp(log10(max(m / Msun, 1e-2)), lM_tab, v_tab) for m in mI]; uI[1] = 0.0
    ρc = rho[1]; Tc0 = Tz[1]
    yz = [(1e-10 * 0 + 1e-7, 0.0, f / 2, 0.0) for f in fz]       # x_e ~ 1e-7 [assumed; ABN02 does not show it]
else
le = lane_emden(1.5); ξs, _, μs = le[end]
le = lane_emden(1.5); ξs, _, μs = le[end]
M = 100Msun; ρc = 2.6e6 * PC.mH; Tc0 = 270.0
a_ = cbrt(M / (4π * ρc * μs))                      # radius unit: M = 4π ρc a³ ξ²|θ'| at the cut
q = q_of(M, nsh)
mI = [0.0; cumsum([dm1 * q^(j - 1) for j in 1:nsh])]; mI[end] = M
# radius of each interface from m(ξ) (linear interpolation in the tabulated profile)
mξ = [4π * ρc * a_^3 * e[3] for e in le]; ξt = [e[1] for e in le]
function ξ_of_m(m)
    m <= 0 && return 0.0
    j = searchsortedfirst(mξ, m); j = clamp(j, 2, length(mξ))
    ξt[j-1] + (m - mξ[j-1]) / (mξ[j] - mξ[j-1]) * (ξt[j] - ξt[j-1])
end
rI = [a_ * ξ_of_m(m) for m in mI]; rI[1] = 0.0
rho = [(mI[j+1] - mI[j]) / (4π / 3 * (rI[j+1]^3 - rI[j]^3)) for j in 1:nsh]
Tz = [Tc0 * (ρ / ρc)^(2 / 3) for ρ in rho]
uI = zeros(nsh + 1); yz = nothing
end
y0 = (1e-10, 0.0, 2.5e-4, 0.0)
k5s = get(ENV, "K5", "pss")
k5f = k5s == "db" ? (T -> PC.k5_db(T)) : k5s == "ra04" ? PC.k5_ra04 : PC.k5_pss
ystore = fill(y0, 2000)
if yz !== nothing
    a0 = 5                                         # ngrs of the captured SOD state used by build!
    for j in 1:nsh; ystore[a0+j-1] = yz[j]; end
    for k in 1:a0-1; ystore[k] = yz[1]; end; for k in a0+nsh:length(ystore); ystore[k] = yz[end]; end
end
mat = PM.PrimordialGas(ystore, PC.rates_rhfc02(; k5 = k5f), cf, parse(Float64, get(ENV, "RHOEQ", "5e-11")))
PS.build!(; rI, rho, T = Tz, u = uI, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = 0.0,
          selfgrav = true, material = mat)
TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TX.PEXT[] = TF.var(:pgn)[b]
TF.setscalar!(:thet, parse(Float64, get(ENV, "THETA", "0.55")))   # 0.5 lets the centre ring (M2c)
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
tff = sqrt(3π / (32Z0.G * ρc))
TF.setscalar!(:dtime, 1e-3tff); TF.setscalar!(:tmax, 1e3tff)
pr("M2c collapse [IC %s]: cooling %s, %d shells (q = %.4f, dm1 = 3e-7 Msun), R = %.3e cm, rho_c/rho_edge = %.1f, T_edge = %.1f K, p_ext = %.3e\n",
   ic, cool, nsh, q, rI[end], rho[1] / rho[end], Tz[end], TX.PEXT[])
pr("theta = %.2f, C_dyn = %g, k5 = %s, rho_eq = %g, central t_ff0 = %.3e yr\n", TF.getscalar(:thet), TC.DT_DYN[], k5s, mat.ρeq, tff / 3.156e7)
# ON98 Fig. 3a (run A), read off: (log n_c, log T_c)
on98 = [(6.3, 2.40), (8.0, 2.60), (10.0, 2.75), (12.0, 3.05), (15.0, 3.30), (16.0, 3.35), (20.0, 3.80)]
on98T(ln) = (j = findfirst(p -> p[1] >= ln, on98); j === nothing ? NaN : j == 1 ? 10^on98[1][2] :
             10^(on98[j-1][2] + (ln - on98[j-1][1]) / (on98[j][1] - on98[j-1][1]) * (on98[j][2] - on98[j-1][2])))
pr("%7s %8s %10s %9s %9s %10s %10s %9s %9s\n", "step", "log n_c", "t [yr]", "T_c", "ON98 T_c", "f_H2,c", "u_c [km/s]", "max dT/T", "dt [s]")
nextln = Ref(log10(rho[1] / PC.mH) + 0.25)
Tstop = parse(Float64, get(ENV, "TSTOP", "6000")); nstop = parse(Float64, get(ENV, "NSTOP", "1e18"))
try
    TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k)
        dT = PM.chem!(mat)
        ρc_ = TF.var(:dn)[a]; Tc_ = TF.var(:tn)[a]; ln = log10(ρc_ / PC.mH)
        if ln >= nextln[] || k % 500 == 0
            umin = minimum(TF.var(:un)[a:b+1]) / 1e5
            pr("%7d %8.3f %10.4e %9.1f %9.1f %10.3e %10.3f %9.2e %9.2e\n", k, ln, TF.getscalar(:timen) / 3.156e7, Tc_, on98T(ln),
               2mat.y[a][3], umin, dT, TF.getscalar(:dtime))
            ln >= nextln[] && (nextln[] += 0.25)
            if get(ENV, "DTDIAG", "0") == "1"           # which variable/zone limits Δt (timstp: kx, dx per equation)
                kx = TF.var(:kx); dx = TF.var(:dx); names = ("r", "m", "rho", "u", "T")
                j = argmax(abs.(dx[1:5])); xs = TF.getscalar(:stol) / TF.getscalar(:smax)
                pr("        dt limit: %s in zone %d (rel. change %.2e, stol/smax = %.2f)\n", names[j], kx[j], dx[j], xs)
            end
        end
        (Tc_ > Tstop || ρc_ / PC.mH > nstop) && TF.setscalar!(:tmax, TF.getscalar(:timen))
    end)
catch e
    pr("stopped: %s\n", sprint(showerror, e)[1:min(end, 400)])
end
close(io)
