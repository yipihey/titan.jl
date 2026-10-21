# S6a-fix: energy accounting of the prescribed (ALE) mass grid. Static, adiabatic, non-burning, gas-only γ = 5/3 polytrope in exact discrete
# hydrostatic equilibrium (radii marched outward so that TITAN's momentum rows vanish at rest), regular centre, free surface at constant P_ext.
#   TITAN_NOLIB=1 julia --project=env s6fix_test.jl [relabel|accrete]
#   relabel: the mass grid drifts through the matter (m_k → m_k(1 + ψ_k Ṁ_g Δt/M), ψ = 1 inside m/M = 0.5, smooth to 0 at 0.9; surface Lagrangian; no inflow):
#            nothing physical happens, every change of E = U + K + W + P_ext V is spurious.
#   accrete: m_k ∝ M(t) for all k, gas enters at the surface at rest with e_in = e of the top zone; expected ΔE = Σ δm (e_in − G M/R) (ledger as S6a).
# Options: N [80], STEPS [200], FRAC [0.2] (total relative grid drift), DTDYN [3] (dt in t_dyn), CADV [TITAN default], GRAVFIX [1] (TitanResidFast.ALE_GRAVFIX).
using Printf
include(joinpath(@__DIR__, "s4_common.jl"))
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
mode = isempty(ARGS) ? "relabel" : ARGS[1]
const G = Z0.G
n = envi("N", 80); gam = 5 / 3; gmu = 0.6
TF.load_state!(joinpath(PS.INIT, "SOD.tstate.gz")); Rg = TF.getscalar(:ck) / (TF.getscalar(:cm0) * gmu)
# zone masses: geometric outward (ratio q1) in the interior, then geometric decreasing (ratio q2) to resolve the surface
ρc = 1.0; Tc = 1e7; pc = Rg * ρc * Tc; Kp = pc / ρc^gam
nin = round(Int, envf("FIN", 0.5) * n); dμ = [1.0e28 * 1.08^(j - 1) for j in 1:nin]; for j in 1:n-nin; push!(dμ, dμ[nin] * envf("Q2", 0.82)^j); end
# entropy profile: K = Kp in the core, × KR [1] in the envelope, smooth step over KW zones centred on zone fraction KJ: a low-entropy core under a
# high-entropy envelope (convectively stable), the density falls by KR^(3/5) across a few zones as at the edge of the compact core of the S6a star
function Kz(k, nn)
    x = clamp((k - envf("KJ", 0.5) * nn) / envf("KW", 3.0) + 0.5, 0.0, 1.0)
    Kp * envf("KR", 1.0)^(x * x * (3 - 2x))
end
function march(dμ)
    nn = length(dμ); rI = zeros(nn + 1); rho = zeros(nn); p = zeros(nn); μ = [0.0; cumsum(dμ)]
    rho[1] = ρc; p[1] = pc; rI[2] = cbrt(3dμ[1] / ρc)
    for k in 2:nn
        p[k] = p[k-1] - 4π * G * μ[k] / rI[k]^4 * 0.5 * (dμ[k-1] + dμ[k])
        p[k] > 0 || return nothing
        rho[k] = (p[k] / Kz(k, nn))^(1 / gam); rI[k+1] = cbrt(rI[k]^3 + 3dμ[k] / rho[k])
    end
    (; rI, rho, p, μ)
end
# scale the zone masses so that the polytrope ends with p_n ≈ 1e-6 p_c
lo, hi = 1e-4, 1e8
for _ in 1:80
    s = sqrt(lo * hi); m_ = march(dμ .* s)
    (m_ === nothing || m_.p[end] < envf("PEND", 1e-6) * pc) ? (global hi = s) : (global lo = s)
end
dμ .*= lo; st = march(dμ)
T0 = st.p ./ (Rg .* st.rho)
@printf("density ratio across the core edge (zones KJ±KW): %.3e; max dln rho per zone %.3f\n", maximum(st.rho[1:end-1] ./ st.rho[2:end]) , maximum(log.(st.rho[1:end-1] ./ st.rho[2:end])))
if mode == "adaptive"       # TITAN's adaptive grid (lgrid = 1, Dorfi–Drury, resolution in x_me and ln ρ) inside rigid walls: the grid relaxes to its own equilibrium and moves mass between zones
    PS.build!(; rI = st.rI, rho = st.rho, T = T0, u = zeros(n + 1), geom = 2, lgrid = 1, inner = :wall, outer = :wall, selfgrav = true, gam, gmu,
              xmext = envf("XMEXT", 0.3) * st.μ[end], gridtau = envf("GRIDTAU", 1e5), lady = (2, 2, 0), xscale = st.rI[end])
else
    PS.build!(; rI = st.rI, rho = st.rho, T = T0, u = zeros(n + 1), geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = 0.0, selfgrav = true, gam, gmu,
              xmext = envf("XMEXT", 0.3) * st.μ[end])
end
TX.REGULAR_CENTRE[] = true; TX.XMINT[] = Float64[]
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
xmint(k) = TF.getscalar(:xmtot) - TF.var(:xmen)[k]
mode == "adaptive" || (TX.PEXT[] = TF.var(:pgn)[b] - 4π * G * xmint(b + 1) / st.rI[end]^4 * 0.5 * TF.var(:dn)[b] * TF.var(:dvoln)[b])
haskey(ENV, "CADV") && TF.setscalar!(:cadv, envf("CADV", 1.0))
TF.setscalar!(:thet, 1.0); TF.setscalar!(:tmax, 1e60); TF.setscalar!(:conv, 1e-10); TF.setscalar!(:niter, 40)
haskey(ENV, "GRAVWORK") && (TX.GRAVWORK[] = ENV["GRAVWORK"] == "1")
M0 = 4π * xmint(b + 1); R0 = st.rI[end]; tdyn = sqrt(R0^3 / (G * M0))
dt = envf("DTDYN", 3.0) * tdyn; nst = envi("STEPS", 200); frac = envf("FRAC", 0.2)
Mdot = frac * M0 / (nst * dt)
rn, dn, un = TF.var(:rn), TF.var(:dn), TF.var(:un)
"E = U + K + W + P_ext V [erg]; W both as the 8-point Gauss quadrature of −G∫m dm/r (r³ linear in m in a zone; S6a ledger) and as
W_d = −Σ_k G m_k m̄_k/r_k (the potential whose rate is the work of the momentum rows' gravity, m̄ = half-zone masses)."
function energies()
    egn, dvn = TF.var(:egn), TF.var(:dvoln)
    xg = (-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363)
    wg = (0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763)
    dmz = [4π * dn[k] * dvn[k] for k in a:b]
    U = 0.0; K = 0.0; W = 0.0; Wd = 0.0; mc = 0.0
    for k in a:b
        j = k - a + 1; U += dmz[j] * egn[k]; K += dmz[j] * 0.25 * (un[k]^2 + un[k+1]^2)
        r13, r23 = rn[k]^3, rn[k+1]^3
        for (x, w) in zip(xg, wg); f = 0.5 * (1 + x); W -= 0.5 * w * dmz[j] * G * (mc + f * dmz[j]) / cbrt(r13 + f * (r23 - r13)); end
        mc += dmz[j]
        mb = k < b ? 0.5 * (dmz[j] + dmz[j+1]) : 0.5 * dmz[j]
        Wd -= G * 4π * xmint(k + 1) * mb / rn[k+1]
    end
    pV = isnan(TX.PEXT[]) ? 0.0 : TX.PEXT[] * 4π * rn[b+1]^3 / 3
    (; U, K, W, Wd, pV, E = U + K + W + pV, Ed = U + K + Wd + pV)
end
if mode == "relabel"
    TD.ALE_INFLOW[] = false
    psi = zeros(length(TF.var(:xmen)))
    for k in a:b+1; q = xmint(k) / xmint(b + 1); q1 = envf("QR1", 0.5); q2 = envf("QR2", 0.9); x = clamp((q2 - q) / (q2 - q1), 0.0, 1.0); psi[k] = x * x * (3 - 2x); end
    psi[b+1] = 0.0; TD.ALE_PSI[] = psi
elseif mode == "accrete"
    TD.ALE_INFLOW[] = true; TD.ALE_PSI[] = Float64[]
end
TD.ALE_MDOT[] = mode == "adaptive" ? 0.0 : Mdot
mode == "adaptive" && (Mdot = M0 / (nst * dt))       # only the normalisation G M Mdot/R of the printout
TC.DT_RULE[] = dt0 -> dt; TF.setscalar!(:dtime, dt)
E0 = energies(); Eexp = 0.0
ttot = 0.0; Eh = 0.0
@printf("S6a-fix test %s: polytrope gamma 5/3, n = %d zones, M = %.4e g, R = %.4e cm, rho_c/rho_s = %.3e, t_dyn = %.3e s, dt = %.3e s, %d steps, grid drift %.2f M, cadv = %.2f, gravwork %s\n",
        mode, n, M0, R0, dn[a] / dn[b], tdyn, dt, nst, frac, TF.getscalar(:cadv), string(TX.GRAVWORK[]))
@printf("scale of the effect: G M Mdot/R = %.4e erg/s; |W| = %.4e erg\n", G * M0 * Mdot / R0, abs(E0.W))
for k in 1:nst
    Ebef = energies(); Mb = 4π * xmint(b + 1); Rb = rn[b+1]
    if mode == "accrete"
        # inflow: gas at rest with the specific internal energy of the top zone
        e_in = TF.var(:egn)[b]
        TF.setscalar!(:phir0, -Mdot / (4π * Rb^2)); TF.setscalar!(:phir2, -Mdot / (4π * Rb^2) * e_in)
    end
    TD.STEPFN[]()
    if get(ENV, "DEBUGGW", "0") == "1" && k == 1
        c = TX.context(); R_ = zeros(c.neqn, c.nb); ws = TX.getws(Float64, c.mgr); TX.residual!(R_, ws, c, zeros(c.neqn * c.nb)); W_ = ws.a
        Gp = 4π * G; dtt = c.p.dt
        for j in a+1:b
            cj = TX.gravwork(j, c, W_, Float64)
            wk = Gp * W_.xm[j] * 0.5 * (W_.d[j] * W_.dvol[j] + W_.d[j-1] * W_.dvol[j-1]) * W_.urel[j] / W_.rmu[j]
            @printf("  j %3d  C/dt %+.4e  work_urel %+.4e  dm/dt %+.4e  urel %+.3e  rho_b/rho_bar %.4f  rho %.3e\n", j - a, cj, wk, (W_.xmen[j] - c.xmeo[j]) / -dtt, W_.urel[j], W_.db[j] / (0.5 * (W_.d[j] + W_.d[j-1])), W_.d[j])
        end
    end
    Eaft = energies(); Ma = 4π * xmint(b + 1); dM = Ma - Mb
    dtk = TF.getscalar(:timen) - TF.getscalar(:timeo); global ttot += dtk
    if mode == "accrete"
        global Eexp += dM * (TF.getscalar(:phir2) / TF.getscalar(:phir0) - G * 0.5 * (Mb + Ma) / (0.5 * (Rb + rn[b+1])))
        global Eh += dM * (gam - 1) * TF.getscalar(:phir2) / TF.getscalar(:phir0)      # flow work P/ρ of the entering gas (not supplied by the e_in inflow)
    end
    if k % (nst ÷ 10) == 0 || k == 1
        @printf("step %4d M/M0 = %.5f (dM/Mdot t = %.4f) R/R0 = %.6f max|u|/c_s = %.2e  dE = %+.4e  dE_d = %+.4e  (expected %+.4e)  [dE - exp]/(G M Mdot/R t) = %+.4e (Gauss) %+.4e (W_d)\n",
                k, Ma / M0, (Ma - M0) / (Mdot * ttot), rn[b+1] / R0, maximum(abs(un[kk]) / TF.var(:asn)[min(kk, b)] for kk in a:b+1), Eaft.E - E0.E, Eaft.Ed - E0.Ed, Eexp,
                (Eaft.E - E0.E - Eexp) / (G * M0 * Mdot / R0 * ttot), (Eaft.Ed - E0.Ed - Eexp) / (G * M0 * Mdot / R0 * ttot))
    end
end
Ef = energies()
@printf("RESULT %s: (dE - expected)/(G M Mdot/R t) = %+.4e (Gauss W), %+.4e (W_d); dU/(GMMdot/R t) = %+.4e dW = %+.4e dpV = %+.4e; rho_c/rho_c0 - 1 = %+.3e; with the flow work P/rho of the inflow in the expectation: %+.4e; (dE - expected)/|W| = %+.3e (Gauss) %+.3e (W_d)\n", mode,
        (Ef.E - E0.E - Eexp) / (G * M0 * Mdot / R0 * ttot), (Ef.Ed - E0.Ed - Eexp) / (G * M0 * Mdot / R0 * ttot),
        (Ef.U - E0.U) / (G * M0 * Mdot / R0 * ttot), (Ef.W - E0.W) / (G * M0 * Mdot / R0 * ttot), (Ef.pV - E0.pV) / (G * M0 * Mdot / R0 * ttot), dn[a] / st.rho[1] - 1, (Ef.E - E0.E - Eexp + Eh) / (G * M0 * Mdot / R0 * ttot), (Ef.E - E0.E - Eexp) / abs(E0.W), (Ef.Ed - E0.Ed - Eexp) / abs(E0.Wd))
