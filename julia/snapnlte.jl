# Non-LTE H I, He I, He II on a recorded snapshot (R2 of docs/POPIII_RADIATION_PLAN.md; post-processing, feeds nothing back).
#
#   julia --project=<env> -e 'include(ARGS[1])' snapnlte.jl <run dir> [step | last] [nmaxH=20] [nmaxHe=20] [out.csv]
#
# Per zone outside the star (r > R★ from perf.csv, T < 3e5 K): statistical equilibrium of the three atoms of ModelAtom at the
# snapshot's T and ρ, with
#   - the radiation field the run itself uses for its photo-rates (PopIIIMaterial.starlight_setup): J_ν = B_ν(T_r),
#     T_r = (E_r/a)^¼, inside the photosphere; J_ν = W B_ν(T_eff), W the geometric dilution, outside;
#   - Sobolev escape probabilities in the snapshot's velocity field (zone-centre v, dv/dr from the interfaces), with a thermal-
#     width floor v_th/Δr on the velocity gradient where the flow is nearly static; Λ-iteration on β;
#   - n_e from H⁺, He⁺, 2He²⁺ of the atoms themselves (iterated).
# Reports per zone the ionisation state and departure coefficients, and line luminosities: escaping from each zone (Sobolev)
# and attenuated by the grey continuum optical depth to the surface, exp(−τ_χ) with χ the run's flux-mean opacity (APPROXIMATE:
# grey; to be replaced by the frequency-dependent continuum transfer of Spectrum.jl).
using Printf
const ROOT = @__DIR__
include(joinpath(ROOT, "runrec.jl"))
isdefined(Main, :Z0Physics) || include(joinpath(ROOT, "Z0Physics.jl"))
include(joinpath(ROOT, "ModelAtom.jl"))
const MA = ModelAtom

function lastperf(dir, step)
    best = nothing
    for ln in eachline(joinpath(dir, "perf.csv"))
        startswith(ln, "step") && continue
        p = split(ln, ","); k = parse(Int, p[1]); M = parse(Float64, p[25]); M > 1000 && continue
        k <= step && (best = (M, parse(Float64, p[26])))
    end
    best
end

const aRAD = 7.5657e-15; const RSUN = 6.957e10; const LSUN = 3.828e33; const mH_ = 1.6735575e-24

"Solve one atom in its rate matrix with Λ-iteration on the Sobolev β; returns (n, nplus, β)."
function solve_atom(a, L, T, ne, ntot, J, v, r, dvdr, Qmin; nit = 20)
    P0 = MA.static_rates(a, T, ne; J); β = ones(length(L)); sol = nothing
    for it in 1:nit
        P = MA.add_lines!(copy(P0), L, β; J); sol = MA.solve_P(P, ntot)
        βn = similar(β)
        for (k, (l, u, ν, A, Bul, Blu, f)) in enumerate(L)
            K = MA.BCOEF * f * (MA.c / ν) * max(sol.n[l] - a.levels[l].g / a.levels[u].g * sol.n[u], 0.0)
            βn[k] = MA.βsob(K, v, r, dvdr; Qmin)
        end
        δ = maximum(abs.(βn .- β)); β = βn
        δ < 1e-4 && break
    end
    (sol.n, sol.nplus, β)
end
emis(a, L, n, β, l, u) = (k = findfirst(x -> x[1] == l && x[2] == u, L); k === nothing ? 0.0 : n[u] * L[k][4] * β[k] * MA.h * L[k][3])

function main(args)
    dir = args[1]; which = length(args) >= 2 ? args[2] : "last"
    nH_ = length(args) >= 3 ? parse(Int, args[3]) : 20; nHe = length(args) >= 4 ? parse(Int, args[4]) : 20
    out = length(args) >= 5 ? args[5] : ""
    V = lineage(dir); s = which == "last" ? V[end] : V[findlast(x -> x.step <= parse(Int, which), V)]
    Mst, Rst = lastperf(dir, s.step)
    H = MA.hydrogenic(1; nmax = nH_, coll = :chianti, chianti = joinpath(ROOT, "atomdata", "chianti", "h_1"))
    He2 = MA.hydrogenic(2; nmax = nH_, coll = :chianti, chianti = joinpath(ROOT, "atomdata", "chianti", "he_2"))
    He1 = MA.heI(; nmax = nHe)
    LH, L2, L1 = MA.lines(H), MA.lines(He2), MA.lines(He1)
    lab1 = [l.label for l in He1.levels]; t1(x) = findfirst(==(x), lab1)
    yHe = 1 / 12; μH = mH_ * (1 + 4yHe)
    nz = length(s.T); kph = s.kph - 4; rph = s.r[kph+1]
    # grey optical depth from each zone's centre to the surface
    τ = zeros(nz); acc = 0.0
    for k in nz:-1:1; dτ = s.chi[k] * s.ρ[k] * (s.r[k+1] - s.r[k]); τ[k] = acc + dτ / 2; acc += dτ; end
    linedefs = (("H Lya 1216", :H, 1, 2), ("H Halpha 6563", :H, 2, 3), ("H Hbeta 4861", :H, 2, 4),
                ("HeII 304", :He2, 1, 2), ("HeII 1640", :He2, 2, 3), ("HeII 4686", :He2, 3, 4),
                ("HeI 584", :He1, t1("1s2 1S"), t1("1s.2p 1P")), ("HeI 10830", :He1, t1("1s.2s 3S"), t1("1s.2p 3P")),
                ("HeI 5876", :He1, t1("1s.2p 3P"), t1("1s.3d 3D")), ("HeI 6678", :He1, t1("1s.2p 1P"), t1("1s.3d 1D")))
    Lloc = zeros(length(linedefs)); Lobs = zeros(length(linedefs))
    io = out == "" ? nothing : open(out, "w")
    io === nothing || println(io, "k,r_Rsun,T,rho,ne,xHp,xHe0,xHe1,xHe2,b2_H,b1_HeI_2t3S,tau_grey")
    @printf("snapshot step %d, t %.6e yr, M* %.3f Msun, R* %.2f Rsun, photosphere zone %d at %.1f Rsun, T_eff %.0f K\n",
            s.step, s.t / 3.156e7, Mst, Rst, kph, rph / RSUN, s.Teff)
    @printf("%4s %9s %9s %9s %9s %8s %8s %8s %8s %9s %8s\n", "zone", "r/Rsun", "T", "n_H", "n_e", "x(H+)", "He0", "He+", "He++", "b2(H)", "tau")
    tw = time()
    for k in 1:nz
        rc = 0.5 * (s.r[k] + s.r[k+1]); T = s.T[k]
        (rc > Rst * RSUN && T < 3e5 && T > 2000) || continue
        nHk = s.ρ[k] / μH; nHek = yHe * nHk
        v = 0.5 * (s.u[k] + s.u[k+1]); dvdr = (s.u[k+1] - s.u[k]) / (s.r[k+1] - s.r[k]); Δr = s.r[k+1] - s.r[k]
        Tr = (max(s.er[k], 1e-30) / aRAD)^0.25
        J = k <= kph ? (ν -> MA.Bν(ν, Tr)) : (W = 0.5 * (1 - sqrt(max(1 - (rph / rc)^2, 0.0))); ν -> W * MA.Bν(ν, s.Teff))
        ne = nHk; nH0 = (0.0, 0.0); res = nothing
        ntHe1 = nHek; ntHe2 = nHek
        for it in 1:6
            vth(m) = sqrt(2 * MA.kB * T / m)
            nh, nhp, βh = solve_atom(H, LH, T, ne, nHk, J, v, rc, dvdr, vth(mH_) / Δr)
            n1, n1p, β1 = solve_atom(He1, L1, T, ne, ntHe1, J, v, rc, dvdr, vth(4mH_) / Δr)
            n2, n2p, β2 = solve_atom(He2, L2, T, ne, ntHe2, J, v, rc, dvdr, vth(4mH_) / Δr)
            f1 = n1p / max(sum(n1), 1e-300); f2 = n2p / max(sum(n2), 1e-300)
            He0 = nHek / (1 + f1 + f1 * f2); He1_ = He0 * f1; He2_ = He1_ * f2
            ntHe1 = He0 + He1_; ntHe2 = He1_ + He2_
            ne_new = nhp + He1_ + 2He2_
            res = (; nh, nhp, βh, n1, β1, n2, β2, He0, He1_, He2_)
            abs(ne_new / ne - 1) < 1e-3 && (ne = ne_new; break)
            ne = sqrt(ne * ne_new)
        end
        dV = 4π / 3 * (s.r[k+1]^3 - s.r[k]^3)
        for (i, (nm, at, l, u)) in enumerate(linedefs)
            e = at === :H ? emis(H, LH, res.nh, res.βh, l, u) : at === :He2 ? emis(He2, L2, res.n2, res.β2, l, u) :
                emis(He1, L1, res.n1, res.β1, l, u)
            Lloc[i] += e * dV; Lobs[i] += e * dV * exp(-τ[k])
        end
        b2 = res.nh[2] / (res.nhp * ne * MA.Φ(H, 2, T))
        @printf("%4d %9.4g %9.4g %9.3e %9.3e %8.2e %8.2e %8.2e %8.2e %9.3e %8.2e\n", k, rc / RSUN, T, nHk, ne, res.nhp / nHk,
                res.He0 / nHek, res.He1_ / nHek, res.He2_ / nHek, b2, τ[k])
        io === nothing || @printf(io, "%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", k, rc / RSUN, T, s.ρ[k], ne, res.nhp / nHk,
                                  res.He0 / nHek, res.He1_ / nHek, res.He2_ / nHek, b2, res.n1[2] / max(sum(res.n1), 1e-300), τ[k])
    end
    io === nothing || close(io)
    Lph = 4π * rph^2 * 5.670374e-5 * s.Teff^4
    @printf("\nline luminosities [L_sun] (zones solved in %.0f s); continuum L(photosphere) = %.3e L_sun\n", time() - tw, Lph / LSUN)
    @printf("%-16s %12s %12s %12s\n", "line", "escaping", "after grey τ", "/ L_ph")
    for (i, (nm, _, _, _)) in enumerate(linedefs)
        @printf("%-16s %12.4e %12.4e %12.3e\n", nm, Lloc[i] / LSUN, Lobs[i] / LSUN, Lobs[i] / Lph)
    end
end
main(ARGS[2:end])
