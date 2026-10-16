# Export the opacity table and reference values for the report's Rust/WASM core (report/core).
#   julia --project=. scripts/export_reference.jl ../report/core/data
# Writes
#   kaptable.bin    "QSKAP1\0\0", u32 nT, u32 nR, f64 lT0, dlT, lR0, dlR, f32 log10 κR [nT][nR] (R fastest), f32 log10 κP [nT][nR]
#   reference.json  values of this (Julia) implementation that the Rust tests must reproduce
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Printf
out = ARGS[1]; mkpath(out)
tab = load_table()

# ---------------------------------------------------------------- table
let nT = length(tab.lT), nR = length(tab.lR)
    open(joinpath(out, "kaptable.bin"), "w") do io
        write(io, b"QSKAP1\0\0", UInt32(nT), UInt32(nR),
              tab.lT[1], tab.lT[2] - tab.lT[1], tab.lR[1], tab.lR[2] - tab.lR[1])
        for M in (tab.lk, tab.lkP), i in 1:nT, j in 1:nR
            write(io, Float32(M[i, j]))
        end
    end
    @printf("table %d × %d  lT %.3f:%.4f:%.3f  lR %.3f:%.4f:%.3f\n", nT, nR, tab.lT[1], tab.lT[2] - tab.lT[1], tab.lT[end],
            tab.lR[1], tab.lR[2] - tab.lR[1], tab.lR[end])
end

# ---------------------------------------------------------------- reference values (hand-written JSON, no package)
js(x::Real) = isfinite(x) ? @sprintf("%.17g", x) : "null"
js(x::AbstractString) = "\"" * x * "\""
js(x::Symbol) = js(String(x))
js(v::AbstractVector) = "[" * join(js.(v), ",") * "]"
js(d::NamedTuple) = "{" * join(["\"$k\":" * js(v) for (k, v) in pairs(d)], ",") * "}"
js(d::Dict) = "{" * join(["\"$k\":" * js(v) for (k, v) in d], ",") * "}"

ρT = [(1e-13, 6000.0), (3.7e-13, 7300.0), (1e-12, 1.0e4), (3e-12, 2.5e4), (1e-10, 1e5), (2e-9, 2.6e5),
      (1e-6, 3e6), (1e-16, 6200.0), (1e-9, 9000.0)]
saha_ref = [(s = Q.saha(ρ, T); (rho = ρ, T = T, ne = s.ne, xH = s.xH, nHe1 = s.nHe1, nHe2 = s.nHe2)) for (ρ, T) in ρT]
eos_ref = [(st = Q.eos_rhoT(ρ, T); (rho = ρ, T = T, P = st.P, e = st.e, s = st.s, nabla_ad = st.∇ad, Gamma1 = st.Γ1,
            delta = st.δ, cP = st.cP, cs = st.cs, beta = st.β)) for (ρ, T) in ρT]
PT = [(10.0, 7000.0), (0.5 + Q.Prad(8000.0), 8000.0), (1e4, 3e4), (1e7, 2.6e5), (100.0, 1.2e4)]
rhoPT_ref = [(P = P, T = T, rho = Q.rho_PT(P, T)) for (P, T) in PT]
kap_ref = [(rho = ρ, T = T, kR = Q.kap(tab, ρ, T), kP = Q.kapP(tab, ρ, T)) for (ρ, T) in vcat(ρT, [(1e-14, 3000.0), (1e-8, 1e7)])]

M = 1e6Msun; L = LEdd_es(M)
p01 = QSParams(M = M, MBH = 0.1M, L = L, tab = tab)
mlt_ref = []
pts = let C = [parse.(Float64, split(l, ",")) for l in readlines(joinpath(@__DIR__, "..", "output", "profiles", "qs_f0.1_md0.csv"))[3:end]]
    [(C[i][3], C[i][4], C[i][1], C[i][2]) for i in (1, 5, 20, 40, 50, 60, 80, 100)]          # (P, T, r, m) on the steady profile
end
for (P, T, r, m) in pts
    q = Q.mlt(P, T, r, m, L, p01)
    length(q) == 1 && (println("mlt: no state at P=$P T=$T"); continue)
    push!(mlt_ref, (P = P, T = T, r = r, m = m, nabla = q.∇, nabla_rad = q.∇rad, nabla_ad = q.∇ad, v = q.v, rho = q.ρ, kappa = q.κ, Gamma = q.Γ))
end
bc_ref = []
for (fBH, md, R) in [(0.1, 0.0, 8.47e15), (0.1, 0.3, 8.44e15), (0.1, 1.0, 8.27e14), (0.01, 3.0, 1.86e15)]
    p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
    bc = outer_bc(R, p; Mdot = md * Msun / yr)
    push!(bc_ref, (fBH = fBH, Mdot_Msun_yr = md, R = R, Pb = bc.Pb, Tb = bc.Tb, tau_in = bc.τ_in, r_ph = bc.r_ph_infall, Gamma_kappa = bc.Γκ))
end
shoot_ref = []
for (fBH, md, R) in [(0.1, 0.0, 8.0e15), (0.1, 0.0, 9.0e15), (0.1, 1.0, 5e15)]
    p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
    bc = outer_bc(R, p; Mdot = md * Msun / yr)
    st, mrin, prof = integrate_inward(R, bc.Pb, bc.Tb, p)
    push!(shoot_ref, (fBH = fBH, Mdot_Msun_yr = md, R = R, status = st, m_rin = mrin))
end
# all steady solutions (the published table, output/branches.csv, ε_in = 1e-3, scan 14:0.025:16.6 then bisection)
branches = [split(l, ",") for l in readlines(joinpath(@__DIR__, "..", "output", "branches.csv"))[2:end]]
branch_ref = [(fBH = parse(Float64, b[1]), Mdot_Msun_yr = parse(Float64, b[2]), R = parse(Float64, b[3]), Tb = parse(Float64, b[5]),
               T_obs = parse(Float64, b[10])) for b in branches]
# the same scan + bisection as scripts/branches.jl, recomputed with the CURRENT code (output/branches.csv may predate
# code changes): roots, number of roots and the scan signs for every (fBH, Ṁ) of branches.jl.
# `tab32` is the table rounded to Float32 exactly as stored in kaptable.bin: results with it isolate the port from
# the table rounding (the Rust core reads the Float32 table).
tab32 = Q.KapTable(tab.lT, tab.lR, Float64.(Float32.(tab.lk)), Float64.(Float32.(tab.lkP)))
function all_branches(tb)
    res = []
    lRs = collect(14.0:0.025:16.6)
    for fBH in (0.1, 0.01), md in [0.0, 0.3, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 30.0, 100.0]
        p = QSParams(M = M, MBH = fBH * M, L = L, tab = tb)
        Mdot = md * Msun / yr
        sg = [shoot_sign(10^lR, p; Mdot = Mdot)[1] for lR in lRs]
        Rs = Float64[]; Tbs = Float64[]; Tobs = Float64[]
        for i in 1:length(lRs)-1
            (sg[i] == 0 || sg[i+1] == 0 || sg[i] == sg[i+1]) && continue
            a, b, sa = lRs[i], lRs[i+1], sg[i]
            for _ in 1:40
                mid = 0.5(a + b); s = shoot_sign(10^mid, p; Mdot = Mdot)[1]
                s == 0 && (b = mid; continue)
                s == sa ? (a = mid) : (b = mid)
            end
            R = 10^(0.5(a + b)); bc = outer_bc(R, p; Mdot = Mdot)
            push!(Rs, R); push!(Tbs, bc.Tb); push!(Tobs, (bc.Ltot / (4π * σSB * bc.r_ph_infall^2))^0.25)
        end
        push!(res, (fBH = fBH, Mdot_Msun_yr = md, nroots = length(Rs), R = Rs, Tb = Tbs, T_obs = Tobs, signs = sg))
        @printf("branches (current code%s) fBH=%.2f Ṁ=%g: %s\n", tb === tab ? "" : ", Float32 table", fBH, md,
                join([@sprintf("%.9e", r) for r in Rs], ", "))
    end
    return res
end
branch_now = all_branches(tab)
branch_now32 = all_branches(tab32)
# full every-step profile with the Float32 table, integrated from the Float64-table root of solve_R below
prof32 = let p = QSParams(M = M, MBH = 0.1M, L = L, tab = tab32), R = 8.435236237964494e15
    bc = outer_bc(R, p); st, mrin, pr = integrate_inward(R, bc.Pb, bc.Tb, p; save = true)
    (R = R, status = st, m_rin = mrin, r = pr[:r], m = pr[:m], P = pr[:P], T = pr[:T], rho = pr[:ρ], nabla = pr[:∇],
     v = pr[:v], Gamma = pr[:Γ], beta = pr[:β])
end
# one full profile for plotting comparisons
s = solve_R(p01; lRmin = 15.0, lRmax = 16.8, nscan = 37)
prof = (r = s.prof[:r][1:4:end], T = s.prof[:T][1:4:end], rho = s.prof[:ρ][1:4:end], Gamma = s.prof[:Γ][1:4:end])

consts = (G = Q.G, c = Q.c, kB = Q.kB, h = Q.h, me = Q.me, mH = Q.mH, mHe = Q.mHe, sigmaSB = Q.σSB, arad = Q.arad,
          Msun = Q.Msun, yr = Q.yr, chiH = Q.χH, chiHe1 = Q.χHe1, chiHe2 = Q.χHe2, X = Q.COMP.X, Y = Q.COMP.Y,
          M = M, L = L, alpha = p01.αMLT, eps_in = p01.ε_in)
open(joinpath(out, "reference.json"), "w") do io
    println(io, "{\"constants\":", js(consts), ",")
    println(io, "\"saha\":", js(saha_ref), ",")
    println(io, "\"eos\":", js(eos_ref), ",")
    println(io, "\"rho_PT\":", js(rhoPT_ref), ",")
    println(io, "\"kap\":", js(kap_ref), ",")
    println(io, "\"mlt\":", js(mlt_ref), ",")
    println(io, "\"outer_bc\":", js(bc_ref), ",")
    println(io, "\"shoot\":", js(shoot_ref), ",")
    println(io, "\"branches\":", js(branch_ref), ",")
    println(io, "\"branches_current\":", js(branch_now), ",")
    println(io, "\"branches_current_f32table\":", js(branch_now32), ",")
    println(io, "\"profile_f01_md0_f32table\":", js(prof32), ",")
    println(io, "\"profile_f01_md0\":", js((R = s.R, r = prof.r, T = prof.T, rho = prof.rho, Gamma = prof.Gamma)), "}")
end
println("reference.json written; R(f=0.1, Ṁ=0) = ", s.R)
