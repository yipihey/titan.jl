# PopIIIRun: the Pop III driver (the former body of popiii_m2c.jl) as a module, so that its function definitions are part of the
# precompiled TitanStack image (docs/POPIII_PRECOMPILE.md). popiii_m2c.jl is the thin script that loads the module files (or the
# image) and calls `PopIIIRun.main(ARGS)`.
#
# What is here: every function / struct definition of the driver, its constant bookkeeping state (Refs, counters), and the
# declarations of the driver's global variables. What is not: the run itself. The statements that read ENV and ARGS, build the
# initial condition and set the solver switches are in popiii_run_body.jl, which `main` evaluates in this module at run time, in
# the original order, so the settings of the run are read when it starts and are never part of an image. The only load-time
# environment read here is HOT (HOTON), a key the TitanStack guard records anyway (the HOT includes are conditional on it).
# Settings that were `const X = get(ENV, ...)` in the script are now typed globals assigned by the body (same values; the
# declared type keeps the reads type-stable).
module PopIIIRun
using Printf
import ..TitanDriver, ..JuliaPhysics, ..TitanResidFast, ..BlockSolver, ..TitanCore, ..Z0Physics, ..TitanMaterials, ..ProblemSetup, ..PrimordialChem
import ..SCVH, ..PopIIIEqEOS, ..ContinuumOpacity, ..G2Opacity, ..ChemNet, ..ChemRates, ..PopIIIMaterial, ..H2Lines, ..StellarOpacity, ..IonRay
import ..Trbdf2, ..NewtonKrylov, ..StageCons, ..H2Levels, ..GridEq
# modules that exist only in some configurations (MGROUPS: JuliaSolver; HOT=1: the stellar modules)
for nm in (:JuliaSolver, :StellarEOS, :NucNet, :StellarMaterial, :StarMaterial)
    isdefined(parentmodule(@__MODULE__), nm) && Core.eval(@__MODULE__, :(const $nm = $(getfield(parentmodule(@__MODULE__), nm))))
end
const HOTON = get(ENV, "HOT", "0") == "1"      # S6-C1: HOT=1: stellar EOS, nuclear network and convective mixing in the hot core (T > HOTT1), see PopIIIMaterial.HOT

# --- the driver's global variables that the functions below read (assigned by popiii_run_body.jl at run time). They are declared with a value:
# a bare `global x` would make the first assignment a binding change that invalidates every cached method reading x. Other names the script
# assigns are created by the body (the functions do not read them).
global tcmb = nothing
global cf = nothing
global DH = nothing
global SPLIT = nothing
global STEP0 = nothing
global CHEMDIAG::Bool = false
global CKDIR::String = ""
global CKEVERY::Int = 0
global CONVDIAG::Bool = false     # CONVDIAG=1 (set by the body): snapshots carry `cd` (convdiag_snap)
global CPLPROF::Bool = false
global CPLPROFEVERY::Int = 0
global DTSTOP::Float64 = 0.0
global TDC::Bool = false          # TDC=1 (set by the body): turbulent velocity w as the last Newton unknown
global ENVCENTRE::Bool = false
global ENVCORE::Vector{Float64} = Float64[]
global ENVEUL = nothing
global ENVLIN::Int = 0
global ENVSINKTAU::Float64 = 0.0
global FLUSHOUT::Bool = false
global HLL = nothing
global HLV = nothing
global LINECORE::Int = 0
global LINEINC::String = ""
global LTAUCUT::Float64 = 0.0
global M = nothing
global Mg = nothing
global Mjoin = nothing
global NANSCAN::Bool = false
global NBIN::Int = 0
global PRINTEVERY::Int = 0
global PROFALLOC::Float64 = 0.0
global PROFSAMPLE::Int = 0
global PROFSAMPLEN::Int = 0
global RAD::Bool = false
global REC = nothing
global RUNREC::String = ""
global SNAPEVERY::Int = 0
global SNAPFILE::String = ""
global SNAPNEXT::Base.RefValue{Float64} = Ref(0.0)
global SNAPS = nothing
global TIMING::Bool = false
global TMARK::Base.RefValue{Float64} = Ref(0.0)
global T_tab = nothing
global Tfin_ = nothing
global Tstop = nothing
global XADV2::Bool = false
global a = nothing
global av = nothing
global b = nothing
global cool = nothing
global dm1 = nothing
global fH2fin_ = nothing
global fHIfin_ = nothing
global ic = nothing
global io = nothing
global lM_tab = nothing
global lf_tab = nothing
global lmxe_fin = nothing
global lxe_fin = nothing
global mat = nothing
global mcore_fin = nothing
global mv1_fin = nothing
global netw = nothing
global nextln = nothing
global nsh = nothing
global nstop = nothing
global on98 = nothing
global rM1fin = nothing
global rg = nothing
global tff = nothing
global v_tab = nothing
global vfin_ = nothing
global vv1_fin = nothing
global xe0 = nothing
global xefin_0 = nothing
global xefin_flat = nothing


# GRIDEQ settings (assigned by grideq_driver.jl, evaluated by the body)
global GRIDEQP::Union{Nothing,Int} = nothing
global GEQW::Vector{Float64} = Float64[]
global GEQC::Vector{Float64} = Float64[]
global GEQX::Vector{Float64} = Float64[]
global GEQGP::Vector{Float64} = Float64[]
global GEQDEC::Vector{Float64} = Float64[]
global GEQFL::Tuple{Float64,Float64} = (0.0, 0.0)
global GEQRAMP::Union{Nothing,Vector{Float64}} = nothing

# --- definitions (in the order of the former script) ---
const TD = TitanDriver
const TF = TD.TF
const TX = TitanResidFast
const PS = ProblemSetup
const TC = TitanCore
const PC = PrimordialChem
const PM = PopIIIMaterial
const Z0 = Z0Physics
using ForwardDiff
const Msun = 1.98847e33
pr(f, a...) = (s = Printf.format(Printf.Format(f), a...); print(s); print(io, s); flush(io))
q_of(M, n) = let lo = 1.0001, hi = 2.0
    f(q) = dm1 * (q^n - 1) / (q - 1) - M
    for _ in 1:200; m_ = (lo + hi) / 2; f(m_) > 0 ? (hi = m_) : (lo = m_); end
    (lo + hi) / 2
end
const NIMP = Ref(0)
# (kept in the checkpoint; the IMPCOOL slope it counted is gone)
"Composition of zone k as the PrimordialChem tuple (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), for either material."
ycomp(k) = mat isa PM.NetGas ? (d = ChemNet.full(mat.x[k], mat.net); (d[:Hp], d[:Hm], d[:H2], d[:H2p])) : mat.y[k]
"Equilibrium-regime flag of zone k."
eqflag(k) = mat isa PM.NetGas ? mat.eqzone[k] : TF.var(:dn)[k] >= mat.ρeq
on98T(ln) = (j = findfirst(p -> p[1] >= ln, on98); j === nothing ? NaN : j == 1 ? 10^on98[1][2] :
             10^(on98[j-1][2] + (ln - on98[j-1][1]) / (on98[j][1] - on98[j-1][1]) * (on98[j][2] - on98[j-1][2])))
using Serialization
"ENV keys any of the julia/ sources reads (only these go into the meta file)."
function recenvkeys()
    ks = Set{String}()
    for f in readdir(@__DIR__; join = true)
        endswith(f, ".jl") || continue
        for m in eachmatch(r"ENV(?:, |\[)\"([A-Z][A-Z0-9_]*)\"", read(f, String)); push!(ks, m[1]); end
    end
    ks
end
function recopen!(k)
    tag = @sprintf("%09d_%d", k, getpid())
    REC.io[] = open(joinpath(RUNREC, "snaps_$(tag).jls"), "w")
    open(joinpath(RUNREC, "meta_$(tag).txt"), "w") do io
        gitc(c) = try strip(read(Cmd(vcat(["git", "-C", TF.ROOT], c)), String)) catch; "?" end
        println(io, "start_step ", k, "\nstart_t_s ", TF.getscalar(:timen), "\nhost ", gethostname(), "\nthreads ", Threads.nthreads(), "\njulia ", VERSION)
        println(io, "git ", gitc(["rev-parse", "HEAD"]), " ", gitc(["rev-parse", "--abbrev-ref", "HEAD"]))
        println(io, "git_dirty_julia ", isempty(gitc(["status", "--short", "julia"])) ? "no" : "yes")
        ks = recenvkeys()
        for key in sort(collect(keys(ENV))); key in ks && println(io, "env ", key, "=", ENV[key]); end
    end
    pf = joinpath(RUNREC, "perf.csv")
    isfile(pf) || open(pf, "w") do io
        println(io, "step,t_yr,dt_s,wall_s,pid,nsteps,t_step,t_chem,t_lines,t_other,gc_s,alloc_GB,ndiv,nrepeat,tr_steps,fail_tr,fail_bdf,be_fallback,nali,nchemfb,nthreads,nz,log_nc,Tc,mstar,rstar,umin_kms,dtlim,etot_Gmax,etot_vmax,etot_zone")
    end
end
"""
CONVDIAG: the diagnostics `cd` of one snapshot, at the accepted state (a NamedTuple of Float64 vectors; zone vectors have length n = b−a+1, interface vectors n+1 with
index j = k−a+1 for the interface k between zones k−1 and k, NaN at j = 1, n+1 and where MLT is not evaluated). The solver is not touched: `TX.conv_diag_pass` is one Float64
residual evaluation (as `TX.step_diag`) with the material caches saved and restored around it. Definitions (docs/POPIII_CONVECTION.md M1):
- nab, nab_raw, hs_w, nad_used, kap, fpot, vpot, fconv, vconv, lmix, HP, w_tau, w_inf: what the MLT block used (nab_raw before the MLT_HS blend, hs_w = its weight sw, 1 = untouched;
  fpot/vpot the flux and velocity with w_tau = w_inf = 1 (0 where ∇ ≤ ∇_ad,used), fconv/vconv with the gate weights; w_tau, w_inf the two factors of cvw;
  fpot_th (-> Fpot_over_Ftot): all gates open, i.e. also the thermodynamic gate: where ∇_ad,used = 1e30 the frozen ∇_ad, c_P, δ of conv_thermo in its place, with the same ∇, ℓ, H_P, κ).
- nad_eq, nad_fr, del_*: `PM.conv_thermo` at the interface mean (ρ, T) with the composition of zone k (the one MLT uses); w_th = 1 if nad_used < 1e29 (network zones: 0).
- nab_rad = 3 κ ρ F_rad H_P/(4 a c T⁴) (≡ 3κLP/(16π a c G m T⁴) with L = 4π r² F_rad, H_P = P/(ρ g)), F_rad the radiative flux of the interface (frn) actually present, κ the Rosseland κ MLT uses.
- N² of a zone-centre difference (no division by Δln P): N²_S = g_eff δ_fr (dlnT/dr − ∇_ad,fr dlnP/dr), N2_eq the same with the equilibrium ∇_ad, δ,
  N²_L = N²_S − g_eff φ dlnμ/dr; g_eff = −(P_k − P_k−1)/(ρ_i Δr), P = p_gas + f E_r; derivatives between the zone centres k−1, k (centre = mean of the faces). nab_mu = Δlnμ/ΔlnP, nab_L = nad_fr + φ/δ_fr nab_mu.
- mach (zone) = |½(u_k+u_k+1)|/c_s, c_s = √(5/3 p_gas/ρ); mach_i at the interface (c_s of the mean of the two zones). t_settle = H_P/|u_i|, t_N = 1/√|N²_L|, t_turn = ℓ/v (potential v where fconv = 0), t_chem (coupled zones) = 1/max|∂S_i/∂x_i|.
- e_*: the terms of the gas-energy row of the zone [erg s⁻¹ sr⁻¹], e_res their sum (+ the pieces not listed), at the accepted state.
- layer: [r_in/R*, r_out/R*, m_layer [Msun], max(−min(N²_L,N²_eq)) t_dyn², max Fconv/Ftot, max Fpot/Ftot, n_interfaces, R* [Rsun]]; the "sub-shock unstable layer" = the run of interfaces just inside the shock
  (R*: the first interface from the centre with u < −2 km/s, as runrec.jl), contiguous inwards, with mach_i < 0.5 and N²_L < 0 or N²_eq < 0; its bounds are the outermost and innermost of those interfaces, the mass is between them.
"""
function convdiag_snap(k, snap)
    mat isa PM.NetGas && TX.CONVECTION[] || return nothing
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    sv = PM.conv_cache_save(); local ps
    try ps = TX.conv_diag_pass(dt) finally PM.conv_cache_restore!(sv) end
    n = b - a + 1; NaNv(m) = fill(NaN, m)
    ρ, T, pg_, er_, fe_, u, r, m_ = snap.ρ, snap.T, snap.pg, snap.er, snap.fedd, snap.u, snap.r, snap.m
    P = pg_ .+ fe_ .* er_; rc = 0.5 .* (r[1:n] .+ r[2:n+1]); cs = sqrt.(5 / 3 .* pg_ ./ ρ)
    mach = [abs(0.5 * (u[j] + u[j+1])) / cs[j] for j in 1:n]
    zth = Vector{Any}(nothing, n); mu = NaNv(n); br = NaNv(n); tch = NaNv(n)
    for j in 1:n
        kk = a + j - 1
        try th = PM.conv_thermo(mat, kk, ρ[j], T[j]); mu[j] = th.mu; br[j] = th.branch catch; end
        if TX.cplz(kk) && !mat.eqzone[kk]
            try
                J = ForwardDiff.jacobian(z -> PM.cpl_g(mat, kk, z), [ρ[j], T[j], TX.XSN[][1][kk], TX.XSN[][2][kk]])
                tch[j] = 1 / max(abs(J[3, 3]), abs(J[4, 4]))
            catch; end
        end
    end
    I(nm) = (v = NaNv(n + 1); for j in 2:n; v[j] = ps[nm][a+j-1]; end; v)       # pass quantity of interface k = a+j−1
    Z(nm) = (v = NaNv(n); for j in 1:n; v[j] = ps[nm][a+j-1]; end; v)
    names = (:nab, :nab_raw, :hs_w, :nad_used, :cp_used, :del_used, :fpot, :vpot, :fconv, :vconv, :lmix, :HP, :kap, :frad, :g, :sup, :gam1b, :gam1_0, :gam1_inf, :wtau, :tchem, :delb, :cpb, :kap_grey)
    d = Dict{Symbol,Vector{Float64}}(nm => I(nm) for nm in names)
    for nm in (:w_tau, :w_inf, :tau); d[nm] = I(nm); end
    for nm in (:e_stor, :e_pdv, :e_av, :e_kin, :e_adv, :e_dfe, :e_rad, :e_conv, :e_chem, :e_src, :e_res); d[nm] = Z(nm); end
    if TX.TDC_ON[]       # M3b: w, D_mix, the exchange Q, a, U, the residual F of the w row (interfaces); P_t/P and the gas-energy source (zones)
        for nm in (:tdc_w, :tdc_Dmix, :tdc_Q, :tdc_a, :tdc_U, :tdc_Fres); d[nm] = I(nm); end
        for nm in (:tdc_PtP, :e_tdc); d[nm] = Z(nm); end
    end
    for nm in (:nad_eq, :nad_fr, :del_eq, :del_fr, :nab_rad, :nab_mu, :nab_L, :N2_S, :N2_eq, :N2_L, :g_eff, :Ftot, :Fconv_over_Ftot, :Fpot_over_Ftot, :fpot_th, :w_th, :mach_i, :t_settle, :t_N, :t_turn, :cp_fr)
        d[nm] = NaNv(n + 1)
    end
    arad_ = 7.565723e-15; cc_ = 2.99792458e10
    for j in 2:n
        kk = a + j - 1; ρi = 0.5 * (ρ[j-1] + ρ[j]); Ti = 0.5 * (T[j-1] + T[j])
        d[:w_th][j] = d[:nad_used][j] < 1e29 ? 1.0 : 0.0
        try
            th = PM.conv_thermo(mat, kk, ρi, Ti)
            d[:nad_eq][j] = th.nad_eq; d[:nad_fr][j] = th.nad_fr; d[:del_eq][j] = th.del_eq; d[:del_fr][j] = th.del_fr; d[:cp_fr][j] = th.cp_fr
            Δr = rc[j] - rc[j-1]; lT = log(T[j] / T[j-1]) / Δr; lP = log(P[j] / P[j-1]) / Δr; lm = log(mu[j] / mu[j-1]) / Δr
            ge = -(P[j] - P[j-1]) / (ρi * Δr); d[:g_eff][j] = ge
            d[:N2_S][j] = ge * th.del_fr * (lT - th.nad_fr * lP)
            d[:N2_eq][j] = ge * th.del_eq * (lT - th.nad_eq * lP)
            d[:N2_L][j] = d[:N2_S][j] - ge * th.phi * lm
            d[:nab_mu][j] = log(mu[j] / mu[j-1]) / log(P[j] / P[j-1])
            d[:nab_L][j] = th.nad_fr + th.phi / th.del_fr * d[:nab_mu][j]
        catch; end
        # all gates open: w_tau = w_inf = 1, and where the thermodynamic gate is closed (network zone, ∇_ad,used = 1e30) the frozen ∇_ad, c_P, δ in its place
        try
            th_ = PM.conv_thermo(mat, kk, ρi, Ti); open_ = d[:w_th][j] == 1.0
            na_ = open_ ? d[:nad_used][j] : th_.nad_fr; cp_ = open_ ? d[:cp_used][j] : th_.cp_fr; dl_ = open_ ? d[:del_used][j] : th_.del_fr
            g_ = d[:g][j]; HP_ = d[:HP][j]; l_ = d[:lmix][j]; sup_ = d[:nab][j] - na_
            if sup_ > 0
                U_ = 3 * arad_ * cc_ * Ti^3 / (cp_ * ρi^2 * d[:kap][j] * l_^2) * sqrt(8HP_ / (g_ * dl_)); ζ_ = sqrt(U_^2 + sup_) - U_
                d[:fpot_th][j] = ρi * cp_ * Ti * sqrt(g_ * dl_) * l_^2 * ζ_^3 / (4sqrt(2.0) * HP_ * sqrt(HP_))
            else
                d[:fpot_th][j] = 0.0
            end
        catch; end
        d[:nab_rad][j] = 3 * d[:kap][j] * ρi * d[:frad][j] * d[:HP][j] / (4 * arad_ * cc_ * Ti^4)
        Ft = d[:frad][j] + d[:fconv][j]; d[:Ftot][j] = Ft
        d[:Fconv_over_Ftot][j] = d[:fconv][j] / Ft; d[:Fpot_over_Ftot][j] = d[:fpot_th][j] / Ft
        cs_i = 0.5 * (cs[j-1] + cs[j]); d[:mach_i][j] = abs(u[j]) / cs_i
        d[:t_settle][j] = d[:HP][j] / abs(u[j]); d[:t_N][j] = 1 / sqrt(abs(d[:N2_L][j]))
        vt = d[:vconv][j] > 0 ? d[:vconv][j] : d[:vpot][j]; d[:t_turn][j] = vt > 0 ? d[:lmix][j] / vt : NaN
    end
    # the sub-shock unstable layer
    js = findfirst(<(-2e5), u); lay = fill(NaN, 8)
    if js !== nothing
        Rst = r[js]; lay[8] = Rst / 6.957e10; jin = 0; jout = 0
        for j in js-1:-1:2
            (d[:mach_i][j] < 0.5 && (d[:N2_L][j] < 0 || d[:N2_eq][j] < 0)) || break
            jout == 0 && (jout = j); jin = j
        end
        if jout > 0
            Gc = 6.674e-8; js_ = jin:jout
            td2(j) = rc[j]^3 / (Gc * max(m_[j], 1e-30))          # t_dyn² = r³/(G m) at the interface's zone centre, m the interface mass
            mx(v) = (w = filter(!isnan, v); isempty(w) ? NaN : maximum(w))
            st = [(w = filter(!isnan, (d[:N2_L][j], d[:N2_eq][j])); isempty(w) ? NaN : -minimum(w) * td2(j)) for j in js_]
            lay .= [r[jin] / Rst, r[jout] / Rst, (m_[jout] - m_[jin]) / 1.989e33, mx(st), mx(d[:Fconv_over_Ftot][js_]), mx(d[:Fpot_over_Ftot][js_]), length(js_), Rst / 6.957e10]
        else
            lay .= [NaN, NaN, 0.0, NaN, NaN, NaN, 0.0, Rst / 6.957e10]
        end
    end
    if TX.TDC_ON[]       # integrated turbulent energy E_t = Σ m_k w_k² [erg] (interface CV masses) and P_t/P at the post-shock zone (the zone inside R*; the max over the three zones inside R*)
        w0 = [isnan(q) ? 0.0 : q for q in d[:tdc_w]]; mz = diff(m_); mk = [0.5 * ((j > 1 ? mz[j-1] : 0.0) + (j <= n ? mz[j] : 0.0)) for j in 1:n+1]
        ptp = d[:tdc_PtP]; jz = js === nothing ? 0 : js - 1
        d[:tdc_scal] = [sum(mk .* w0 .^ 2), jz >= 1 ? ptp[jz] : NaN, jz >= 1 ? maximum(filter(!isnan, ptp[max(1, jz-2):jz]); init = -Inf) : NaN]
    end
    d[:layer] = lay; d[:branch] = br; d[:mach] = mach; d[:mu] = mu; d[:t_chem] = tch
    NamedTuple{Tuple(sort(collect(keys(d))))}(Tuple(d[nm] for nm in sort(collect(keys(d)))))
end
"One snapshot into the RUNREC stream: everything needed to redraw any profile panel offline."
function recsnap!(k)
    REC.io[] === nothing && recopen!(k)
    d = [ChemNet.full(mat.x[j], mat.net) for j in a:b]
    isempty(REC.spn[]) && (REC.spn[] = sort(collect(keys(d[1]))))
    sp = REC.spn[]
    sd = PM.SLDIAG[]
    snap = (; step = k, t = TF.getscalar(:timen), dt = TF.getscalar(:dtime), a, b,
            r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1], m = 4π .* (ENVEUL === nothing ? TX.XMINT[][a:b+1] : envmenc()),
            ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], pg = TF.var(:pgn)[a:b],
            er = TF.var(:ern)[a:b], fr = TF.var(:frn)[a:b+1], chi = TF.var(:chifn)[a:b], kp = TF.var(:xkpn)[a:b],
            fedd = TF.var(:fedd)[a:b], qf = TF.var(:qf)[a:b], eq = [eqflag(j) for j in a:b], Λext = copy(mat.Λext[a:b]),
            spn = sp, Y = Float64[Float64(ForwardDiff.value(d[j][s])) for s in sp, j in eachindex(d)],
            gasq = copy(TX.GASQ[]), hotY = HOTON ? SG[].Y[:, a:b] : nothing, Lnuc = LNUC[], Llines = LLINE[],
            nchemfb = PM.NCHEMFB[], kph = sd === nothing ? 0 : sd.kph, Teff = sd === nothing ? NaN : Float64(sd.Te))
    mat isa PM.NetGas && PM.HDCHEM[] && (snap = merge(snap, (; hd = PM.hd_snap(mat, a, b))))     # HDCHEM=1: x_D+, x_HD, HD cooling
    if CONVDIAG       # diagnostic only: never touches the solution (see convdiag_snap)
        cd = try convdiag_snap(k, snap) catch e; @warn "CONVDIAG: snapshot diagnostics failed at step $k" exception = (e, catch_backtrace()); nothing end
        cd === nothing || (snap = merge(snap, (; cd)))
    end
    serialize(REC.io[], snap); flush(REC.io[])
end
# ETOTMON=1 (diagnostic, default off; no influence on the solution): the total-energy monitor of every converged implicit solve, TX.etot_monitor, via TD.STAGE_HOOK.
# `buf`: the stages since the last trace line (DTERRTRACE), `pG/pv/pz`: the largest |G|, the largest per-zone value and its zone since the last perf.csv row.
mutable struct EtotMon; on::Bool; buf::Vector{Any}; pG::Float64; pv::Float64; pz::Int; end
const EM = EtotMon(false, Any[], 0.0, 0.0, 0)
"The ETOT text of the DTERR trace line: the stages since the last line (G = global total-energy residual / energy throughput of a solve, see TX.etot_monitor)."
function etot_text()
    isempty(EM.buf) && return ""
    b_ = EM.buf; EM.buf = Any[]
    Gl = b_[end][4].G; Σi = sum(x[4].imb for x in b_); Σt = sum(x[4].thr for x in b_)
    iG = argmax([abs(x[4].G) for x in b_]); iv = argmax([x[4].vmax for x in b_]); is = argmax([x[4].vshare for x in b_])
    @sprintf(" | ETOT %d solves: G(last %s) %+.2e, sum %+.2e, max|G| %.2e (step %d %s); worst zone/own throughput z%d %.2e (step %d %s), largest share z%d %.2e (step %d %s)",
             length(b_), b_[end][2], Gl, Σi / max(Σt, 1e-300), abs(b_[iG][4].G), b_[iG][1], b_[iG][2], b_[iv][4].zmax, b_[iv][4].vmax, b_[iv][1], b_[iv][2],
             b_[is][4].zshare, b_[is][4].vshare, b_[is][1], b_[is][2])
end
"One perf.csv row (called at the TIMING report, before the accumulators are reset)."
function recperf!(k, tnow)
    REC.io[] === nothing && recopen!(k)
    gct = Base.gc_num().total_time; gcb = Base.gc_bytes()
    st = Trbdf2.STATS; un_ = TF.var(:un)[a:b+1]
    oth = time() - tnow + TOTH[] - TCHEM[] - TLINE[]
    ns = REC.klast[] < 0 ? 0 : k - REC.klast[]
    open(joinpath(RUNREC, "perf.csv"), "a") do io
        @printf(io, "%d,%.8e,%.4e,%.1f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%.6e,%.6g,%.6g,%.3f,%s,%.3e,%.3e,%d\n",
                k, TF.getscalar(:timen) / 3.156e7, TF.getscalar(:dtime), time() - REC.wall0, getpid(), ns, TSTEP[], TCHEM[], TLINE[], oth,
                (gct - REC.gct[]) / 1e9, (gcb - REC.gcb[]) / 2^30, TD.NDIV[], TD.NREP[], get(st, :steps, 0), get(st, :fail_tr, 0), get(st, :fail_bdf, 0),
                get(st, :be_fallback, 0), NALICALL[], PM.NCHEMFB[], Threads.nthreads(), b - a + 1, log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a],
                REC.mstar[], REC.rstar[], minimum(un_) / 1e5, REC.dtlim[], EM.pG, EM.pv, EM.pz)
    end
    EM.pG = 0.0; EM.pv = 0.0; EM.pz = 0
    REC.klast[] = k; REC.gct[] = gct; REC.gcb[] = gcb
end
# H2 line transfer: ALI only (statistical equilibrium of the levels coupled to the formal solution, H2Lines.cooling_nlte; the
# escape-probability and hybrid modes are gone). Every transfer call starts from the previous call's converged populations (POPS0,
# aligned from the outside on the fixed mesh; cold start at the first call and after a remap).
const LEDGERZ = Ref(false)
const NALICALL = Ref(0)
const POPS0 = Ref(zeros(0, 0))
"The previous call's populations on a reduced grid of nz_ zones (columns that have no predecessor are NaN: cold start), or nothing."
function pops0_for(nz_)
    P = POPS0[]; n0 = size(P, 2)
    n0 == 0 && return nothing
    n0 == nz_ && return P
    Q = fill(NaN, size(P, 1), nz_)
    for j in 1:nz_; jo = j + n0 - nz_; jo >= 1 && (Q[:, j] .= view(P, :, jo)); end
    Q
end
function lineinc()
    LINEINC in ("1", "cmb") || return nothing
    Tb = parse(Float64, get(ENV, "TCMB", "0")); sd = PM.SLDIAG[]
    (LINEINC == "cmb" || sd === nothing || !(sd.Te > 0)) && return (Tbg = Tb, Tstar = 0.0, Rph = 0.0)
    (Tbg = Tb, Tstar = Float64(sd.Te), Rph = Float64(TF.var(:rn)[sd.kph+1]))
end
# LINEIMP (always on since the clean-up; the switch is gone): the line cooling of each zone is a function of the zone's state in the chemistry and in
# the hydro Newton (PopIIIMaterial.LINEIMP: populations re-solved at the iterate's T with the radiative rates β, J̄ of the last transfer
# call), not Λ held from the last call. Retires as devices: the IMPCOOL slope (second transfer call) and its [−10, 30] clamp, the
# Δt ≤ 0.1 (or DTCOOL) t_cool cap, TCOOLTH (only the cap used it). LINEFRAC still skips transfer calls (the model is then kept).
const LINEIMP = true
const LIJCMB = Ref{Vector{Float64}}(Float64[])
function nanscan(where, k)
    n = 0
    for nm in (:tn, :dn, :egn, :ern, :un, :frn, :rn, :pgn, :xkpn, :chifn)
        v = TF.var(nm)
        for kk in a:b+1
            isfinite(v[kk]) && continue
            n += 1; n <= 12 && pr("NANSCAN step %d %s: %s[%d] (zone %d) = %s; T %.4g rho %.4g eq %s x %s\n", k, where, String(nm), kk, kk - a + 1, string(v[kk]),
                                  TF.var(:tn)[min(kk, b)], TF.var(:dn)[min(kk, b)], string(mat.eqzone[min(kk, b)]), string(round.(mat.x[min(kk, b)], sigdigits = 3)))
        end
    end
    for kk in a:b
        all(isfinite, mat.x[kk]) && continue
        n += 1; n <= 12 && pr("NANSCAN step %d %s: composition of zone %d not finite: %s (T %.4g)\n", k, where, kk - a + 1, string(mat.x[kk]), TF.var(:tn)[kk])
    end
    nothing
end
"(LINEIMP) zones remapped: drop the line models (the held Λext of the remap is used until the next transfer call)."
li_inval!() = (PM.li_invalidate!(); POPS0[] = zeros(0, 0); nothing)
"(LINEIMP) per TITAN zone, the line model from the transfer call `Lc` (H2Lines.LAST[]) on the reduced grid (see LINEIMP above)."
function li_store!(Lc, dd, nh2, nl, zk, ncz, nz_)
    hasproperty(Lc, :β) || return nothing
    nk = length(mat.x); length(PM.LIMODEL[]) == nk || (PM.LIMODEL[] = Any[nothing for _ in 1:nk])
    for kk in a:b; PM.LIMODEL[][kk] = nothing; end
    if isempty(LIJCMB[])
        Tb = parse(Float64, get(ENV, "TCMB", "0"))
        LIJCMB[] = [Tb > 0 ? H2Lines.Bν(ν_, Tb) : 0.0 for ν_ in HLL.ν]
    end
    jc_ = LIJCMB[]
    for j in ncz+1:nz_
        β_ = Lc.β[:, j]; J_ = max.(Lc.Jeff[:, j], β_ .* jc_)
        PM.li_set!(zk(j), HLV, β_, J_, nl.nH[j], nh2[j], nl.nHe[j], nl.nHp[j], dd[j], Lc.Labs[j])
    end
    nothing
end
"""
S6: H₂ line transfer outside the continuum photosphere only (RAD runs, LTAUCUT > 0 [10]): the zones deeper than continuum
optical depth LTAUCUT from the surface are merged into one transparent central zone without H₂ for the line transfer (their
line cooling is 0: the continuum carries the energy there and H₂ is dissociated); the cost ∝ zones² × bins drops (the hot core
also set the velocity span of the frequency grid). Returns (Λ, Lr) on the full grid.
"""
function linecool(rr, uu, dd, tt, nh2, k)
    jc = 0
    if RAD && LTAUCUT > 0
        ch = TF.var(:chifn); τ = 0.0
        for j in nsh:-1:1
            τ += ch[a+j-1] * dd[j] * (rr[j+1] - rr[j])
            τ > LTAUCUT && (jc = j; break)
        end
        jc = jc >= nsh - 4 ? 0 : jc                  # (keep at least a few zones)
    end
    # continuum extinction at the line frequencies, grey: χ_c = κ_P ρ per zone (TITAN's xkpn, the G2Opacity Planck mean: CIA, H⁻ bf/ff, H ff/bf; not frequency dependent);
    # the merged core is opaque to the rays (χ_c Δr = LTAUCUT across its radius)
    kpv = TF.var(:xkpn)
    jc < 2 && return linecool_full(rr, uu, dd, tt, nh2, k, 0, 0, [kpv[a+j-1] * dd[j] for j in 1:nsh])
    # the merged core is NCORE transparent zones with interfaces uniform in p² (LINECORE, S6): impact parameters through the core, so the
    # inner nodes have rays over all μ (with one core zone the inner edge saw a single ray through the centre); u linear in r inside
    nc = LINECORE; rin = rr[jc+1]
    rcz = [rin * sqrt((i - 1) / nc) for i in 1:nc]
    r2 = [rcz; rr[jc+1:end]]; u2 = [uu[jc+1] .* rcz ./ rin; uu[jc+1:end]]
    m_c = sum(dd[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:jc); d2 = [fill(m_c / (4π / 3 * rin^3), nc); dd[jc+1:end]]
    t2 = [fill(tt[jc], nc); tt[jc+1:end]]; h2 = [zeros(nc); nh2[jc+1:end]]
    chic = [fill(LTAUCUT / rin, nc); [kpv[a+j-1] * dd[j] for j in jc+1:nsh]]
    Λ2, L2, x2 = linecool_full(r2, u2, d2, t2, h2, k, jc, nc, chic)
    Λ = zeros(nsh); Λ[jc+1:end] .= Λ2[nc+1:end]; Lr = zeros(nsh + 1); Lr[jc+1:end] .= L2[nc+1:end]
    # the line power the opaque merged core absorbs (the continuum there is thick: LTAUCUT) heats the gas of the zones it stands for, uniformly in volume
    Vc = [4π / 3 * (r2[j+1]^3 - r2[j]^3) for j in 1:nc]; Λ[1:jc] .= -sum(H2Lines.LAST[].Labs[1:nc] .* Vc) / (4π / 3 * rin^3)
    Λ, Lr, x2
end
"""
Line cooling per zone (and Lr) for zone temperatures tt (ALI line transfer, H2Lines.cooling_nlte, with the grey continuum absorption chic);
step k; jc = number of core zones merged away. Λ is the net divergence of the line flux −∇·F: collisional line cooling minus the line power
the continuum absorbed, which the gas takes back (PopIIIMaterial.li_cool does the same for the zone's cooling at the iterate).
"""
function linecool_full(rr, uu, dd, tt, nh2, k, jc, ncz = jc > 0 ? 1 : 0, chic = nothing)
    nz_ = length(dd); zk(j) = j <= ncz ? a + jc - 1 : a + (j - ncz) + jc - 1   # TITAN zone of reduced zone j (ncz merged-core zones first)
    yc = [ycomp(zk(j)) for j in 1:nz_]; nHn = dd ./ PC.μH
    nl = (LV = HLV, nH = [PC.xH(yc[j]) * nHn[j] for j in 1:nz_], nHe = PC.yHe .* nHn, nHp = [yc[j][1] * nHn[j] for j in 1:nz_])
    tq = time()
    out_ = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, inc = lineinc(), chic, pops0 = k >= 0 ? pops0_for(nz_) : nothing)
    Lc = H2Lines.LAST[]
    if k >= 0
        POPS0[] = copy(Lc.pops); NALICALL[] += 1
        li_store!(Lc, dd, nh2, nl, zk, ncz, nz_)
    end
    if TIMING && k >= 0
        # closed line-energy ledger: collisional line cooling of the levels = line luminosity leaving the surface + line power absorbed by the continuum
        vol = [4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_]
        Λcl = [H2Levels.collcool(HLV, tt[j], view(Lc.pops, :, j), nl.nH[j], nh2[j], nl.nHe[j], nl.nHp[j]) for j in ncz+1:nz_]
        Sc = sum(Λcl .* vol[ncz+1:nz_]; init = 0.0); Sa = sum(Lc.Labs .* vol); Sk = sum(Lc.Labs[1:ncz] .* vol[1:ncz]; init = 0.0); Le = out_[2][end]
        pr("        ALI call %.1f s (ALI iterations %d, converged %s, warm start %s)\n", time() - tq, Lc.iters, Lc.converged, pops0_for(nz_) === nothing ? "no" : "yes")
        pr("        LEDGER k %d: sum V Lambda_coll %.6e = L_esc %.6e + absorbed by continuum %.6e (merged core %.3e); residual %.3e (relative to sum V Lambda_coll %.2e)\n",
           k, Sc, Le, Sa, Sk, Sc - Le - Sa, (Sc - Le - Sa) / max(abs(Sc), 1e-300))
        if !LEDGERZ[]               # (once) the zones with continuum depth from the surface 0.1 .. 10 at the first call: collisional cooling, continuum-absorbed power, net
            LEDGERZ[] = true; τz = 0.0; kpz = TF.var(:xkpn)
            for j in nz_:-1:ncz+1
                τz += (chic === nothing ? 0.0 : chic[j]) * (rr[j+1] - rr[j])
                0.1 <= τz <= 10 && pr("        LEDGERZ zone %d (TITAN %d) r %.3e R_sun tau_cont(surface) %.3f T %.0f nH2 %.2e: Lambda_coll %.4e Lambda_cont %.4e (%.2e of coll) net %.4e flux-based %.4e\n",
                                      j, zk(j), rr[j] / 6.957e10, τz, tt[j], nh2[j], Λcl[j-ncz], Lc.Labs[j], Lc.Labs[j] / Λcl[j-ncz], Λcl[j-ncz] - Lc.Labs[j], out_[1][j])
            end
        end
    end
    out_
end
const LLINE = Ref(0.0)
const JCMIN = Ref(1)
const KOFF = Ref(0)
const SG = Ref{Any}(nothing)
# LINEIMP line models as plain data (the level data LV is the constant HLV, the PopWS workspace is rebuilt): a checkpoint stays readable without the modules
lim_plain(m) = m === nothing ? nothing : (; β = copy(m.β), J = copy(m.J), nH = m.nH, nH2 = m.nH2, nHe = m.nHe, nHp = m.nHp, ρ = m.ρ, hc = m.hc)
lim_make(p) = p === nothing ? nothing : (; LV = HLV, β = copy(p.β), J = copy(p.J), nH = p.nH, nH2 = p.nH2, nHe = p.nHe, nHp = p.nHp, ρ = p.ρ,
                                         hc = hasproperty(p, :hc) ? p.hc : 0.0, ws = H2Levels.PopWS(HLV))
# (a checkpoint from before the continuum return: hc = 0 until the next call)
matstate() = mat isa PM.NetGas ? (; x = deepcopy(mat.x), eqzone = copy(mat.eqzone), Λext = copy(mat.Λext), Tref = copy(mat.Tref), slope = copy(mat.slope), Y = HOTON ? copy(SG[].Y) : nothing,
                                     star = HOTON ? (; fscale = copy(SG[].fscale), fslope = copy(SG[].fslope), fcur = copy(SG[].fcur), Dmix = copy(SG[].Dmix), dtcap = SG[].dtcap, capE = SG[].capE, ledger = SG[].ledger, Yprev = copy(SG[].Yprev)) : nothing,
                                     # lagged over one step (set in chemistry, used by the next hydro step): without them a restart is not exact
                                     lag = (; NLOPC = deepcopy(PM.NLOPC[]), GASQ = copy(TX.GASQ[]), RAYQ = copy(TX.RAYQ[]), CQB = copy(TX.CQB[]), YQB = copy(TX.YQB[]), SLDIAG = deepcopy(PM.SLDIAG[]), CHEMB = copy(TX.CHEMB[]),
                                                 LIMODEL = [lim_plain(m_) for m_ in PM.LIMODEL[]], NLOST = deepcopy(PM.NLOST[]))) :   # (lag.LIMODEL: LINEIMP's H2 line models, set by the last transfer call, used by the next hydro step)
             error("checkpoint: only the NetGas material (CHEM=net*) is supported")
# jinc: 1 for a checkpoint written inside the step loop's output callback, where TITAN's step counter jstep has not yet been advanced for the
# step just done (run_julia! increments it after output()); the restart then continues with jstep + 1, as the uninterrupted run would
# HELDDUMP=1 (diagnostic, default off): with every checkpoint also write <name>.held, the numeric content of every global (arrays, Refs, tuples, dicts) of the
# solver modules and of this driver, to find the Julia-side state that a restart does not reproduce (tools/ckheld.jl compares two such files)
function held_flat!(out, name, v, depth = 0)
    depth > 6 && return
    if v isa Number; v isa Bool || v isa Integer || v isa AbstractFloat ? push!(out, name => Float64[v]) : nothing
    elseif v isa AbstractArray{<:Number}; length(v) <= 20_000_000 && push!(out, name => Float64.(vec(collect(v))))
    elseif v isa Ref; isassigned(v) && held_flat!(out, name, v[], depth + 1)
    elseif v isa AbstractArray; length(v) <= 200_000 && for (i, e) in enumerate(v); held_flat!(out, "$name[$i]", e, depth + 1); end
    elseif v isa Union{Tuple,NamedTuple}; for (i, e) in enumerate(v); held_flat!(out, "$name.$i", e, depth + 1); end
    elseif v isa AbstractDict; length(v) <= 10_000 && for kk in sort(collect(keys(v)); by = string); held_flat!(out, "$name{$kk}", v[kk], depth + 1); end
    end
    nothing
end
function held_dump(path)
    out = Pair{String,Vector{Float64}}[]
    mods = Pair{String,Module}["Main" => @__MODULE__]      # the driver's own globals (were Main's) and every module that was include()d into Main (not the registered packages)
    for nm in names(Main; all = true)
        isdefined(Main, nm) || continue; v = getfield(Main, nm)
        v isa Module && v !== Main && v !== Base && v !== Core && pathof(v) === nothing && push!(mods, string(nm) => v)
    end
    for (mn, m) in mods
        for nm in names(m; all = true)
            s_ = string(nm); (startswith(s_, "#") || nm in (:ans, :ARGS, :Core, :Base, :Main, :ENV, :TF)) && continue
            isdefined(m, nm) && !isconst(m, nm) && m === (@__MODULE__) && continue
            v = try getfield(m, nm) catch; continue end
            (v isa Module || v isa Function || v isa Type) && continue
            try held_flat!(out, "$mn.$s_", v) catch; end
        end
    end
    serialize(path, Dict(out))
end
function checkpoint(name, k; jinc = 0)
    get(ENV, "HELDDUMP", "0") == "1" && held_dump(joinpath(CKDIR, "m2c_$(name).held"))
    mkpath(CKDIR)
    TF.save_state(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    ax = (; k, jnext = Int(TF.getscalar(:jstep)) + jinc, mat = matstate(), POPS0 = copy(POPS0[]),
        NALICALL = NALICALL[], NIMP = NIMP[], LLINE = LLINE[], nextln = nextln[], SNAPNEXT = SNAPNEXT[], SNAPS,
        trb = (; STATS = copy(Trbdf2.STATS), BE_LEFT = Trbdf2.BE_LEFT[]),
        # solver / controller / held-input state that is not in the TITAN COMMON (without it the first steps after a restart differ: no FSAL error estimate,
        # no PI history, no mesh feature tracking, a restarting table-to-network ramp)
        ctl = (; HNEXT = Trbdf2.HNEXT[], ERRPREV = Trbdf2.ERRPREV[], HEUR = Trbdf2.HEUR[], LASTACC = Trbdf2.LASTACC[], DSAL = deepcopy(Trbdf2.DSAL[]),
               ESTATS = copy(Trbdf2.ESTATS), LIMHIST = copy(Trbdf2.LIMHIST)),
        edd = (; HEDD = copy(TD.HEDD), HXEDD = copy(TD.HXEDD), FXEDD = copy(TD.FXEDD), TIEW = copy(TD.TIEW), HPLUS = copy(TC.HPLUS), GXEDD = TD.GXEDD[]),   # Eddington closure arrays held between steps
        geq = (; J0 = GEQJ0[], PREV = deepcopy(GEQPREV[]), HELD = deepcopy(TX.GRIDEQ[])), pswramp = (; J0 = PM.RAMPJ0[], K0 = PM.RAMPK0[]), PEXT = TX.PEXT[], env = Dict(k_ => v_ for (k_, v_) in ENV if k_ in CKENV))
    TX.MG[] && (ax = merge(ax, (; mg = (; edges = copy(TX.MGEDGE[]), E = deepcopy(TX.MGEN[]), F = deepcopy(TX.MGFN[])))))   # MGROUPS: the group fields
    TX.TDC_ON[] && (ax = merge(ax, (; tdc = (; WTN = copy(TX.WTN[]), WTO = copy(TX.WTO[]), a = Int(TF.getscalar(:ngrs)), b = Int(TF.getscalar(:ngre))))))   # TDC: the turbulent velocity
    serialize(joinpath(CKDIR, "m2c_$(name).aux"), ax)
end
const GEQHELD = Ref{Any}(nothing)
const JNEXT = Ref(-1)
const CKAUX = Ref{Any}(nothing)
function edd_restore!(e)
    for (dst, src) in ((TD.HEDD, e.HEDD), (TD.HXEDD, e.HXEDD), (TD.FXEDD, e.FXEDD), (TD.TIEW, e.TIEW), (TC.HPLUS, e.HPLUS))
        resize!(dst, length(src)); copyto!(dst, src)
    end
    TD.GXEDD[] = e.GXEDD
end
const MGAUX = Ref{Any}(nothing)
# (MGROUPS) the group fields of the checkpoint a run restarted from
# settings that must agree between a run and its restart (the set-up and the physics); the rest may change
const CKENV = ("IC", "TCMB", "AVZONE", "AVQL1", "AVCQ2", "CHEM", "NBIN", "K5", "MCLOUD", "LRJOIN", "DLNJOIN", "DM", "PSWITCH",
               "INTEGRATOR", "THETA", "XE0", "RAD")
# (not CDYN: the accretion phase relaxes the dynamical Δt limit at a restart)
function restart!(name; allow = ())
    aux = deserialize(joinpath(CKDIR, "m2c_$(name).aux"))
    for (k_, v_) in aux.env
        (k_ in allow || !(k_ in CKENV)) && continue        # (a key whose switch was removed in the clean-up: its value is the method now)
        get(ENV, k_, nothing) == v_ || error("restart: ENV $k_ = $(repr(get(ENV, k_, nothing))) differs from the checkpoint's $(repr(v_))")
    end
    for k_ in CKENV; k_ in allow && continue; haskey(ENV, k_) && !haskey(aux.env, k_) && error("restart: ENV $k_ is set but was not in the checkpoint's run"); end
    TF.load_state!(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1   # (zone count may differ after rezoning)
    s = aux.mat; mat.x .= deepcopy(s.x); mat.eqzone .= s.eqzone; mat.Λext .= s.Λext; mat.Tref .= s.Tref; mat.slope .= s.slope
    ChemNet.HECHEM && mat isa PM.NetGas && PM.he_pad!(mat)                    # (restart from a checkpoint written without HECHEM)
    HOTON && hasproperty(s, :Y) && s.Y !== nothing && (SG[].Y .= s.Y)      # (a checkpoint from before HOT=1 keeps the initial abundances)
    if HOTON && hasproperty(s, :star) && s.star !== nothing
        t_ = s.star; g_ = SG[]; g_.fscale .= t_.fscale; g_.fslope .= t_.fslope; g_.fcur .= t_.fcur; g_.Dmix = t_.Dmix; g_.dtcap = t_.dtcap; g_.capE = t_.capE; g_.ledger = t_.ledger; g_.Yprev .= t_.Yprev
    end
    if hasproperty(s, :lag)                        # (older checkpoints: these start empty, so the first step after a restart differs)
        l_ = s.lag; PM.NLOPC[] = deepcopy(l_.NLOPC); PM.SLDIAG[] = deepcopy(l_.SLDIAG)
        PM.NLOIT && PM.nlost_restart!(mat)       # (older checkpoints: only the zones with non-LTE factors (NLOPC) get their composition and three-level parameters back)
        # NLOST (composition and three-level parameters at the end of the last chemistry = the populations of the next step's iterate-based opacity): the
        # rebuild above misses every zone when NLOPC is empty (NLOIT runs), so the first step after a restart had the plain table opacity in all network zones
        hasproperty(l_, :NLOST) && (PM.NLOST[] = Any[deepcopy(e_) for e_ in l_.NLOST])
        length(l_.GASQ) == length(TX.GASQ[]) && (TX.GASQ[] .= l_.GASQ); length(l_.RAYQ) == length(TX.RAYQ[]) && (TX.RAYQ[] .= l_.RAYQ)
        hasproperty(l_, :CQB) && length(l_.CQB) == length(TX.CQB[]) && (TX.CQB[] .= l_.CQB; TX.YQB[] .= l_.YQB)
        if !hasproperty(l_, :CQB) && !isempty(TX.CQB[]) && PM.SLDIAG[] !== nothing      # checkpoint from before CQB/YQB: cq, yq were booked inside RAYQ/GASQ; separate them (SLDIAG keeps them per zone)
            nconv = 0
            for (k_, t_) in enumerate(PM.SLDIAG[].tl3)
                t_ === nothing && continue
                TX.CQB[][k_] = t_.chemq; TX.YQB[][k_] = t_.cycq; TX.RAYQ[][k_] -= t_.chemq; TX.GASQ[][k_] -= t_.cycq; nconv += 1
            end
            pr("restart: checkpoint without CQB/YQB; chemistry and cycle heat of %d zones taken out of RAYQ/GASQ and booked as CQB/YQB\n", nconv)
        end
        if hasproperty(l_, :LIMODEL)
            PM.LIMODEL[] = Any[lim_make(m_) for m_ in l_.LIMODEL]
        else
            pr("WARNING: checkpoint %s has no H2 line models (older code): the first step after the restart has no line cooling until the first transfer call\n", name)
        end
        hasproperty(l_, :CHEMB) && length(l_.CHEMB) == length(TX.CHEMB[]) && (TX.CHEMB[] .= l_.CHEMB)
    end
    POPS0[] = hasproperty(aux, :POPS0) ? copy(aux.POPS0) : zeros(0, 0); NALICALL[] = aux.NALICALL; NIMP[] = aux.NIMP
    LLINE[] = aux.LLINE; nextln[] = aux.nextln; SNAPNEXT[] = aux.SNAPNEXT
    SNAPS === nothing || aux.SNAPS === nothing || (empty!(SNAPS); append!(SNAPS, aux.SNAPS))
    empty!(Trbdf2.STATS); merge!(Trbdf2.STATS, aux.trb.STATS); Trbdf2.BE_LEFT[] = aux.trb.BE_LEFT
    if hasproperty(aux, :ctl)
        c_ = aux.ctl; Trbdf2.HNEXT[] = c_.HNEXT; Trbdf2.ERRPREV[] = c_.ERRPREV; Trbdf2.HEUR[] = c_.HEUR; Trbdf2.LASTACC[] = c_.LASTACC; Trbdf2.DSAL[] = deepcopy(c_.DSAL)
        empty!(Trbdf2.ESTATS); merge!(Trbdf2.ESTATS, c_.ESTATS); empty!(Trbdf2.LIMHIST); merge!(Trbdf2.LIMHIST, c_.LIMHIST)
    else
        pr("WARNING: checkpoint %s has no solver-controller state (written by older code): the first step after the restart has no TR-BDF2 error estimate and the controller history is empty\n", name)
    end
    if hasproperty(aux, :geq)
        GEQJ0[] = aux.geq.J0; GEQPREV[] = deepcopy(aux.geq.PREV)
        # the grid equation's held monitor (tauv, sx) was computed in the checkpointed step before chem!'s EOS/opacity refresh: not reproducible from the final state
        GEQHELD[] = aux.geq.HELD
         PM.RAMPJ0[] = aux.pswramp.J0; PM.RAMPK0[] = aux.pswramp.K0
    else
        pr("WARNING: checkpoint %s has no grid-equation / PSWRAMP state (older code): the mesh feature tracking restarts empty and the ramps restart\n", name)
    end
    TX.PEXT[] = aux.PEXT; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    hasproperty(aux, :mg) && (MGAUX[] = aux.mg)
    TF.setscalar!(:tmax, 1e3tff)                   # a run stopped by TSTOP/NSTOP has tmax = its last time
    CKAUX[] = aux
    JNEXT[] = hasproperty(aux, :jnext) ? aux.jnext : -1     # (set after the grideq_step! below; older checkpoints: jstep stays that of the last step, one behind the uninterrupted run)
    KOFF[] = aux.k
    pr("restarted from checkpoint %s at step %d, t = %.6e yr, log n_c = %.3f, T_c = %.1f K\n", name, aux.k, TF.getscalar(:timen) / 3.156e7,
       log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a])
end
struct ThTag end
function thchain(f, ρ, T, k)
    (ρ isa ForwardDiff.Dual || T isa ForwardDiff.Dual) || return f(ρ, T, k)
    Dout = typeof(ρ + T)
    (ForwardDiff.valtype(Dout) === Float64 && (!(ρ isa ForwardDiff.Dual) || typeof(ρ) === Dout) && (!(T isa ForwardDiff.Dual) || typeof(T) === Dout)) || return f(ρ, T, k)   # nested / mixed tags: direct
    ρv, Tv = ForwardDiff.value(ρ), ForwardDiff.value(T)
    D2 = ForwardDiff.Dual{typeof(ForwardDiff.Tag(ThTag(), Float64)),Float64,2}
    o = f(D2(ρv, ForwardDiff.Partials((1.0, 0.0))), D2(Tv, ForwardDiff.Partials((0.0, 1.0))), k)
    pρ = ρ isa ForwardDiff.Dual ? ForwardDiff.partials(ρ) : zero(ForwardDiff.partials(T))
    pT = T isa ForwardDiff.Dual ? ForwardDiff.partials(T) : zero(ForwardDiff.partials(ρ))
    map(x -> Dout(ForwardDiff.value(x), ForwardDiff.partials(x, 1) * pρ + ForwardDiff.partials(x, 2) * pT), o)
end
"""
S6-B: settings of a radiating run (after a RAD rebuild and after restarting a RAD checkpoint): StellarOpacity continuum for
NetGas, CMB intensity at the outer boundary, Eddington-factor mode, thin-envelope E_r out of the Δt rule, Δt cap.
"""
const XMEXT = Ref{Any}(nothing)
function radsettings!()
    X = 1 / (1 + 4PC.yHe)
    # OPACSENS="fcond,fext" (sensitivity test, default 1,1): κ_cond × fcond, and κ_rad × fext beyond OPAL's log R = 1 (over one dex in log R)
    os_ = parse.(Float64, split(get(ENV, "OPACSENS", "1,1"), ","))
    PM.OPAC[] = all(==(1.0), os_) ? ((ρ, T) -> (StellarOpacity.kappa_R(ρ, T, X), StellarOpacity.kappa_P(ρ, T, X))) : ((ρ, T) -> begin
        κ = StellarOpacity.kappa_R(ρ, T, X); κr = StellarOpacity.kappa_R_rad(ρ, T, X); ic = max(1 / κ - 1 / κr, zero(κ))     # 1/κ_cond,eff
        lR = log10(ρ) - 3 * (log10(T) - 6); w = clamp(lR - 1, 0.0, 1.0); w = w * w * (3 - 2w)
        (1 / (1 / (κr * os_[2]^w) + ic / os_[1]), StellarOpacity.kappa_P(ρ, T, X))
    end)
    Tb = parse(Float64, get(ENV, "TCMB", "0"))
    TF.setscalar!(:ximr, Z0.σSB * Tb^4 / π)
    TF.setscalar!(:ltran, parse(Int, get(ENV, "LTRAN", "1")))      # 1: formal solution for f, g; 2: Eddington approximation f = 1/3, g = 1/2
    TC.ER_TAUMIN[] = parse(Float64, get(ENV, "ERTAUMIN", "1.0"))   # E_r of the optically thin envelope does not limit Δt
    if HOTON                                          # convection (hot branch only: PM.material_thermo), as the S5 stars
        TX.CONVECTION[] = true; TX.ALPHA_MLT[] = parse(Float64, get(ENV, "ALPHAMLT", "1.5")); TX.MLT_HSDISC[] = true; TX.MLT_RCAP[] = true
        PM.COLDMLT[] = parse(Float64, get(ENV, "COLDMLT", "0")); PM.COLDMLT[] > 0 && (TX.MLT_NOINFALL[] = true)
        # MLTBUOY=1 (M2, provisional, default off = bit-identical; docs/POPIII_CONVECTION.md M2): MLT superadiabaticity from the blob's density contrast with the zones' thermodynamics at the iterate (TX.MLT_BUOY)
        TX.MLT_BUOY[] = get(ENV, "MLTBUOY", "0") == "1"
        TX.MLT_BUOY[] && atexit(() -> println(TX.buoy_summary()))
    end
    TD.MSGIO[] = get(ENV, "MSG", "0") == "1" ? stdout : devnull
end
"""
Rebuild the TITAN state from a zone profile (interfaces rI, u, Fr; zones ρ, T, Er; `matz` = per-zone material state, NamedTuple
of vectors x, eqzone, Λext, Tref, slope) with PS.build! (rad = true), keeping time, Δt, step number, the solver controls and
PEXT. Used by the RAD rebuild (S6-B) and by the rezoning (S6-C2).
"""
function rebuild!(rI, ρz, Tz, uI, Er, Fr, matz)
    t0, dt0, js = TF.getscalar(:timen), TF.getscalar(:dtime), TF.getscalar(:jstep)
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    pext = TX.PEXT[]
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lgrid = 3, inner = :wall, outer = :free, pext = 0.0,
              selfgrav = true, material = mat, av, xmext = XMEXT[])
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    for j in 1:nsh
        k = a + j - 1
        mat.x[k] = copy(matz.x[j]); mat.eqzone[k] = matz.eqzone[j]; mat.Λext[k] = matz.Λext[j]; mat.Tref[k] = matz.Tref[j]; mat.slope[k] = matz.slope[j]
    end
    if HOTON
        g_ = SG[]
        for j in 1:nsh; g_.Y[:, a+j-1] .= matz.Y[j]; end
        SM.ghosts!(g_); g_.Yprev .= g_.Y; g_.acc.nstage = 0
        fill!(g_.fscale, 1.0); fill!(g_.fslope, 0.0); fill!(g_.fcur, 1.0); g_.Dmix = zeros(length(TF.var(:dn)))
        if hasproperty(matz, :fslope)                 # (rezoning: level-factor slopes and the mixed regions follow the zones)
            for j in 1:nsh; g_.fscale[a+j-1] = matz.fscale[j]; g_.fslope[a+j-1] = matz.fslope[j]; g_.Dmix[a+j-1] = matz.Dmix[j]; end
        end
    end
    for k in 1:a-1; mat.x[k] = copy(mat.x[a]); mat.eqzone[k] = mat.eqzone[a]; end
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end; LINEIMP && li_inval!()
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    TX.TDC_ON[] && (fill!(TX.WTN[], 0.0); fill!(TX.WTO[], 0.0))      # TDC: a rebuilt mesh (RAD rebuild, rezoning) starts from w = 0
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    PS.install_julia!(BlockSolver)
    TX.PEXT[] = pext
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, t0); TF.setscalar!(:timeo, t0); TF.setscalar!(:dtime, dt0); TF.setscalar!(:jstep, js); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    POPS0[] = zeros(0, 0)                             # line populations: cold start at the next line-transfer call
end
"""
S6 validation: continue an S5 star (s5_run.jl checkpoint `STARFROM`, directory S5DIR) on the NetGas hot path: same profile, abundances, solver
controls and time; then T_c, L, X_c(t) can be compared with S5's own continuation.
"""
function starfrom!(name)
    dir = get(ENV, "S5DIR", "/tmp/s5"); TF.load_state!(joinpath(dir, "s5_$(name).tstate.gz"))
    aux = deserialize(joinpath(dir, "s5_$(name).aux"))
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre))
    rI = TF.var(:rn)[a:b+1]; ρz = TF.var(:dn)[a:b]; Tz = TF.var(:tn)[a:b]; uI = TF.var(:un)[a:b+1]; Er = TF.var(:ern)[a:b]; Fr = TF.var(:frn)[a:b+1]
    n = b - a + 1; XMEXT[] = TF.getscalar(:xmext); TX.PEXT[] = 0.0; aold = a
    matz = (; x = [PM.x_from_table(mat, ρz[j], min(Tz[j], 1e6)) for j in 1:n], eqzone = trues(n), Λext = zeros(n), Tref = ones(n), slope = zeros(n),
              Y = [aux.Y[:, aold+j-1] for j in 1:n], fscale = ones(n), fslope = [aux.fslope[aold+j-1] for j in 1:n], Dmix = [aux.Dmix[aold+j-1] for j in 1:n])
    rebuild!(rI, ρz, Tz, uI, Er, Fr, matz)
    SG[].dtcap = aux.dtcap; SG[].capE = aux.capE; TF.setscalar!(:tmax, 1e30)
    pr("STARFROM %s: %d zones, t = %.6e yr, T_c = %.5e K, rho_c = %.3f, X_c = %.5f\n", name, nsh, TF.getscalar(:timen) / 3.156e7, TF.var(:tn)[a], TF.var(:dn)[a], SG[].Y[SM.IH, a])
end
const Rsun = 6.957e10
"Grey optical depth from the outer edge at the interfaces i0..b+1 (the τ term of the GRIDEQ monitor)."
function gm_tau(i0)
    ch, dd, rr = TF.var(:chifn), TF.var(:dn), TF.var(:rn); t = zeros(b + 2 - i0); τ = 0.0
    for k in b:-1:i0; τ += ch[k] * dd[k] * (rr[k+1] - rr[k]); t[k-i0+1] = τ; end
    t
end
xmi(k) = ENVEUL === nothing ? TX.XMINT[][k] : TF.getscalar(:xmtot) - TF.var(:xmen)[k]
const XMEPREV = Ref(Float64[])
const ENVMS0 = Ref(0.0)
const ENVMDOT = Ref(0.0)
envmenc() = TF.getscalar(:xmtot) .- TF.var(:xmen)[a:b+1]
# interior mass per steradian at the interfaces
function envsurface()
    j = findfirst(<(-2e5), TF.var(:un)[a:b+1]); j === nothing ? nsh + 1 : j
end
function envsink!()
    ENVCENTRE && return
    xme = TF.var(:xmen); j = envsurface(); Ms = xme[a] - xme[a+j-1]
    Md = max(0.0, Ms - ENVMS0[]) / ENVSINKTAU; ENVMDOT[] = Md                  # per steradian
    φ0 = -Md / TF.var(:rn)[a]^2
    TF.setscalar!(:phil0, φ0); TF.setscalar!(:phil2, φ0 * (TF.var(:egn)[a] + TF.var(:pgn)[a] / TF.var(:dn)[a]))
end
const NXADV1 = Ref(0)
"Second-order correction to the implicit upwind result `xnew` (vectors per zone, from the start-of-step `X`): see the XADV2 note."
function xadv2_correct!(xnew, X, F, M, n)
    cadv = TF.getscalar(:cadv); nf = length(X[1])
    S = [zeros(nf) for _ in 1:n]
    for j in 2:n-1, i in 1:nf
        da = X[j][i] - X[j-1][i]; db = X[j+1][i] - X[j][i]
        S[j][i] = da * db > 0 ? cadv * da * db / (da + db) : 0.0
    end
    # antidiffusive flux through interior interface j (between zones j-1, j): F (x_face − x_donor), only where |F| < the donor's mass
    G = [zeros(nf) for _ in 1:n+1]
    for j in 2:n
        if F[j] >= 0
            F[j] < M[j-1] && (G[j] .= F[j] .* 0.5 .* S[j-1])
        else
            -F[j] < M[j] && (G[j] .= F[j] .* (-0.5) .* S[j])
        end
    end
    ok = trues(nf); xc = [xnew[j] .+ (G[j] .- G[j+1]) ./ M[j] for j in 1:n]
    for j in 1:n, i in 1:nf; xc[j][i] < 0 && (ok[i] = false); end
    for j in 1:n
        for i in 1:nf; ok[i] || (xc[j][i] = xnew[j][i]); end
        d = ChemNet.full(xc[j], mat.net)
        (d[:H] >= 0 && d[:e] >= 0) || (xc[j] = copy(xnew[j]))
    end
    all(ok) || (NXADV1[] += 1)
    for j in 1:n; xnew[j] = xc[j]; end
    nothing
end
"""
COUPLED (the method; the switch is gone): x_H⁺ and x_H₂ of every network zone are unknowns of the hydro Newton (TitanResidFast.COUPLED,
PopIIIMaterial cpl_*; implies CHEMCK). The chemistry after the step keeps only what the hydro solve does not carry (HD, the three-level
bookkeeping, the table switch); table and hot zones keep the split step.
"""
const COUPLED = true
"""
Implicit upwind transport of mat.x with the mass that crossed each interface in the converged step (change of the exterior mass m_e), so
that it is the discrete continuity equation applied to the mass fractions: a uniform composition stays uniform and elements are conserved.
"""
function advect_species!()
    mat isa PM.NetGas && PM.HDCHEM[] && PM.hd_pad!(mat)                       # (first step after a restart from a run without HD)
    xme = TF.var(:xmen); old = XMEPREV[]; n = nsh
    F = [xme[a+j-1] - old[a+j-1] for j in 1:n+1]                              # > 0: mass moved outwards through interface j
    M = [xme[a+j-1] - xme[a+j] for j in 1:n]
    lo = zeros(n); up = zeros(n); di = zeros(n)
    for j in 1:n
        di[j] = M[j] + max(F[j+1], 0.0) + max(-F[j], 0.0)
        j > 1 ? (lo[j] = -max(F[j], 0.0)) : (di[j] -= max(F[j], 0.0))         # (no inflow through a boundary: own composition)
        j < n ? (up[j] = -max(-F[j+1], 0.0)) : (di[j] -= max(-F[j+1], 0.0))
    end
    M0 = [old[a+j-1] - old[a+j] for j in 1:n]
    cp_ = zeros(n); den = zeros(n)
    den[1] = di[1]; cp_[1] = up[1] / den[1]
    for j in 2:n; den[j] = di[j] - lo[j] * cp_[j-1]; cp_[j] = up[j] / den[j]; end
    xnew = [M0[j] .* mat.x[a+j-1] for j in 1:n]
    xnew[1] = xnew[1] ./ den[1]
    for j in 2:n; xnew[j] = (xnew[j] .- lo[j] .* xnew[j-1]) ./ den[j]; end
    for j in n-1:-1:1; xnew[j] = xnew[j] .- cp_[j] .* xnew[j+1]; end
    XADV2 && mat isa PM.NetGas && xadv2_correct!(xnew, [copy(mat.x[a+j-1]) for j in 1:n], F, M, n)
    for j in 1:n; mat.x[a+j-1] = xnew[j]; end
    if HOTON                                                                 # the nuclear abundances move with the same fluxes
        Y = SG[].Y; ynew = M0[1] .* Y[:, a]; ys = Vector{Vector{Float64}}(undef, n); ys[1] = ynew ./ den[1]
        for j in 2:n; ys[j] = (M0[j] .* Y[:, a+j-1] .- lo[j] .* ys[j-1]) ./ den[j]; end
        for j in n-1:-1:1; ys[j] = ys[j] .- cp_[j] .* ys[j+1]; end
        for j in 1:n; Y[:, a+j-1] .= ys[j]; end
        SM.ghosts!(SG[])
    end
    copyto!(old, xme)
    nothing
end
function envflags!()
    TX.REGULAR_CENTRE[] = ENVCENTRE; empty!(TX.XMINT[]); TX.GASROW[] = true; TX.INNER_SINK[] = !ENVCENTRE
    PM.CHEMFALLBACK[] = true
    XMEPREV[] = copy(TF.var(:xmen)::Vector{Float64})
end
"ln r of the fixed mesh's interfaces from the first interface (r_in, or ENVCORE's first interface for a regular centre) to rout (ENVEUL/ENVCORE spacing)."
function envlr(rout)
    rin, dl1, dl2, rb = ENVEUL[1] * Rsun, ENVEUL[2], ENVEUL[3], ENVEUL[4] * Rsun
    rcb = ENVCENTRE ? ENVCORE[3] * Rsun : 0.0
    lr = [log(ENVCENTRE ? ENVCORE[1] * Rsun : rin)]; dl = ENVCENTRE ? ENVCORE[2] : dl1
    gout_ = length(ENVEUL) >= 5 ? ENVEUL[5] : 1.02; gin_ = length(ENVEUL) >= 6 ? ENVEUL[6] : 0.95
    while lr[end] < log(rout)
        r_ = exp(lr[end])
        r_ > rb ? (dl = min(dl2, gout_ * dl)) : r_ > rcb && (dl = max(dl1, gin_ * dl))
        push!(lr, lr[end] + dl)
    end
    (length(lr) + (ENVCENTRE ? ENVLIN : 0)) % 2 == 0 && pop!(lr)                    # (an even zone count: an odd number of interfaces)
    lr[1] .+ (lr .- lr[1]) .* (log(rout) - lr[1]) / (lr[end] - lr[1])
end
# ICMLOG=m1 (default unset = unchanged): the IC mesh is logarithmic in enclosed mass instead of envlr's / the rule's (GRIDEQ) mesh:
# N = GEQN (default nsh) zones: the N interfaces above the centre at m_k = m1 (M_cloud/m1)^((k-1)/(N-1)), k = 1..N (m1 in Msun; the centre zone holds m1),
# r(m) the exact inverse of the IC's M(r) interpolation. Only for the start from the IC. The first interface is Lagrangian under GRIDEQ
# (TitanResidFast.GEQ_LAGFIRST, set from the presence of ICMLOG in the environment: restarts of the run follow).
# The start (every run without RESTART/STARFROM; formerly ICEUL=1): the initial condition directly on the fixed log mesh of
# ENVEUL/ENVCORE (lgrid = 2, regular centre), at t = 0. The mesh is `envlr`'s, from ENVCORE's first interface to the outer radius of the IC. Mass-conserving: the interior mass at
# every interface is that of the IC's tabulated M(r) (piecewise-uniform density between its nodes, so the zone masses add up to the
# cloud mass exactly); T and the H₂ fraction are the mass-weighted zone means of the IC's T(M), f_H2(M); u at the interfaces from the
# IC's v(M); x_e as the IC (XE0). Radiation: E_r = a [T_CMB⁴ + (T⁴ − T_CMB⁴)(1 − e^{−τ_R})], F_r flux-limited diffusion.
# All fixed-mesh settings (GASROW, REGULAR_CENTRE, no INNER_SINK, xmen/xmtot interior mass, CHEMFALLBACK, empty
# ALI factors). Only the abn02* initial conditions (they carry the tabulated M(r)).
"""
Zero point of the mass coordinate at the centre (fixed mesh with a regular centre; used by ICEUL, ICXMZ=0 turns it off). TITAN's
mass unknown is the exterior mass x_m,e = xmtot − m (per steradian), with xmtot ≥ the grid mass (PS.build!: grid + xmext), and the
zone masses are its differences. At the IC the central zones of the log mesh hold ~1e-17 of the cloud: zone 39 of the 100-zone mesh
(r ≈ 1.5e3 R☉) has 4.5e20 g sr⁻¹ against an ulp of 5.9e20 of x_m,e ≈ 3.2e36, so its mass row is round-off and Newton diverged at
step 1 whatever Δt (T, E_r corrections −9.4 there). Here xmtot := ε (half the central zone's mass) and x_m,e = ε − m: negative outside
the first interface, but the residual only uses differences of x_m,e, xmtot − x_m,e and relative (sign-preserving) corrections,
so the interior mass is carried at full relative precision at every radius. Recomputed from the interior masses m (summed from the
centre by PS.build!), not by shifting x_m,e.
"""
function xmzero!()
    xmen = TF.var(:xmen); xmeo = TF.var(:xmeo); dn = TF.var(:dn); dv = TF.var(:dvoln); xmt0 = TF.getscalar(:xmtot)
    TF.var(:rn)[a] == 0 || error("xmzero!: needs a regular centre (r = 0 at the first interface)")
    m = zeros(length(xmen)); for k in a+1:b+2; m[k] = m[k-1] + dn[k-1] * dv[k-1]; end      # interior mass, summed from the centre
    ε = 0.5 * m[a+1]
    for k in a:b+2; xmen[k] = ε - m[k]; end
    xmen[a-1] = 2xmen[a] - xmen[a+1]                                                   # (the regular centre's ghost, as TitanCore's bndry)
    for k in a-1:b+2; xmeo[k] = xmen[k]; end
    TF.setscalar!(:xmtot, ε); TF.setscalar!(:xmext, ε - m[b+2])
    pr("XMZERO: mass zero point at the centre, xmtot %.4e (was %.4e) g/sr; x_m,e at the interfaces %d..%d: %s; grid mass %.6f Msun\n", ε, xmt0, a - 1, a + 2,
       join((@sprintf("%.4e", xmen[k]) for k in a-1:a+2), " "), 4π * m[b+1] / Msun)
    nothing
end
function iceulbuild!(; rIset = nothing, icmlog = nothing)     # rIset: the interfaces (centre 0 … outer edge) in place of envlr's mesh (GRIDEQ's starting mesh); icmlog = (m1 [Msun], N): ICMLOG's mass mesh
    (RAD && ENVCENTRE) || error("the start needs RAD=1 and ENVEUL=0,... (regular centre, ENVCORE)")
    radsettings!()
    rout = TF.var(:rn)[b+1]; Mtot = M                                          # the IC's outer radius and cloud mass
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    rg3 = rg .^ 3
    Mof(r) = (r <= 0 && return 0.0; r >= rg[end] && return Mg[end]; j = clamp(searchsortedlast(rg, r), 1, length(rg) - 1);
              Mg[j] + (r^3 - rg3[j]) / (rg3[j+1] - rg3[j]) * (Mg[j+1] - Mg[j]))
    # exact inverse of Mof (r^3 linear in m between the table nodes)
    rofM(m) = (j = clamp(searchsortedlast(Mg, m), 1, length(Mg) - 1); cbrt(rg3[j] + (m - Mg[j]) / (Mg[j+1] - Mg[j]) * (rg3[j+1] - rg3[j])))
    if icmlog !== nothing       # ICMLOG: the N interfaces above the centre at log-uniform enclosed mass from m1 to Mtot (N - 1 zones above the centre zone, which holds m1)
        m1_, N_ = icmlog[1] * Msun, icmlog[2]
        mk_ = [m1_ * (Mtot / m1_)^((k - 1) / (N_ - 1)) for k in 1:N_]
        rk_ = [rofM(m) for m in mk_]; rk_[end] = rout
        rk_[1] > 0 || error("ICMLOG: r(m1) <= 0")
        all(>(0), diff(rk_)) || error("ICMLOG: the interfaces are not monotone")
        for k in 1:N_-1; abs(Mof(rk_[k]) - mk_[k]) <= 1e-10 * mk_[k] || error("ICMLOG: Mof(r(m)) does not reproduce m at interface $(k+1)"); end
        rIset = vcat(0.0, rk_)
    end
    rI_ = rIset === nothing ? (r_ = exp.(envlr(rout)); r_[end] = rout; vcat([r_[1] * i / ENVLIN for i in 0:ENVLIN-1], r_)) : copy(rIset); n = length(rI_) - 1
    mI_ = [Mof(r) for r in rI_]; mI_[end] = Mtot
    all(>(0), diff(mI_)) || error("ICEUL: a zone without mass")
    if icmlog !== nothing
        dl_ = diff(log.(rI_[2:end])); i2_ = searchsortedfirst(mI_, 1e-2Msun) - 1       # (zone j spans interfaces j, j+1 of rI_; i2_ = the zone containing 1e-2 Msun)
        pr("ICMLOG: IC starting mesh, %d zones: first interface (Lagrangian) at m %.3g Msun, r %.4g R_sun, then log in enclosed mass to %.4g Msun (r %.4g R_sun), dln r %.3f .. %.3f, dln m %.4f\n",
           n, mI_[2] / Msun, rI_[2] / Rsun, Mtot / Msun, rI_[end] / Rsun, minimum(dl_), maximum(dl_), log(Mtot / mI_[2]) / (n - 1))
        pr("ICMLOG: 1e-2 Msun at r %.4g R_sun, in zone %d (of %d); interfaces 1,2,3,4 r [R_sun] %s\n", rofM(1e-2Msun) / Rsun, i2_, n,
           join((@sprintf("%.4g", rI_[k] / Rsun) for k in 2:5), " "))
    end
    ρz = [(mI_[j+1] - mI_[j]) / (4π / 3 * (rI_[j+1]^3 - rI_[j]^3)) for j in 1:n]
    # mass-weighted zone means of T(M), f_H2(M) (64 midpoints in m per zone; the IC tables are in log M [M☉], constant below their first node)
    zmean(f, j) = (ms = range(mI_[j], mI_[j+1], length = 65); sum(f(0.5 * (ms[i] + ms[i+1]) / Msun) for i in 1:64) / 64)
    fin_ = ic in ("abn02fin", "abn02red")
    Tz = fin_ ? [zmean(Tfin, j) for j in 1:n] : [zmean(m -> interp(log10(m), lM_tab, T_tab), j) for j in 1:n]
    fz = fin_ ? [zmean(m -> finy(m)[3], j) for j in 1:n] : [zmean(m -> 10^interp(log10(m), lM_tab, lf_tab), j) for j in 1:n]   # (fin: x_H2 per H)
    uI_ = fin_ ? [vfin(mI_[j] / Msun, rI_[j]) for j in 1:n+1] : ic == "abn02z19" ? [m < Mjoin ? 0.0 : 1e5 * interp(log10(m / Msun), lM_tab, v_tab) * min(1.0, (m - Mjoin) / Mjoin) for m in mI_] :
          [1e5 * interp(log10(max(m / Msun, 1e-2)), lM_tab, v_tab) for m in mI_]
    uI_[1] = 0.0; uI_[end] = 0.0                                               # regular centre; fixed outer wall (as the remap)
    xz = fin_ ? [xinit((zmean(m -> finy(m)[1], j), 0.0, fz[j], 0.0)) for j in 1:n] : [xinit((xe0, 0.0, f / 2, 0.0)) for f in fz]; eqz = falses(n)
    # initial radiation field: diffusion inside, the CMB and the free-streaming photospheric light outside
    Tb = parse(Float64, get(ENV, "TCMB", "0")); ar = Z0.arad
    χ = [PM.OPAC[](ρz[j], Tz[j])[1] * ρz[j] for j in 1:n]
    τ = zeros(n); for j in n-1:-1:1; τ[j] = τ[j+1] + 0.5 * (χ[j] * (rI_[j+1] - rI_[j]) + χ[j+1] * (rI_[j+2] - rI_[j+1])); end
    Er = [ar * (Tb^4 + (Tz[j]^4 - Tb^4) * (-expm1(-τ[j]))) for j in 1:n]
    rc = [0.5 * (rI_[j] + rI_[j+1]) for j in 1:n]; Fr = zeros(n + 1)
    for j in 2:n
        gE = (Er[j] - Er[j-1]) / (rc[j] - rc[j-1]); Em = 0.5 * (Er[j] + Er[j-1])
        Fr[j] = -Z0.c * gE / (3 * 0.5 * (χ[j] + χ[j-1]) + abs(gE) / Em) * (-expm1(-τ[j]))
    end
    Fr[n+1] = Fr[n]
    Mnew = sum(4π / 3 * ρz[j] * (rI_[j+1]^3 - rI_[j]^3) for j in 1:n)
    TX.GASROW[] = true
    for k in eachindex(mat.x); j = clamp(k - a + 1, 1, n); mat.x[k] = copy(xz[j]); mat.eqzone[k] = eqz[j]; end
    PS.build!(; rI = rI_, rho = ρz, T = Tz, u = uI_, geom = 2, rad = true, Er, Fr, lum = nothing, lgrid = 2, inner = :wall, outer = :wall, Mc = 0.0,
              selfgrav = true, material = mat, av, xmext = Mnew / 4π)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    get(ENV, "ICXMZ", "1") == "1" && xmzero!()
    for j in 1:nsh
        k = a + j - 1; mat.x[k] = xz[j]; mat.eqzone[k] = eqz[j]; mat.Λext[k] = 0.0; mat.Tref[k] = Tz[j]; mat.slope[k] = 0.0
    end
    for k in 1:a-1; mat.x[k] = copy(mat.x[a]); mat.eqzone[k] = mat.eqzone[a]; end
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    PS.install_julia!(BlockSolver)
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, 0.0); TF.setscalar!(:timeo, 0.0); TF.setscalar!(:jstep, 0)
    TF.setscalar!(:dtime, haskey(ENV, "ENVDT0") ? parse(Float64, ENV["ENVDT0"]) : 1e-3tff); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    envflags!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    POPS0[] = zeros(0, 0); SNAPS === nothing || empty!(SNAPS)
    xme = TF.var(:xmen); ENVMS0[] = xme[a] - xme[a+envsurface()-1]
    pr("ICEUL: initial condition on the fixed mesh: %d zones, r 0, %.3g R_sun .. %.4g pc (dln r %.4f), grid mass %.2f Msun (cloud %.2f), rho_c %.4e T_c %.1f, tau_R(centre) %.3e, %d interfaces falling, dt0 %.3e s, M_shell0 %.3e Msun\n",
       n, rI_[2] / Rsun, rI_[end] / 3.0857e18, log(rI_[3] / rI_[2]), Mnew / Msun, Mtot / Msun, ρz[1], Tz[1], τ[1], count(<(0), uI_), TF.getscalar(:dtime), 4π * ENVMS0[] / Msun)
end
const LNUC = Ref(0.0)
const TSTEP = Ref(0.0)
const TCHEM = Ref(0.0)
const TLINE = Ref(0.0)
const TOTH = Ref(0.0)
const TOUT = Ref(0.0)
const KLAST = Ref(0)
const LASTOK = Ref{Any}(nothing)
# NWFD (diagnostic): for the unknowns j of every zone z in zlo:zhi, the central difference (h = 1e-5) of (a) the residual with the Newton context frozen (what
# JACCHK does: everything the iterate's update! refreshes outside residual! -- EOS/opacity fills, Eddington factors f, g, h, the flux-tie weights -- held),
# and (b) the full nonlinear map: the state is perturbed with update! (EOS, opacity, formal solution refreshed), a fresh context is built and the residual
# evaluated at x = 0; the state is restored by the inverse update. Both against the assembled Jacobian of the iterate. Prints every (row, zone, var) with
# |FD - AD| > 1e-3 x (largest diagonal entry of the row) for rows in zones z-3:z+3, for (a) and (b) separately, and the (b)-(a) difference.
function nwfd!(zlo, zhi, iter)
    a_ = Int(TF.getscalar(:ngrs)); c = TX.context(); ne = c.neqn; nb = c.nb; k0 = c.k0; ws = TX.getws(Float64, c.mgr)
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2)); AD = [copy(e) for e in E]
    rhs = TF.var(:rhs); h = parse(Float64, get(ENV, "NWFDH", "1e-5"))
    pr("NWFD: step dt %.4e theta %.2f iteration %d, neqn %d, zones %d..%d (relative numbering; absolute = +%d), code %s\n", TF.getscalar(:dtime), TF.getscalar(:thet), iter, ne, zlo, zhi, a_ - 1, string(c.code))
    base_ = (; f = copy(TF.var(:fedd)), g = TF.getscalar(:geddr), h = copy(TD.HEDD), hx = copy(TD.HXEDD), fx = copy(TD.FXEDD), gx = TD.GXEDD[], tie = copy(TD.TIEW))
    FZ = Ref(String[])
    function resid_now()
        for n_ in FZ[]
            n_ == "f" && copyto!(TF.var(:fedd), base_.f); n_ == "g" && TF.setscalar!(:geddr, base_.g)
            n_ == "h" && copyto!(TD.HEDD, base_.h); n_ == "hx" && copyto!(TD.HXEDD, base_.hx); n_ == "fx" && copyto!(TD.FXEDD, base_.fx)
            n_ == "gx" && (TD.GXEDD[] = base_.gx); n_ == "tie" && copyto!(TD.TIEW, base_.tie)
        end
        cc = TX.context(); Rr = zeros(ne, nb); TX.residual!(Rr, ws, cc, zeros(ne * nb)); Rr
    end
    function bump!(za, j, δ)
        fill!(rhs, 0.0); rhs[j, za] = δ; TC.update!()
    end
    revert(δ, j) = (j <= 3 || j == 5 || j == 6) ? -δ / (1 + δ) : -δ
    Rf = zeros(ne, nb)
    R0 = resid_now()
    sets_ = [s_ == "none" ? String[] : String.(split(s_, ",")) for s_ in split(get(ENV, "NWFDSETS", "none"), ";")]
    for (iset, fz_) in enumerate(sets_)
    FZ[] = fz_; ACC = Dict{Tuple{Int,Int},Tuple{Float64,Int,Float64}}()
    pr("NWFD set %d: closure quantities held at the pre-perturbation values in the full FD: %s\n", iset, isempty(fz_) ? "none (full refresh)" : join(fz_, ","))
    for z in zlo:zhi
        za = a_ + z - 1; kbz = za - k0 + 1
        for j in 1:ne
            xv = zeros(ne * nb); ix = (kbz - 1) * ne + j
            Rp = zeros(ne, nb); Rm = zeros(ne, nb)
            xv[ix] = h; TX.residual!(Rp, ws, c, xv); xv[ix] = -h; TX.residual!(Rm, ws, c, xv)
            FDa = (Rp .- Rm) ./ 2h
            bump!(za, j, h); Rp2 = resid_now(); bump!(za, j, revert(h, j))
            bump!(za, j, -h); Rm2 = resid_now(); bump!(za, j, revert(-h, j))
            FDb = (Rp2 .- Rm2) ./ 2h
            for kb in max(1, kbz - 3):min(nb, kbz + 3), i in 1:ne
                o = kbz - kb; ad = abs(o) <= 2 ? AD[o + 3][i, j, k0 + kb - 1] : 0.0
                scale = maximum(abs, (AD[3][i, jj, k0 + kb - 1] for jj in 1:ne)) + 1e-300
                ea = abs(FDa[i, kb] - ad) > 1e-3 * scale; eb = abs(FDb[i, kb] - ad) > 1e-3 * scale; ec = abs(FDb[i, kb] - FDa[i, kb]) > 1e-3 * scale
                let v_ = abs(FDb[i, kb] - ad) / scale; o_ = get(ACC, (i, j), (0.0, 0, 0.0)); ACC[(i, j)] = (max(o_[1], v_), o_[2] + (v_ > 1e-3 ? 1 : 0), max(o_[3], abs(FDa[i, kb] - ad) / scale)); end
                (iset == 1 && (ea || eb || ec)) && pr("   row %d zone %d | col var %d zone %d (offset %d): FDfrozen %.5e FDfull %.5e AD %.5e diagscale %.3e | %s%s%s\n", i, kb + k0 - a_, j, z, o, FDa[i, kb], FDb[i, kb], ad, scale,
                                       ea ? "frozen-vs-AD " : "", eb ? "full-vs-AD " : "", ec ? "full-vs-frozen" : "")
            end
        end
    end
    for ((i_, j_), (mx_, n_, mf_)) in sort(collect(ACC)); (mx_ > 1e-3 || mf_ > 1e-3) && pr("   SET%d row %d col var %d: max |FDfull-AD|/diag %.4g (%d entries > 1e-3), max |FDfrozen-AD|/diag %.4g\n", iset, i_, j_, mx_, n_, mf_); end
    end
    pr("NWFD: done; state drift check, max |R(after) - R(before)| / (max |R|) = %.3e\n", maximum(abs, resid_now() .- R0) / (maximum(abs, R0) + 1e-300))
    nothing
end
const LINECAPV = Ref(NaN)
const JACSAVE = Ref{Any}(nothing)
const RUNEND = Ref{Any}(nothing)
"""
TABAUDIT=1 (diagnostic): for every table zone below the hot branch (and, with TABAUDIT=2, also every network zone) compare at the zone's (ρ, T)
the equilibrium table (PopIIIEqEOS composition, u, p; StellarOpacity κ_R, κ_P) with the network's own LTE state (ChemNet.steady with J = B(T),
which by detailed balance is also the thick-limit three-level state): x_e, x_H2, x_H-, He stages, u, p, c_v and Γ1 (derivatives at
equilibrium composition, finite differences of the steady state), and the continuum means of PopIIIMaterial.nlte_means for the network
composition and for the table composition. Prints one line per zone, then exits.
"""
function tabaudit()
    B(T) = ν -> 2 * ChemNet.hP * ν^3 / ChemNet.c^2 / expm1(ChemNet.hP * ν / (ChemNet.kB * T))
    net = mat.net; ph = filter(p -> true, mat.photos); dn_ = TF.var(:dn); tn_ = TF.var(:tn); rn_ = TF.var(:rn); h = PM.HOT[]
    lte(ρ, T, x0) = ChemNet.steady(x0, ρ, T, net; env = (Trad = T, J = B(T), static = true), photos = ph)
    E21 = PM.E21; kB = ChemNet.kB
    function netstate(ρ, T, x0)
        x = lte(ρ, T, x0); u, p = ChemNet.eos(ρ, T, x, net); x, u, p
    end
    function means(ρ, T, d)
        nH = ρ / ChemNet.μH(); n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
        nHe = nH * max(get(d, :He, PC.yHe), 0.0); nHm = nH * max(get(d, :Hm, 0.0), 0.0)
        PM.nlte_means(ρ, T; n1, n2 = 4n1 * exp(-E21 / (kB * T)), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, hbfP = true)
    end
    @printf("TABAUDIT: zone r[Rsun] T[K] nH | x_e tab/net | x_H2 tab/net | x_H- tab/net | x_He+ x_He++ tab/net | u tab/net (rel diff) | p rel diff | c_v rel diff | Γ1 tab/net | κR_tab κP_tab | κ(E,P,R) nlte net comp | κ(E,P,R) nlte tab comp\n")
    if get(ENV, "TABAUDIT", "1") == "4"          # opacity audit of every zone below the hot branch (table and network zones)
        X = 1 / (1 + 4 * PC.yHe); dlr(ρ, T) = log10(ρ) - 3 * (log10(T) - 6)
        @printf("TABAUD4: zone r[Rsun] T nH eq | logT logR flag(1 OPAL 2 ext 3 blend 4 own) | OPAC κR κP | run κR κE κP | nlte(net comp) κE κP κR | nlte(LTE comp) κE κP κR | x_e x_H2 x_H- net\n")
        for k in a:b
            T = tn_[k]; (h !== nothing && T > h.T1) && continue
            ρ = dn_[k]; nH = ρ / ChemNet.μH(); lT = log10(T); lR = dlr(ρ, T)
            try
                fl = StellarOpacity.coverage(lT, lR, X); κR, κP = PM.OPAC[](ρ, T)
                TX.ZONE[] = k; rr = TX.material_opac(mat, ρ, T)
                dn2 = ChemNet.full(mat.x[k], net); dt_ = ChemNet.full(PM.x_from_table(mat, ρ, T), net)
                mn = means(ρ, T, dn2); mt = means(ρ, T, dt_); g(d, s_) = Float64(get(d, s_, 0.0))
                @printf("TABAUD4 %3d %9.4g %8.0f %.2e %d | %.3f %.2f %d | %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3e %.2e %.2e\n",
                        k - a + 1, rn_[k] / 6.957e10, T, nH, mat.eqzone[k], lT, lR, fl, κR, κP, rr..., mn..., mt..., g(dn2, :e), g(dn2, :H2), g(dn2, :Hm))
            catch e
                println("TABAUD4 $(k - a + 1) failed: ", first(sprint(showerror, e), 120))
            end
        end
        return nothing
    end
    # TABAUDIT=3: a (n_H, T) grid instead of the zones (where the ideal network departs from the table's non-ideal SCVH)
    pts = get(ENV, "TABAUDIT", "1") == "3" ? [(0, nH_ * ChemNet.μH(), T_) for nH_ in 10.0 .^ (parse(Float64, get(ENV, "TABAUDNH", "17")):0.5:24) for T_ in (6e3, 8e3, 1e4, 1.5e4, 2e4, 3e4, 5e4, 1e5, 2e5, 2.9e5)] :
          [(k, dn_[k], tn_[k]) for k in a:b if (mat.eqzone[k] || get(ENV, "TABAUDIT", "1") == "2") && !(h !== nothing && tn_[k] > h.T1)]
    for (k, ρ, T) in pts
        nH = ρ / ChemNet.μH()
        try
            xt = PM.x_from_table(mat, ρ, T); dt_ = ChemNet.full(xt, net)
            ut, pt = PopIIIEqEOS.eos(ρ, T)
            fT(t) = PopIIIEqEOS.eos(ρ, t); fR(r) = PopIIIEqEOS.eos(r, T)
            uT = ForwardDiff.derivative(t -> fT(t)[1], T); pT = ForwardDiff.derivative(t -> fT(t)[2], T)
            uR = ForwardDiff.derivative(r -> fR(r)[1], ρ); pR = ForwardDiff.derivative(r -> fR(r)[2], ρ)
            Γt = ρ / pt * (pR + pT * (pt / ρ^2 - uR) / uT)
            xn, un, pn = netstate(ρ, T, xt); dn2 = ChemNet.full(xn, net)
            δ = 1e-3
            _, u1, p1 = netstate(ρ, T * (1 + δ), xn); _, u0, p0 = netstate(ρ, T * (1 - δ), xn)
            _, u3, p3 = netstate(ρ * (1 + δ), T, xn); _, u2, p2 = netstate(ρ * (1 - δ), T, xn)
            nuT = (u1 - u0) / (2δ * T); npT = (p1 - p0) / (2δ * T); nuR = (u3 - u2) / (2δ * ρ); npR = (p3 - p2) / (2δ * ρ)
            Γn = ρ / pn * (npR + npT * (pn / ρ^2 - nuR) / nuT)
            κR, κP = PM.OPAC[](ρ, T); mn = means(ρ, T, dn2); mt = means(ρ, T, dt_)
            g(d, s) = Float64(get(d, s, 0.0))
            @printf("TABAUDIT %3d %8.3g %8.0f %.2e | %.3e %.3e | %.2e %.2e | %.2e %.2e | %.2e %.2e  %.2e %.2e | %.4e %.4e (%+.2e) | %+.2e | %+.2e | %.4f %.4f | %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e\n",
                    k == 0 ? 0 : k - a + 1, k == 0 ? ρ : rn_[k] / 6.957e10, T, nH, g(dt_, :e), g(dn2, :e), g(dt_, :H2), g(dn2, :H2), g(dt_, :Hm), g(dn2, :Hm),
                    g(dt_, :Hep), g(dt_, :Hepp), g(dn2, :Hep), g(dn2, :Hepp), ut, un, (un - ut) / abs(ut), (pn - pt) / pt, (nuT - uT) / uT,
                    Γt, Γn, κR, κP, mn..., mt...)
        catch e
            println("TABAUDIT $(k == 0 ? 0 : k - a + 1) T $T nH $nH: failed: ", first(sprint(showerror, e), 120))
        end
    end
end
const PROFT = Dict{Symbol,Float64}()
const PROFN = Dict{Symbol,Int}()
const PROFSTEP = Ref(-1)
const PROFB = Dict{Symbol,Float64}()
# (PERF) bytes allocated inside each routine (all threads, Base.gc_bytes)
"""
MGROUPS (docs/POPIII_MULTIGROUP.md): G radiation groups as Newton unknowns 10 … 9 + 2G (TitanResidFast MGRad.jl), the group opacities and
emissivities from the material at the Newton iterate (PopIIIMaterial material_opacg!, MGOpac.jl), the star-light rates from the group
fields (starlight_setup), the pivoted banded solver. The group fields come from the checkpoint (MGAUX) or, the first time, from the grey field (split by the local emission balance of each group, docs §25).
"""
function mg_install!()
    COUPLED || error("MGROUPS needs COUPLED=1")
    Int(TF.getscalar(:lrad)) == 1 || error("MGROUPS needs the radiation field (RAD=1)")
    (TX.COMPTON[] || TX.INNER_SINK[] || !isnan(TX.PEXT[])) && error("MGROUPS: COMPTON, INNER_SINK, PEXT are not supported")
    edges = PM.mg_edges(ENV["MGROUPS"]); G = length(edges) + 1
    TF.PARAMS[:meqn] >= 9 + 2G || error("MGROUPS: meqn = $(TF.PARAMS[:meqn]) < $(9 + 2G)")
    tcmb = parse(Float64, get(ENV, "TCMB", "0"))
    TX.mg_setup!(edges, PM.EVERG / ChemNet.kB, tcmb, 9, length(TF.var(:dn)))
    TF.setscalar!(:neqn, 9 + 2G)
    PM.G2O.FASTON || error("MGROUPS needs G2OPAC=2 (the group means come from G2Opacity.group_fast)")
    PM.G2B === nothing || error("MGROUPS does not support G2BLEND")
    PM.SLIMP || error("MGROUPS needs SLIMP=2 (the chemistry's radiation is the group fields at the Newton iterate)")
    PM.mg_init!(mat.net.free); r_ = PM.MGROLE[]; TX.MGQG[] = r_.Lya; TX.MGSLG[] = r_.Hm
    mat.net.free[ChemNet.ckinfo(mat.net).slow[2]] === :H2 || error("MGROUPS: the second coupled species must be H2 (the LW closure reads it)")
    TX.MGXH[] = PM.Z0Physics.χH / PM.EVERG; TX.MGGH[] = r_.H; TX.MGLG[] = 1
    TX.mg_valloc!(length(TF.var(:dn)))
    TX.MGSLQ[] = PM.sl_qt_mg; TX.MGH2FN[] = PM.mg_h2iter!
    TX.MGZONEQ[] = (k, ρ, T, a, b, qv, t, g, dp, xn) -> PM.mg_hasrec(PM.CPLTERMS[][k]) ? PM.cpl_zone_qg(mat, k, ρ, T, a, b, qv, t, g, dp, xn) :
                                                   PM.cpl_zone_q(mat, k, ρ, T, a, b, one(t), t, zero(t), g, dp, xn)      # (a zone without a star-light record: the grey path)
    TD.ROUTINES[:matslv] = JuliaSolver.matslv!        # banded LU with partial pivoting (the block solver pivots only inside a block: diverges for neqn = 9 + 2G, t_mg23)
    let f0 = TC.EDDFAC[]; TC.EDDFAC[] = () -> mg_eddfac!(f0) end      # Eddington factors per group
    h0 = TC.CPLH[]
    TC.CPLH[] = (; newstep = () -> (h0.newstep(); TX.mg_newstep!()), reset = () -> (h0.reset(); TX.mg_reset!()),
                 update = (rhs, a_, b_) -> (h0.update(rhs, a_, b_); TX.mg_update!(rhs, a_, b_)),
                 dmax = (rhs, a_, b_) -> max(h0.dmax(rhs, a_, b_), TX.mg_dmax!(rhs, a_, b_)), predicted = TX.mg_predicted!)
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    ax = MGAUX[]
    if ax !== nothing && ax.edges == TX.MGEDGE[]
        for g in 1:G; copyto!(TX.MGEN[][g], ax.E[g]); copyto!(TX.MGFN[][g], ax.F[g]); copyto!(TX.MGEO[][g], ax.E[g]); copyto!(TX.MGFO[][g], ax.F[g]); end
        TX.mg_newstep!(); src = "the checkpoint"
    else
        PM.nloc_reset!(length(mat.x)); PM.mg_fill!(mat)        # (the group records the balance split reads)
        TX.mg_init_from_grey!(); src = "the grey field (split by the local emission balance of each group, then solved: docs §26)"
        let f1 = TD.STEPFN[], first = Ref(true)       # the radiation-row solve runs once, inside the first step call: the chemistry's per-step preparation (cpl_prepare!) must have run
            TD.STEPFN[] = function ()
                first[] && (first[] = false; PM.cpl_prepare!(mat); mg_radsolve!())
                f1()
            end
        end
    end
    PM.nloc_reset!(length(mat.x))                    # (a restart evaluated the opacities once before the groups existed: no group records in that cache)
    PM.mg_fill!(mat); copyto!(TX.MGVO[], TX.MGVN[])
    pr("MGROUPS: %d groups, edges %s eV, unknowns 10..%d (neqn %d, meqn %d), q from group %d, chemistry energy from group %d; fields from %s\n",
       G, join((@sprintf("%.4g", e) for e in edges), ","), 9 + 2G, 9 + 2G, TF.PARAMS[:meqn], TX.MGSLG[], TX.MGQG[], src)
    nothing
end
"""
(MGROUPS, restart from a checkpoint without group fields; docs §26) Solve the radiation rows alone at the held gas state: the unknowns E_g, F_g of all groups and zones, the gas
state (ρ, T, u, r, m, x_H⁺, x_H₂, E_r, F_r) held, to the true-residual tolerance of the Newton iteration. The residual is the group rows of the step (`mg_rows!`) with their own time
derivatives removed (`TX.MGSTEADY`: Δt → ∞ in the E_g time term, the flux inertia and the tie's grid term), at θ = 1; the band Jacobian is that of the full system with the rows and columns of the
non-radiation unknowns removed (`TX.mg_restrict_band!`), the Krylov iteration (exact J·v by AD) runs on the group unknowns only (`NewtonKrylov.RADONLY`). The split of `mg_init_from_grey!` is the Newton
initial guess. Afterwards the old level is the new one (groups and the grey E_r, F_r), and the report line gives ΣE_g − E_r against the gas and the total energy.
"""
function mg_radsolve!()
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); G = TX.MGN[]
    ern = TF.var(:ern); frn = TF.var(:frn); dvn = TF.var(:dvoln); dn = TF.var(:dn); un = TF.var(:un); egn = TF.var(:egn); cc = TF.getscalar(:cc)
    Er0 = [ern[k] for k in 1:length(ern)]; Fr0 = copy(frn)
    split_ = [copy(TX.MGEN[][g]) for g in 1:G]; fsplit = [copy(TX.MGFN[][g]) for g in 1:G]
    thet0 = TF.getscalar(:thet); dt0 = TF.getscalar(:dtime); nit0 = TF.getscalar(:niter)
    pr("MGRADSOLVE: start; timeo %.6e timen %.6e dtime %.4e thet %.3f conv %.3g niter %d | max|xmen-xmeo|/xmen %.3e\n", TF.getscalar(:timeo), TF.getscalar(:timen), dt0, thet0, TF.getscalar(:conv), Int(nit0),
       maximum(abs((TF.var(:xmen)[k] - TF.var(:xmeo)[k]) / max(abs(TF.var(:xmen)[k]), 1e-300)) for k in a_:b_))
    TX.mg_predicted!()                                 # scales at the state the solve starts from (records after cpl_prepare!)
    g0 = TD.ROUTINES[:matgen]; ok = false; t0 = time()
    TX.MGSTEADY[] = true; NewtonKrylov.RADONLY[] = true
    TD.ROUTINES[:matgen] = () -> (g0(); TX.mg_restrict_band!())
    chk0 = NewtonKrylov.OPT.check; NewtonKrylov.OPT.check = 0           # (NKCHECK belongs to the first step's iteration, not to this solve)
    TF.setscalar!(:thet, 1.0); TF.setscalar!(:niter, typeof(nit0)(parse(Int, get(ENV, "MGRSNITER", "1000"))))
    try
        ok = Trbdf2.newton0!()
    finally
        TX.MGSTEADY[] = false; NewtonKrylov.RADONLY[] = false; TD.ROUTINES[:matgen] = g0; NewtonKrylov.OPT.check = chk0
        TF.setscalar!(:thet, thet0); TF.setscalar!(:dtime, dt0); TF.setscalar!(:niter, nit0)
    end
    pr("MGRADSOLVE: newton %s in %d iterations, true residual %.3e (conv %.3g), last dmax %.3e, %.1f s\n", ok ? "converged" : "FAILED", Int(TF.getscalar(:iter)), Trbdf2.TRLAST[], TF.getscalar(:conv), TF.getscalar(:dmax), time() - t0)
    ok || error("MGRADSOLVE: the steady radiation solve did not converge")
    for g in 1:G; copyto!(TX.MGEO[][g], TX.MGEN[][g]); copyto!(TX.MGFO[][g], TX.MGFN[][g]); end
    copyto!(TF.var(:ero), ern); copyto!(TF.var(:fro), frn)
    size(TX.MGVO[]) == size(TX.MGVN[]) && copyto!(TX.MGVO[], TX.MGVN[])
    TX.mg_newstep!()
    # report
    ΔEr = sum((ern[k] - Er0[k]) * dvn[k] for k in a_:b_); Er = sum(Er0[k] * dvn[k] for k in a_:b_)
    Eg = sum(dn[k] * dvn[k] * egn[k] for k in a_:b_); Ek = sum(dn[k] * dvn[k] * 0.25 * (un[k]^2 + un[k+1]^2) for k in a_:b_)
    pr("MGRADSOLVE: sum_k dV (sum_g E_g - E_r) = %.6e erg of E_r(grey) %.6e (%.3e); gas internal energy %.6e (%.3e), gas internal + kinetic + radiation %.6e (%.3e)\n", ΔEr, Er, ΔEr / Er,
       Eg, ΔEr / Eg, Eg + Ek + Er, ΔEr / (Eg + Ek + Er))
    let (mz, mg_, mk) = (0.0, 0, 0), (rz, rg, rk) = (0.0, 0, 0), (sz, sg, sk) = (0.0, 0, 0)
        for k in a_:b_
            for g in 1:G
                d_ = abs(TX.MGEN[][g][k] - split_[g][k]) / max(abs(split_[g][k]), 1e-300); d_ > mz && (mz = d_; mg_ = g; mk = k - a_ + 1)
            end
        end
        for k in a_+1:b_, g in 1:G
            q = abs(TX.MGFN[][g][k]) / (cc * max(0.5 * (TX.MGEN[][g][k-1] + TX.MGEN[][g][k]), 1e-300)); q > rz && (rz = q; rg = g; rk = k - a_ + 1)
            q = abs(TX.MGFN[][g][k]) / (cc * max(TX.MGEN[][g][k-1], TX.MGEN[][g][k], 1e-300)); q > sz && (sz = q; sg = g; sk = k - a_ + 1)
        end
        pr("MGRADSOLVE: largest |F_g|/(c E_g) (E_g the mean of the interface's two zones) %.4g (group %d, interface %d); with the larger of the two zones %.4g (group %d, interface %d); largest relative change of E_g against the split %.3e (group %d, zone %d)\n",
           rz, rg, rk, sz, sg, sk, mz, mg_, mk)
        for g in 1:G
            pr("MGRADSOLVE: group %d: sum E dV %.6e (split %.6e) erg; max|F|/cE %.3e\n", g, sum(TX.MGEN[][g][k] * dvn[k] for k in a_:b_), sum(split_[g][k] * dvn[k] for k in a_:b_),
               maximum(abs(TX.MGFN[][g][k]) / (cc * max(0.5 * (TX.MGEN[][g][k-1] + TX.MGEN[][g][k]), 1e-300)) for k in a_+1:b_))
        end
    end
    nothing
end
"""
(MGROUPS) Eddington factors per group: after the grey formal solution `f0`, the same formal solution (TitanCore.eddfac!,
tangent rays, linear-in-τ source, VEF_BG background, EDDSCAT) once per group on the group's arrays, swapped into the state arrays it reads:
χ → χ_g, absorption κ_P → κ_E,g, B → emission/κ_E,g of the zone's group record (TX.mg_level on MGVN: the same values as the rows;
the source is (emission + (χ_g − κ_E,g) J_g)/χ_g: absorption κ_E,g, the rest of χ_g scattering), E_r → E_g, F_r → F_g, I_cmb → I_cmb b_g(T_cmb).
Out: f_g, h_g, the tie weights of the group's own depth, g_g. Then the grey state is restored and fedd := Σ f_g E_g / Σ E_g.
"""
function mg_eddfac!(f0)
    f0()
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); G = TX.MGN[]
    ins = (:chifn, :xkpn, :plfn, :ern, :frn, :dtau); sv = Dict(s_ => copy(TF.var(s_)) for s_ in (ins..., :fedd))
    hsv = copy(TD.HEDD); wsv = copy(TD.TIEW); fxsv = copy(TD.FXEDD); gxsv = TD.GXEDD[]
    gr0, gl0, xi0 = TF.getscalar(:geddr), TF.getscalar(:geddl), TF.getscalar(:ximr)
    chi, kp, pl, er, fr = (TF.var(s_) for s_ in (:chifn, :xkpn, :plfn, :ern, :frn)); xke = TF.var(:xken); tn = TF.var(:tn)
    dn, dvol, rmu = TF.var(:dn), TF.var(:dvoln), TF.var(:rmun); V = TX.MGVN[]
    c4pi = TF.getscalar(:cc) / (4π)
    κh, κhg, _ = PM.mg_h2now(); gH2 = TX.MGH2G[]     # H₂ Lyman–Werner absorption and its fluorescence (stage 2c)
    try
        for g in 1:G
            j = TX.mgv(g, 1)
            for k in a_-1:b_+1
                m = TX.mg_zcol(dn, dvol, rmu, k)
                kE, chi[k], kP, a1 = TX.mg_level(V, j, k, sv[:chifn][k], xke[k], sv[:xkpn][k], m)
                kE = max(kE, 1e-300); em = kP * sv[:plfn][k] * TX.mg_bfrac(g, tn[k])
                (g >= TX.MGGH[] && a1 > 0) && (em += a1 * sv[:plfn][k] * TX.mg_xfrac(g, tn[k], TX.MGXH[]))
                # Raman pairs (TX.MGRAMP): absorption and source of the group's share of the transfer
                for (p_, pr_) in enumerate(TX.MGRAMP[])
                    kc = V[TX.mgvp(p_, 1), k]; kw = V[TX.mgvp(p_, 2), k]; (kc > 0 || kw > 0) || continue
                    pr_[1] == g && (kE += kw; em += kc * c4pi * TX.MGEN[][pr_[2]][k])
                    pr_[2] == g && (kE += kc; em += kw * c4pi * TX.MGEN[][pr_[1]][k])
                end
                if gH2 > 0 && κh[k] > 0
                    g == gH2 && (kE += κh[k]; chi[k] += κh[k])
                    g == TX.MGH2T[] && (em += (κh[k] - κhg[k]) * c4pi * TX.MGEN[][gH2][k])
                end
                kp[k] = kE; pl[k] = em / kE
            end
            copyto!(er, TX.MGEN[][g]); copyto!(fr, TX.MGFN[][g]); TF.setscalar!(:ximr, xi0 * TX.MGBCMB_[g])
            TC.HPON[] = true; TD.SWGRP[] = g
            try f0() finally TC.HPON[] = false; TD.SWGRP[] = 0 end
            hp = TX.MGHPLUS[][g]; @inbounds for k in eachindex(hp); hp[k] = k <= length(TC.HPLUS) ? TC.HPLUS[k] : 0.0; end
            copyto!(TX.MGFEDD[][g], TF.var(:fedd))
            length(TD.HEDD) >= b_ + 1 && (h = TX.MGHEDD[][g]; @inbounds for k in eachindex(h); h[k] = k <= length(TD.HEDD) ? TD.HEDD[k] : 0.0; end)
            length(TD.TIEW) >= b_ + 1 && (w = TX.MGTIEW[][g]; @inbounds for k in eachindex(w); w[k] = k <= length(TD.TIEW) ? TD.TIEW[k] : 0.0; end)
            TX.MGGEDDR[][g] = TF.getscalar(:geddr)
        end
    finally
        for s_ in (ins..., :fedd); copyto!(TF.var(s_), sv[s_]); end
        copyto!(TD.HEDD, hsv); copyto!(TD.TIEW, wsv); copyto!(TD.FXEDD, fxsv); TD.GXEDD[] = gxsv
        TF.setscalar!(:geddr, gr0); TF.setscalar!(:geddl, gl0); TF.setscalar!(:ximr, xi0)
    end
    fe = TF.var(:fedd)
    for k in a_-1:b_+2
        sE = 0.0; sF = 0.0
        for g in 1:G; e_ = TX.MGEN[][g][k]; sE += e_; sF += TX.MGFEDD[][g][k] * e_; end
        sE > 0 && (fe[k] = sF / sE)
    end
    TX.MGVEFOK[] = true
    nothing
end
"(MGROUPS) a record per snapshot in POPIII_OUT/mg_rec.jls: groups, grey field, state, ratios, star-light rates."
function mg_record(k)
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); rg = a_-1:b_+2; sd = PM.SLDIAG[]; V = TX.MGVN[]; G = TX.MGN[]
    rat(a, b) = [b[i] > 0 ? a[i] / b[i] : 0.0 for i in eachindex(a)]
    rec = (; step = k, t = TF.getscalar(:timen), dt = TF.getscalar(:timen) - TF.getscalar(:timeo), r = TF.var(:rn)[rg], ρ = TF.var(:dn)[rg], T = TF.var(:tn)[rg],
           E = TF.var(:ern)[rg], F = TF.var(:frn)[rg], chi = TF.var(:chifn)[rg], kE = TF.var(:xken)[rg], kP = TF.var(:xkpn)[rg], fedd = TF.var(:fedd)[rg],
           Eg = [TX.MGEN[][g][rg] for g in 1:TX.MGN[]], Fg = [TX.MGFN[][g][rg] for g in 1:TX.MGN[]],
           rE = [V[TX.mgv(g, 1), rg] for g in 1:G], rP = [V[TX.mgv(g, 2), rg] for g in 1:G], rR = [V[TX.mgv(g, 3), rg] for g in 1:G],
           rEc = [V[TX.mgv(g, 4), rg] for g in 1:G], k1 = [V[TX.mgv(g, 5), rg] for g in 1:G], a1 = [V[TX.mgv(g, 6), rg] for g in 1:G],
           edges = copy(TX.MGEDGE[]), kram = [V[TX.mgv(g, 7), rg] for g in 1:G], rampairs = copy(TX.MGRAMP[]),
           ramC = [rat(V[TX.mgvp(p_, 1), rg], V[TX.mgv(TX.MGRAMP[][p_][1], 7), rg]) for p_ in eachindex(TX.MGRAMP[])],
           ramW = [rat(V[TX.mgvp(p_, 2), rg], V[TX.mgv(TX.MGRAMP[][p_][1], 7), rg]) for p_ in eachindex(TX.MGRAMP[])],
           tw = V[end, rg], lcl = TX.MGLCO[][rg],
           sumE = [sum(TX.mg_bfrac(g, V[end, z] > 0 ? V[end, z] : TF.var(:tn)[z]) * V[TX.mgv(g, 1), z] for g in 1:G) for z in rg],
           u = TF.var(:un)[rg],
           sumP = [sum(TX.mg_bfrac(g, TF.var(:tn)[z]) * V[TX.mgv(g, 2), z] for g in 1:G) for z in rg],
           tl3 = sd === nothing ? nothing : [sd.tl3[z] === nothing ? nothing : (; β = sd.tl3[z].β, RL = sd.tl3[z].RL, ne = sd.tl3[z].ne, P12 = sd.tl3[z].P12) for z in rg],
           fg = TX.MGVEFOK[] ? [TX.MGFEDD[][g][rg] for g in 1:TX.MGN[]] : nothing, gr = copy(TX.MGGEDDR[]),
           hplus = TX.MGVEFOK[] ? [TX.MGHPLUS[][g][rg] for g in 1:TX.MGN[]] : nothing,
           n2 = V[end-1, rg],
           sl = sd === nothing ? nothing : (; kph = sd.kph, Te = sd.Te, L = sd.L, kHm = sd.kHm[rg], kH = sd.kH[rg], kLW = sd.kLW[rg]),
           xHp = [ChemNet.full(mat.x[z], mat.net)[:Hp] for z in a_:b_], xH2 = [ChemNet.full(mat.x[z], mat.net)[:H2] for z in a_:b_],
           wrec = [z <= length(PM.MGWREC[]) ? PM.MGWREC[][z] : nothing for z in rg],
           h2 = (h_ = PM.mg_h2now(); (; κ = h_[1][rg], κg = h_[2][rg], s = h_[3][rg])))
    open(io -> serialize(io, rec), joinpath(get(ENV, "POPIII_OUT", "."), "mg_rec.jls"), "a")
    nothing
end
_isrepoframe(f) = (s = string(f.file); endswith(s, ".jl") && isfile(joinpath(@__DIR__, basename(s))))
function alloc_report(io, rate)
    r = Profile.Allocs.fetch(); byf = Dict{String,Float64}(); byc = Dict{String,Float64}(); tot = 0.0
    for a in r.allocs
        fr = [f for f in a.stacktrace if _isrepoframe(f)]
        key1 = isempty(fr) ? "(no repository frame)" : string(basename(string(fr[1].file)), ":", fr[1].line, " ", fr[1].func)
        keyc = isempty(fr) ? key1 : join((string(basename(string(f.file)), ":", f.line, " ", f.func) for f in fr[1:min(3, end)]), " < ")
        byf[key1] = get(byf, key1, 0.0) + a.size / rate; byc[keyc] = get(byc, keyc, 0.0) + a.size / rate; tot += a.size / rate
    end
    @printf(io, "allocation profile: %.1f GB estimated (%d samples, rate %g)\n\nby innermost repository frame:\n", tot / 1e9, length(r.allocs), rate)
    for (k, v) in sort(collect(byf); by = last, rev = true)[1:min(60, end)]; @printf(io, "%10.3f GB  %5.1f %%  %s\n", v / 1e9, 100v / tot, k); end
    println(io, "\nby the three innermost repository frames:")
    for (k, v) in sort(collect(byc); by = last, rev = true)[1:min(80, end)]; @printf(io, "%10.3f GB  %5.1f %%  %s\n", v / 1e9, 100v / tot, k); end
end
function mainthread_report(io)
    data = Profile.fetch(include_meta = true); dict = Profile.getdict(data)
    nt = Threads.nthreads(:default)
    blocks = Tuple{Int,Int,UInt64,Vector{UInt}}[]; st = 1
    for i in 1:length(data)
        if Profile.is_block_end(data, i)
            tid = Int(data[i - Profile.META_OFFSET_THREADID])
            slp = Int(data[i - Profile.META_OFFSET_SLEEPSTATE])
            clk = UInt64(data[i - Profile.META_OFFSET_CPUCYCLECLOCK])
            ips = data[st:(i - Profile.nmeta - 2)]
            push!(blocks, (tid, slp, clk, ips))
            st = i + 1
        end
    end
    # Group blocks into timer-interrupt ticks (consecutive blocks with non-repeating tid or within 2 ms of cycle clock)
    ticks = Vector{Tuple{Int,Int,Vector{UInt}}}[]
    cur = Tuple{Int,Int,Vector{UInt}}[]; seen = Set{Int}(); last_clk = UInt64(0)
    for (tid, slp, clk, ips) in blocks
        tid <= nt || continue
        if !isempty(cur) && (tid in seen || (clk > last_clk && clk - last_clk > 5_000_000))
            push!(ticks, cur); cur = Tuple{Int,Int,Vector{UInt}}[]; empty!(seen)
        end
        push!(cur, (tid, slp, ips)); push!(seen, tid); last_clk = clk
    end
    isempty(cur) || push!(ticks, cur)
    n_ticks = length(ticks); n_main_only = 0; n_parallel = 0; n_worker_busy = 0
    by_inner = Dict{String,Int}(); by_chain = Dict{String,Int}(); by_func = Dict{String,Int}()
    for tk in ticks
        m1 = nothing; workers_busy = 0
        for (tid, slp, ips) in tk
            rep = StackTraces.StackFrame[]
            for ip in ips
                for f in get(dict, ip, StackTraces.StackFrame[])
                    _isrepoframe(f) && push!(rep, f)
                end
            end
            busy = slp == 1 && !isempty(rep)
            if tid == 1
                busy && (m1 = rep)
            elseif busy
                workers_busy += 1
            end
        end
        workers_busy > 0 && (n_worker_busy += workers_busy)
        if m1 !== nothing && workers_busy == 0
            n_main_only += 1
            k1 = string(basename(string(m1[1].file)), ":", m1[1].line, " ", m1[1].func)
            kc = join((string(basename(string(f.file)), ":", f.line, " ", f.func) for f in m1[1:min(4, end)]), " < ")
            by_inner[k1] = get(by_inner, k1, 0) + 1
            by_chain[kc] = get(by_chain, kc, 0) + 1
            seenf = Set{String}()
            for f in m1
                fn = string(basename(string(f.file)), ":", f.func)
                fn in seenf || (push!(seenf, fn); by_func[fn] = get(by_func, fn, 0) + 1)
            end
        elseif m1 !== nothing || workers_busy > 0
            n_parallel += 1
        end
    end
    @printf(io, "Thread activity summary (%d default threads, %d ticks):\n", nt, n_ticks)
    @printf(io, "  Main-thread-only busy ticks: %d (%.1f %% of ticks)\n", n_main_only, 100 * n_main_only / max(n_ticks, 1))
    @printf(io, "  Parallel (worker-active) ticks: %d (%.1f %% of ticks)\n", n_parallel, 100 * n_parallel / max(n_ticks, 1))
    @printf(io, "  Average worker utilization: %.1f %%\n\n", 100 * n_worker_busy / max(n_ticks * max(nt - 1, 1), 1))
    println(io, "Main-thread-only samples by inclusive repository function (present on stack):")
    for (k, v) in sort(collect(by_func); by = last, rev = true)[1:min(40, end)]
        @printf(io, "%6d  %5.1f %%  %s\n", v, 100v / max(n_main_only, 1), k)
    end
    println(io, "\nMain-thread-only samples by innermost repository frame:")
    for (k, v) in sort(collect(by_inner); by = last, rev = true)[1:min(40, end)]
        @printf(io, "%6d  %5.1f %%  %s\n", v, 100v / max(n_main_only, 1), k)
    end
    println(io, "\nMain-thread-only samples by 4 innermost repository frames:")
    for (k, v) in sort(collect(by_chain); by = last, rev = true)[1:min(50, end)]
        @printf(io, "%6d  %5.1f %%  %s\n", v, 100v / max(n_main_only, 1), k)
    end
end
function prof_sample(js)
    if js == PROFSAMPLE
        Base.invokelatest(() -> (Profile.init(n = 10^7, delay = 0.003); Profile.clear(); PROFALLOC > 0 ? (Profile.Allocs.clear(); Profile.Allocs.start(sample_rate = PROFALLOC)) : Profile.start_timer()))
    elseif js == PROFSAMPLE + PROFSAMPLEN
        PROFALLOC > 0 && return Base.invokelatest() do
            Profile.Allocs.stop()
            open(io -> alloc_report(io, PROFALLOC), joinpath(get(ENV, "POPIII_OUT", "."), "cplallocs_$(js).txt"), "w")
        end
        Base.invokelatest() do
            Profile.stop_timer()
            open(io -> mainthread_report(io), joinpath(get(ENV, "POPIII_OUT", "."), "cplprofile_mainonly_$(js).txt"), "w")
            open(joinpath(get(ENV, "POPIII_OUT", "."), "cplprofile_$(js).txt"), "w") do io
                Profile.print(IOContext(io, :displaysize => (100000, 250)); format = :flat, sortedby = :count, mincount = 50)
            end
            # (PERF512) also the call tree, all threads merged, with the module-qualified frames of the repository code
            open(joinpath(get(ENV, "POPIII_OUT", "."), "cplprofile_tree_$(js).txt"), "w") do io
                Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :tree, mincount = 100, maxdepth = 60, noisefloor = 2.0, C = false, combine = true)
            end
        end
    end
end
function prof_report()
    js = Int(TF.getscalar(:jstep)); PROFSAMPLE >= 0 && prof_sample(js)
    (CPLPROF && js % CPLPROFEVERY == 0 && js != PROFSTEP[]) || return nothing
    PROFSTEP[] = js
    @printf("CPLPROF step %d:", js); for s_ in (:matgen, :matslv, :maxdel, :update, :peq, :timstp); @printf(" %s %.1f s/%d %.1f GB", s_, get(PROFT, s_, 0.0), get(PROFN, s_, 0), get(PROFB, s_, 0.0) / 1e9); end
    empty!(PROFB)
    @printf("  cpl_zone jacobians %d (%.1f s), reused %d\n", PM.CPLPN[], PM.CPLPT[] * 1e-9, PM.CPLPH[]); flush(stdout)
    empty!(PROFT); empty!(PROFN); PM.CPLPN[] = 0; PM.CPLPH[] = 0; PM.CPLPT[] = 0
    nothing
end

# --- the initial-condition builder's helper functions (were defined inside the script's IC blocks; they use only globals, and iceulbuild! calls them) ---
# monotone cubic Hermite (Fritsch–Carlson) through the read-off points: no kinks (a piecewise-linear log n(r) put
# spurious pressure jumps of up to 40× gravity at the knots, M2e)
function pchip(xs, ys)
    n = length(xs); d = diff(ys) ./ diff(xs); m = zeros(n)
    m[1] = d[1]; m[n] = d[n-1]
    for i in 2:n-1; m[i] = d[i-1] * d[i] <= 0 ? 0.0 : (d[i-1] + d[i]) / 2; end
    for i in 1:n-1
        d[i] == 0 && (m[i] = 0.0; m[i+1] = 0.0; continue)
        a, b = m[i] / d[i], m[i+1] / d[i]; h = a^2 + b^2
        h > 9 && (t = 3 / sqrt(h); m[i] = t * a * d[i]; m[i+1] = t * b * d[i])
    end
    x -> begin
        x <= xs[1] && return ys[1] + m[1] * (x - xs[1]) * 0   # constant outside
        x >= xs[n] && return ys[n] + m[n] * (x - xs[n])        # linear (log-log power law) beyond the last point
        j = searchsortedlast(xs, x); h_ = xs[j+1] - xs[j]; t = (x - xs[j]) / h_
        (2t^3 - 3t^2 + 1) * ys[j] + (t^3 - 2t^2 + t) * h_ * m[j] + (-2t^3 + 3t^2) * ys[j+1] + (t^3 - t^2) * h_ * m[j+1]
    end
end
interp(x, xs, ys) = pchip(xs, ys)(x)
Tfin(m) = 10^Tfin_(log10(max(m, mcore_fin)))                         # m in M☉ (constant below 0.01 M☉)
xefin(m) = xefin_flat ? xefin_0 : 10^(let x = log10(m); x <= lmxe_fin[1] ? lxe_fin[1] : x >= lmxe_fin[end] ? lxe_fin[end] :
    (j = searchsortedlast(lmxe_fin, x); lxe_fin[j] + (x - lmxe_fin[j]) / (lmxe_fin[j+1] - lmxe_fin[j]) * (lxe_fin[j+1] - lxe_fin[j])) end)
function finy(m)                                                     # (x_H⁺, x_H⁻, x_H₂, x_H₂⁺) per H nucleus at m [M☉]
    m = max(m, mcore_fin); f2 = 10^fH2fin_(log10(m)); f1 = 10^fHIfin_(log10(m))
    (xefin(m), 0.0, 0.5 * f2 / (f2 + f1), 0.0)
end
vfin(m, r) = m >= mv1_fin ? 1e5 * vfin_(log10(m)) : 1e5 * vv1_fin * r / rM1fin   # [cm/s]; m [M☉], r [cm]
ytup(d) = (d[:Hp], d[:Hm], d[:H2], d[:H2p])
xinit(y) = [s === :Hp ? y[1] : s === :Hm ? y[2] : s === :H2 ? y[3] : s === :H2p ? y[4] : s === :Lip ? 0.0 : 0.0 for s in netw.free]

# the network material's cooling and background radiation (were closures of the script: named here so that the specializations of the material code
# on them are part of the image; they read the globals tcmb and cf exactly as the closures did)
Bν_(ν, T) = 2ChemNet.hP * ν^3 / ChemNet.c^2 / expm1(ChemNet.hP * ν / (ChemNet.kB * T))
J_cmb(ν) = Bν_(ν, tcmb)
J_zero(ν) = 0.0
coolnet(ρ, T, d) = cf(ρ, T, d === nothing ? PopIIIEqEOS.composition(ρ, T) : ytup(d))

# --- HOT=1 (S6-C1): the stellar core's helpers (were defined by the script's HOT block; the HOT key is guarded, so these exist exactly when the image was built with it) ---
if HOTON
    const SM = StarMaterial; const NN = NucNet; const SE = StellarEOS
    const HX = 1 / (1 + 4PC.yHe)
    const LASTSPLIT = Ref((NaN, 0, NaN))
    "e_gas of the blended material for an explicit stellar composition (the zone's, or a trial one in mix!)."
    function e_blend(ρ, T, comp)
        h = PM.HOT[]; T <= h.T1 && return PopIIIEqEOS.eos(ρ, T)[1]
        es = SE.gas_state(ρ, T, comp).e - h.Δ0; T >= h.T2 && return es
        u = PopIIIEqEOS.eos(ρ, T)[1]; u + PM.hotw(T, h.T1, h.T2) * (es - u)
    end
    "T with e_blend(ρ, T, comp) = e (safeguarded Newton in ln T), and the relative error of e."
    function T_blend(ρ, e, T0, comp)
        lT = log(T0); lo = lT - 1.5; hi = lT + 1.5; f = Inf
        for _ in 1:60
            e0 = e_blend(ρ, exp(lT), comp); f = e0 - e
            abs(f) <= 1e-12 * abs(e) && break
            f > 0 ? (hi = lT) : (lo = lT)
            d = (e_blend(ρ, exp(lT + 1e-6), comp) - e0) / 1e-6
            lN = lT - f / d; lT = (d > 0 && lo < lN < hi) ? lN : 0.5 * (lo + hi)
        end
        exp(lT), abs(f) / abs(e)
    end
    const HOTR = Ref{Any}((; Ebur = 0.0, corr = 0.0)); const HOTM = Ref{Any}((; nconv = 0)); const EBURN = Ref(0.0)
    "S6-C1: after a converged step and chem!: burn! and mix! of the hot zones (S5)."
    function hotburn!(dt)
        sg = SG[]
        if !any(TF.var(:tn)[k] > PM.HOT[].T1 for k in a:b); sg.acc.nstage = 0; return 0.0; end
        sd = TX.step_diag(dt)
        HOTR[] = SM.burn!(sg; sd); EBURN[] += HOTR[].Ebur
        HOTM[] = SM.mix!(sg; sd, capmix = false, bridge = 4, tmax = 0.02)
        dt > 0 ? HOTR[].Ebur / dt : 0.0      # a step of zero duration (Δt below the Float64 spacing of t) burns nothing
    end
end

# --- the GRIDEQ moving mesh (was grideq_driver.jl: its settings are assigned by that file at run time) ---
geqv(key, defs) = (v = haskey(ENV, key) ? parse.(Float64, split(ENV[key], ",")) : Float64[]; [i <= length(v) ? v[i] : defs[i] for i in eachindex(defs)])
const GEQJ0 = Ref(-1)                  # jstep at the first grideq_step! of this process (GEQRAMP; -2: no ramp, the mesh is close to the rule)
"max_j(M_j Δln r_j) / median over the zones above the centre: how far the mesh is from equidistribution (GEQRAMP)."
function geq_quality(s)
    D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL, mix = GEQDEC[3])
    q = [D.M[i] * log(s.r[i+2] / s.r[i+1]) for i in eachindex(D.M)]; qs = sort(q)
    (; ratio = maximum(q) / qs[(length(qs) + 1) ÷ 2], k = argmax(q) + 1)
end
const GEQPREV = Ref{Any}(nothing)      # (t, x, M) of the last step: the monitor function's own motion (feature tracking)
const GEQLAST = Ref{Any}(nothing)
"State on the present mesh (interfaces centre first; m per steradian)."
function geq_state()
    r = TF.var(:rn)[a:b+1]; u = TF.var(:un)[a:b+1]; m = envmenc()
    (; r, u, m, ρ = TF.var(:dn)[a:b], P = TF.var(:pgn)[a:b], T = TF.var(:tn)[a:b], E = TF.var(:ern)[a:b], cs = TF.var(:asn)[a:b], τ = gm_tau(a),
       X = GEQX[1] > 0 ? geq_lx() : nothing) |> geq_xd
end
"The extra logarithmic fields of the explicit rule (target, quality, feature speed): √wx ln(x + xfl) rows (GEQX) and √wχ ln(χρ) (GEQW entry 8)."
function geq_xd(s)
    rows = Vector{Vector{Float64}}()
    s.X === nothing || append!(rows, [sqrt(GEQX[1]) .* s.X[q, :] for q in axes(s.X, 1)])
    GEQW[8] > 0 && push!(rows, sqrt(GEQW[8]) .* log.(TF.var(:chifn)[a:b] .* s.ρ))
    merge(s, (; Xd = isempty(rows) ? nothing : permutedims(reduce(hcat, rows))))
end
"GEQX: ln(x_s + xfl) of every network species (rows) in the zones a…b (columns)."
function geq_lx()
    mat isa PM.NetGas || return nothing
    d = [ChemNet.full(mat.x[k], mat.net) for k in a:b]; sp = sort(collect(keys(d[1])))
    Float64[log(max(Float64(ForwardDiff.value(d[j][q])), 0.0) + GEQX[2]) for q in sp, j in eachindex(d)]
end
geq_target(s; N = nsh - 1) = GridEq.target(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ, N; w = Tuple(GEQW), gmax = GEQC[1], alph = GEQC[2], τlo = GEQC[6], τhi = GEQC[7],
                                            X = s.Xd, wx = 1.0, fl = GEQFL, mix = GEQDEC[3])
"Grid totals (erg, g): internal (e incl. the chemical offsets), kinetic (zone-centred u), radiation, gravitational (−G m dm / r_c)."
function geq_energy()
    r = TF.var(:rn)[a:b+1]; u = TF.var(:un)[a:b+1]; m = 4π .* envmenc(); e = TF.var(:egn)[a:b]; E = TF.var(:ern)[a:b]; n = nsh
    dm = diff(m); ei = sum(dm .* e); ek = sum(0.5 * dm[j] * (0.5 * (u[j] + u[j+1]))^2 for j in 1:n)
    er = sum(4π / 3 * (r[j+1]^3 - r[j]^3) * E[j] for j in 1:n); eg = -sum(6.674e-8 * 0.5 * (m[j] + m[j+1]) * dm[j] / (0.5 * (r[j] + r[j+1])) for j in 1:n)
    (; ei, ek, er, eg, M = m[end])
end
"Signal-crossing time of every zone, r_c/(|u| + c_s)."
geq_tcross(s) = [0.5 * (s.r[j] + s.r[j+1]) / (0.5 * (abs(s.u[j]) + abs(s.u[j+1])) + s.cs[j]) for j in eachindex(s.ρ)]
"""
Resolution diagnostics: zones across the accretion shock (from the interface of fastest infall inwards to where the infall has
dropped below 0.2 of it) and across the photosphere (interfaces with 0.1 ≤ τ ≤ 10), narrowest zone, largest neighbour ratio.
"""
function geq_diag(s)
    ku = argmin(s.u); nsk = 0
    if s.u[ku] < 0
        j = ku; while j > 1 && s.u[j-1] < 0.2 * s.u[ku]; j -= 1; end; nsk = ku - j
    end
    nph = count(t -> 0.1 <= t <= 10, s.τ); nph5 = count(t -> 1 <= t <= 1e5, s.τ)
    dT = maximum(abs(log(s.T[j+1] / s.T[j])) for j in 2:length(s.T)-1)      # largest jump of ln T across one interface
    dl = [log(s.r[j+1] / s.r[j]) for j in 2:length(s.ρ)]
    rat = maximum(max(dl[j] / dl[j-1], dl[j-1] / dl[j]) for j in 2:length(dl))
    Dn = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL, mix = GEQDEC[3])
    Cq = Dn.If; Ipq = Dn.Ip                             # I_f, I_p of the explicit monitors (guarantee: Δlog10 r <= I_f/(a N pr))
    (; nsk, rsk = s.r[ku], uk = s.u[ku], nph, nph5, dT, dlmin = minimum(dl), rat, Cq, Ipq)
end
function grideq_step!()
    s = geq_state(); n = nsh; dt = TF.getscalar(:dtime); tc = geq_tcross(s); κ = GEQC[3]; κf = GEQC[9]
    js_ = Int(TF.getscalar(:jstep))
    if GEQJ0[] == -1                                              # first call of this process: ramp only a mesh far from the rule (GEQRAMP)
        Q = geq_quality(s); f_ = GEQRAMP === nothing ? Inf : length(GEQRAMP) >= 3 ? GEQRAMP[3] : 10.0
        GEQJ0[] = Q.ratio > f_ ? js_ : -2
        GEQRAMP === nothing || pr("GEQRAMP: max M dln r / median %.3g (threshold %.3g) at zone %d: %s\n", Q.ratio, f_, Q.k,
                                  GEQJ0[] >= 0 ? @sprintf("kappa ramps from %.3g to %.3g over %d steps", GEQRAMP[1], κ, GEQRAMP[2]) : "no ramp")
    end
    (GEQRAMP === nothing || GEQJ0[] < 0) || (κ = GEQRAMP[1] * (κ / GEQRAMP[1])^clamp((js_ - GEQJ0[]) / GEQRAMP[2], 0.0, 1.0))
    fvr = (GEQRAMP !== nothing && length(GEQRAMP) >= 4) ? GEQRAMP[4] : Inf      # GEQRAMP entry 4: cap of the ramped grid speed in units of |u| + c_s
    tfe = fill(Inf, n); nfast = 0                                 # feature-crossing time of each zone (feature tracking)
    if κf > 0
        D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL, mix = GEQDEC[3])
        t = TF.getscalar(:timen); pv = GEQPREV[]
        if pv !== nothing && t > pv[1]
            Mp = [GridEq.lin(pv[2], pv[3], x) for x in D.x]; nd = length(D.x); v = zeros(nd)
            for i in 1:nd
                il, ih = max(i - 1, 1), min(i + 1, nd); Md = (D.M[ih] - D.M[il]) / (D.x[ih] - D.x[il])
                v[i] = -(D.M[i] - Mp[i]) / (t - pv[1]) * Md / (Md^2 + (2D.M[i])^2)
            end
            for i in 1:nd
                il, ih = max(i - 1, 1), min(i + 1, nd); j = i + 1      # (D covers zones 2…n)
                # Coherent translation filter: opposite-sign flanks around a stationary breathing peak (v[il]*v[ih] <= 0) cancel;
                # also bound by the local donor-preserving advective rate f_v * min(|u_j|, |u_{j+1}|) / r_c so feature tracking
                # never drives a quasi-hydrostatic front faster than the gas velocity |u| (which would flip the upwind donor u - rdot).
                vcoh = v[il] * v[ih] > 0 ? abs(sum(v[il:ih]) / (ih - il + 1)) : 0.0
                rc = 0.5 * (s.r[j] + s.r[j+1])
                vadv = GEQC[4] * min(abs(s.u[j]), abs(s.u[j+1])) / max(rc, 1e-300)
                vm = min(vcoh, vadv)
                vm > 0 && (tfe[j] = log(s.r[j+1] / s.r[j]) / vm)
            end
        end
        GEQPREV[] = (t, D.x, D.M)
    end
    if GRIDEQP == 1
        prev_tv = (TX.GRIDEQ[] !== nothing && hasproperty(TX.GRIDEQ[], :tauv) && length(TX.GRIDEQ[].tauv) == length(TF.var(:rn))) ? TX.GRIDEQ[].tauv : nothing
        tauv = zeros(length(TF.var(:rn))); tfl = zeros(length(TF.var(:rn))); sx = zeros(length(TF.var(:rn)))
        for j in 1:n
            # Implicit-grid relaxation time (mode 1): tau_j = smax(smin(kappa t_cross, f_dyn t_dyn), tfl) is formed in TitanResidFast.grideq_rows!
            # from the Newton iterate. Held here: kappa (GEQRAMP), the speed floor tfl (start-of-step, as before) and tauv = kappa t_cross (diagnostic).
            # The explicit feature tracker (kappa_f) and the 0.5 tau_old ratchet no longer enter: both are lagged and keep the mesh from following a collapse.
            dxj = j >= 2 ? log(s.r[j+1] / s.r[j]) : 1.0
            tfl[a+j-1] = dxj * tc[j] / min(GEQC[4] * GEQC[3] / κ, fvr); tauv[a+j-1] = κ * tc[j]
            jm, jp = max(j - 1, 2), min(j + 1, n); h = jp - jm == 2 ? 0.5 : 1.0
            sx[a+j-1] = 0.0                  # (the optical-depth term w_tau is formed from the Newton iterate in grideq_rows!)
            # velocity jump in units of the signal speed, (r ∂u/∂r)/(|u| + c_s): resolves shocks (the zone's own Δu, per zone index)
            GEQW[7] > 0 && (sx[a+j-1] += GEQW[7] * ((s.u[j+1] - s.u[j]) / (0.5 * (abs(s.u[j]) + abs(s.u[j+1])) + s.cs[j]))^2)
            # GEQX: composition, ln(x + xfl) of every species, adjacent interface jumps
            s.X === nothing || (sx[a+j-1] += GEQX[1] * h * sum((s.X[q, jp] - s.X[q, j])^2 + (s.X[q, j] - s.X[q, jm])^2 for q in axes(s.X, 1)))
        end
        TX.GRIDEQ[] = (; k1 = a + 1, k2 = b + 1, w = (GEQW[1], GEQW[2], GEQW[3], GEQW[4], GEQW[5]), alph = GEQC[2], gmax = GEQC[1], tauv, tfl, flc = 1 / min(GEQC[4] * GEQC[3] / κ, fvr), wtau = GEQW[6], tlo = GEQC[6], thi = GEQC[7], kap = κ, fdyn = GEQC[10], sx, wχ = GEQW[8], gp = (GEQGP[1], GEQGP[2] * 3.828e33 / (4π)), fr = GEQFL[1], fm = GEQFL[2], amix = GEQDEC[3])
    else
        rt = geq_target(s); w = TX.GRIDVEL[]; fill!(w, 0.0)
        for j in 3:n                                              # interfaces above the first one, below the outer edge
            tci = 0.5 * (tc[j-1] + tc[j]); trel = max(3dt, κ * tci)
            wk = (rt[j-1] - s.r[j]) / trel
            vm = GEQC[4] * (abs(s.u[j]) + 0.5 * (s.cs[j-1] + s.cs[j])); abs(wk) > vm && (wk = sign(wk) * vm)
            wd = wk > 0 ? s.r[j+1] - s.r[j] : s.r[j] - s.r[j-1]; q = abs(wk) * dt / wd; q > GEQC[5] && (wk *= GEQC[5] / q)
            w[a+j-1] = isfinite(wk) ? wk : 0.0
        end
        for k in 2:n-1                                            # no zone shrinks below half its narrower neighbour
            d = s.r[k+1] - s.r[k]; dn = d + (w[a+k] - w[a+k-1]) * dt; nb = min(s.r[k] - s.r[k-1], s.r[k+2] - s.r[k+1]) / 2
            if dn < d && dn < nb
                need = (min(d, nb) - dn) / dt; ml, mr = k >= 3, k + 1 <= n
                ml && mr ? (w[a+k-1] -= need / 2; w[a+k] += need / 2) : ml ? (w[a+k-1] -= need) : mr ? (w[a+k] += need) : nothing
            end
        end
    end
    GEQLAST[] = merge(geq_diag(s), (; nfast, κ))
    nothing
end
function grideq_line()
    d = GEQLAST[]; d === nothing && return ""
    rn = TF.var(:rn); ro = TF.var(:ro); un = TF.var(:un); cs = TF.var(:asn); dt = TF.getscalar(:timen) - TF.getscalar(:timeo); vmax = 0.0
    if dt > 0
        for k in a+2:b
            v = abs(rn[k] - ro[k]) / dt / (abs(un[k]) + 0.5 * (cs[k-1] + cs[k])); vmax = max(vmax, v)
        end
    end
    @sprintf("        grid eq: shock at %.4g R_sun (u %.1f km/s) across %d zones, photosphere (0.1 <= tau <= 10) across %d interfaces (1 <= tau <= 1e5: %d), largest |dln T| per interface %.3f; dln r min %.3e, max neighbour ratio %.2f; grid speed max %.3f (|u| + c_s); %d zones relaxing at their feature speed; kappa %.3g; I_f %.3g I_p %.3g (a %.2f)\n",
             d.rsk / Rsun, d.uk / 1e5, d.nsk, d.nph, d.nph5, d.dT, d.dlmin, d.rat, vmax, d.nfast, d.κ, d.Cq, d.Ipq, GEQDEC[3])
end

"""The step-loop output callback of `TD.run_julia!` (was an anonymous function in the driver script; it uses only globals)."""
function output_step(k0)
        k = k0 + KOFF[]; KLAST[] = k
        tnow = time(); TSTEP[] += tnow - TMARK[] - TOUT[]; TOUT[] = 0.0
        if haskey(ENV, "ZDIAG")                     # S6: energy terms of the hydro step per gas mass in zones z1..z2 (before chem!)
            zz = parse.(Int, split(ENV["ZDIAG"], ",")); dt_ = TF.getscalar(:timen) - TF.getscalar(:timeo)
            for z in zz[1]:zz[2]
                ρz = TF.var(:dn)[z]; Tz = TF.var(:tn)[z]; kp = TF.var(:xkpn)[z]; Ez = TF.var(:ern)[z]
                ke = TF.var(:xken)[z]
                xch = (ke * ChemNet.c * Ez - kp * 4 * 5.670374e-5 * Tz^4) * dt_            # (κ_E c E − κ_P 4σT⁴) Δt, the hydro exchange [erg/g]
                gq = isempty(TX.GASQ[]) ? 0.0 : TX.GASQ[][z] * dt_ / ρz; rq = isempty(TX.RAYQ[]) ? 0.0 : TX.RAYQ[][z] * dt_ / ρz
                pr("        ZDIAG %d: T %.1f (To %.1f) e %.5e (eo %.5e) rho %.3e (do %.3e) | per g over dt %.2e: kE cE-kP 4sT4 %.3e, GASQ %.3e, RAYQ %.3e, Lambda_ext %.3e, -p d(1/rho) %.3e; kP %.3e kE %.3e chi %.3e eq %s, cbar %.3e\n",
                   z, Tz, TF.var(:to)[z], TF.var(:egn)[z], TF.var(:ego)[z], ρz, TF.var(:do)[z], dt_, xch, gq, rq, mat.Λext[z] * dt_ / ρz,
                   -TF.var(:pgn)[z] * (1 / ρz - 1 / TF.var(:do)[z]), kp, ke, TF.var(:chifn)[z], string(mat.eqzone[z]), isempty(TX.CHEMB[]) ? 0.0 : TX.CHEMB[][z] * dt_)
                if ENVEUL !== nothing && !isempty(XMEPREV[])     # Eulerian budget: Δ(M e) per g vs exchange + advection (F e_b through both faces) + p dV; rest = MLT, diffusion, viscosity, …
                    xn = TF.var(:xmen); xo = XMEPREV[]; Mn = xn[z] - xn[z+1]; Mo = xo[z] - xo[z+1]; eb = TF.var(:egb)
                    Fi = xn[z] - xo[z]; Fo = xn[z+1] - xo[z+1]
                    adv = (Fi * eb[z] - Fo * eb[z+1]) / Mn; dME = (Mn * TF.var(:egn)[z] - Mo * TF.var(:ego)[z]) / Mn
                    rr_ = TF.var(:rn); uu_ = TF.var(:un); pdv = -TF.var(:pgn)[z] * (rr_[z+1]^2 * uu_[z+1] - rr_[z]^2 * uu_[z]) * dt_ / Mn   # (fixed mesh: p ∇·u, new u)
                    pr("        ZBUD %d: d(Me)/M %.3e = exch %.3e + adv %.3e (e_b in %.4e out %.4e, e %.4e, F in/M %.2e out/M %.2e) + GASQ %.3e + pdV %.3e + rest %.3e | x_H2 %.3e x_p %.3e\n",
                       z, dME, xch, adv, eb[z], eb[z+1], TF.var(:egn)[z], Fi / Mn, Fo / Mn, gq, pdv, dME - xch - adv - gq - pdv, ChemNet.full(mat.x[z], mat.net)[:H2], ChemNet.full(mat.x[z], mat.net)[:Hp])
                end
            end
        end
        tc0 = time(); Tb4 = CHEMDIAG ? copy(TF.var(:tn)) : nothing; ub4 = CHEMDIAG ? copy(TF.var(:egn)) : nothing; ENVEUL !== nothing && advect_species!()
        EB = isempty(TX.CHEMB[]) ? nothing : copy(TF.var(:egn))
        SLZ_ = haskey(ENV, "SLZ") ? parse.(Int, split(ENV["SLZ"], ",")) : nothing
        if SLZ_ !== nothing; xb_ = [copy(mat.x[z]) for z in SLZ_[1]:SLZ_[2]]; Tb_ = TF.var(:tn)[SLZ_[1]:SLZ_[2]]; end
        FLUSHOUT && flush(stdout)                    # (FLUSHOUT=1: the step's diagnostics reach the output file now)
        NANSCAN && nanscan("after the hydro step", k)
        dT = PM.cpl_chem!(mat)
        NANSCAN && nanscan("after the chemistry", k)
        prof_report()
        GRIDEQP === nothing || grideq_step!()
        if SLZ_ !== nothing && PM.SLDIAG[] !== nothing      # S6: starlight and chemistry of zones z1..z2 over this step's chem! (molecular-edge onset)
            sd = PM.SLDIAG[]; dtq = TF.getscalar(:timen) - TF.getscalar(:timeo)
            pr("        SLZ k %d dt %.3e: photosphere zone %d T_eff %.0f\n", k, dtq, sd.kph, sd.Te)
            ld_ = PM.LYADIAG[]
            ld_ === nothing || pr("        SLZ Lya: zones %d..%d (%d layers; photosphere rebuilt in zones %s, top T %.0f), %d bins, %d ALI iterations (dS %.1e), %.2f s\n", ld_.kin, ld_.kin + ld_.n - 1, ld_.nlayer, string(ld_.ksub), ld_.Tsurf, ld_.nb, ld_.it, ld_.δ, ld_.time)
            for (i, z) in enumerate(SLZ_[1]:SLZ_[2])
                d0_ = ChemNet.full(xb_[i], mat.net); d1_ = ChemNet.full(mat.x[z], mat.net)
                pr("        SLZ %d: n %.3e T %.1f -> %.1f  xH2 %.3e -> %.3e  xp %.3e -> %.3e  xHm %.2e | kLW %.3e fsh %.3e kHm %.3e eq %s E_r %.3e\n", z, TF.var(:dn)[z] / PC.μH,
                   Tb_[i], TF.var(:tn)[z], d0_[:H2], d1_[:H2], d0_[:Hp], d1_[:Hp], d1_[:Hm], sd.kLW[z], sd.fsh[z], sd.kHm[z], string(mat.eqzone[z]), TF.var(:ern)[z])
                t3_ = sd.tl3[z]
                (t3_ !== nothing && hasproperty(t3_, :qbb)) && pr("        SLZ %d tl3: T_c %.0f P12 %.2e RL %.2e beta %.2e ne %.2e | chemq %.3e cycq %.3e Lya-bb %.3e erg/cm3/s, Lya-bb dt/(rho e) %.2e\n", z,
                    t3_.Tc, t3_.P12, t3_.RL, t3_.β, t3_.ne, t3_.chemq, t3_.cycq, t3_.qbb, t3_.qbb * dtq / (TF.var(:dn)[z] * TF.var(:egn)[z]))
                if ld_ !== nothing && t3_ !== nothing && ld_.kin <= z < ld_.kin + ld_.n && !isnan(ld_.nJ[z-ld_.kin+1])
                    i_ = z - ld_.kin + 1; Tx_(nb) = nb > 0 ? 10.19881 * 1.602176634e-12 / (1.380649e-16 * log1p(1 / nb)) : 0.0
                    pr("        SLZ %d Lya: nbar(Jbar) %.3e (T_exc %.0f) nbar(S_L) %.3e (T_exc %.0f) nbar(J_ext) %.3e | R_Lya %.3e (Sobolev %.3e) P12 %.3e (Sobolev %.3e)\n", z,
                       ld_.nJ[i_], Tx_(ld_.nJ[i_]), ld_.nS[i_], Tx_(ld_.nS[i_]), ld_.Jext[i_], 3 * 6.2648e8 * ld_.oml[i_], ld_.Rsob[i_], t3_.P12, ld_.Psob[i_])
                end
            end
        end
        if EB !== nothing                            # BALSPLIT: c̄ ← c̄ + Δe_chem/Δt where the zone stayed on the network
            eg_ = TF.var(:egn); dtc_ = TF.getscalar(:timen) - TF.getscalar(:timeo); cb_ = TX.CHEMB[]; on_ = TX.CHEMBON[]
            for kk in eachindex(cb_)
                cb_[kk] = (on_[kk] && dtc_ > 0) ? cb_[kk] + (eg_[kk] - EB[kk]) / dtc_ : 0.0
            end
        end
        if CHEMDIAG && dT > 0.01                      # the zones chem! moved most
            tt_ = TF.var(:tn); jj = sort(collect(a:b); by = kk -> -abs(tt_[kk] / Tb4[kk] - 1))[1:3]
            for kk in jj
                pr("        chem! zone %d: T %.1f -> %.1f, m %.4e Msun, n %.3e, p %.3e, eq %s, e_before %.5e e_eos_after %.5e, dt %.2e\n", kk - a + 1, Tb4[kk], tt_[kk], 4π * xmi(kk+1) / Msun,
                   TF.var(:dn)[kk] / (1.22PC.mH), TF.var(:pgn)[kk], string(mat.eqzone[kk]), ub4[kk], TF.var(:egn)[kk], TF.getscalar(:timen) - TF.getscalar(:timeo))
            end
        end
        HOTON && (LNUC[] = hotburn!(TF.getscalar(:timen) - TF.getscalar(:timeo)))
        TCHEM[] += time() - tc0
        RAD && TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())   # T changed in chem!: Planck function, opacity, e_g + E_r/ρ
        REC.on && SNAPEVERY > 0 && (k <= 2 || k % SNAPEVERY == 0) && recsnap!(k)
        PM.MGON && SNAPEVERY > 0 && (k <= 2 || k % SNAPEVERY == 0) && mg_record(k)
        if SNAPS !== nothing                         # profiles for post-processing: every 0.25 dex of n_c, and every SNAPEVERY steps
            lnc = log10(TF.var(:dn)[a] / PC.mH)
            if lnc >= SNAPNEXT[] || (SNAPEVERY > 0 && (k <= 2 || k % SNAPEVERY == 0))
                push!(SNAPS, (; step = k, t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1],
                    ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], y = [ycomp(a + j - 1) for j in 1:nsh],
                    m = 4π .* (ENVEUL === nothing ? TX.XMINT[][a:b+1] : envmenc()), Λext = copy(mat.Λext[a:b]), eq = [eqflag(a + j - 1) for j in 1:nsh],
                    er = TF.var(:ern)[a:b], fr = TF.var(:frn)[a:b+1], pg = TF.var(:pgn)[a:b], chi = TF.var(:chifn)[a:b],
                    fedd = TF.var(:fedd)[a:b], qf = TF.var(:qf)[a:b], xp = [ycomp(a + j - 1)[1] for j in 1:nsh], xh2 = [ycomp(a + j - 1)[3] for j in 1:nsh], L_lines = LLINE[]))
                lnc >= SNAPNEXT[] && (SNAPNEXT[] = lnc + 0.25)
                serialize(SNAPFILE, SNAPS)
            end
        end
        if cool == "lines"
            # H₂ line transfer on the new state (RHFC02 A2.4); Λ held fixed over the next step; Δt ≤ 0.1 u ρ/|Λ| (A2.2)
            rr = TF.var(:rn)[a:b+1]; uu = TF.var(:un)[a:b+1]; dd = TF.var(:dn)[a:b]; tt = TF.var(:tn)[a:b]
            nh2 = [dd[j] / PC.μH * ycomp(a + j - 1)[3] for j in 1:nsh]
            tl0 = time(); Λl, Lr, _ = linecool(rr, uu, dd, tt, nh2, k); TLINE[] += time() - tl0
            tcmin = Inf
            for j in 1:nsh
                mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                ue = 1.5 * TF.var(:pgn)[a+j-1] / dd[j]          # thermal energy scale of the cooling time (diagnostic: the line cap of DTDIAG and DTERR's trace)
                Λl[j] != 0 && abs(ue * dd[j] / Λl[j]) < tcmin && (tcmin = abs(ue * dd[j] / Λl[j]); JCMIN[] = j)
            end
            LLINE[] = Lr[end]
            if get(ENV, "DTDIAG", "0") == "1"
                j_ = JCMIN[]; pr("        line cap: t_cool %.3e s at zone %d (m %.4e Msun, T %.1f, n %.3e, e+U %.3e, Lambda %.3e); dtime before %.3e, smax %.3e\n", tcmin, j_, 4π * xmi(a+j_) / Msun, tt[j_], dd[j_] / (1.22PC.mH),
                                TF.var(:egn)[a+j_-1] + PC.UOFF_TH, Λl[j_], TF.getscalar(:dtime), TF.getscalar(:smax))
                pr("        LINEIMP: %d zone line models (ALI rates); no t_cool cap\n", count(!isnothing, PM.LIMODEL[]))
            end
            LINECAPV[] = tcmin                       # (diagnostic only: DTERR sets Δt; LINEIMP has no t_cool cap)
        end
        if get(ENV, "STEPKX", "0") == "1"             # which equation and zone set the step's relative change (timstp!)
            pr("STEP %d dt %.3e smax %.3e: %s\n", k, TF.getscalar(:dtime), TF.getscalar(:smax),
               join(("$(TF.var(:kx)[i]):$(round(TF.var(:dx)[i], sigdigits = 2))" for i in 1:Int(TF.getscalar(:neqn))), " "))
        end
        if RAD && haskey(ENV, "ERDIAG") && k in parse.(Int, split(ENV["ERDIAG"], ","))     # radiation field across the grid
            ern = TF.var(:ern); frn = TF.var(:frn); fe = TF.var(:fedd); ch = TF.var(:chifn); dd = TF.var(:dn); tt = TF.var(:tn); rr = TF.var(:rn)
            Tb = parse(Float64, get(ENV, "TCMB", "0")); τ = 0.0; τs = zeros(b)
            for kk in b:-1:a; τ += ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]); τs[kk] = τ; end
            pr("  radiation at step %d (geddr %.3f):\n", k, TF.getscalar(:geddr))
            for kk in (haskey(ENV, "ERDIAGZ") ? (zr = parse.(Int, split(ENV["ERDIAGZ"], ",")); (a - 1 + zr[1]):zr[2]:(a - 1 + zr[3])) : (a:6:b))     # ERDIAGZ=lo,step,hi (zones)
                pr("    z %3d m %.3e r %.3e tau %.3e dtau %.2e E/aT4 %.3e E/aTcmb4 %.3e F/cE %+.3e f %.3f T %.1f\n", kk - a + 1, 4π * xmi(kk+1) / Msun, rr[kk+1],
                   τs[kk], ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]), ern[kk] / (Z0.arad * tt[kk]^4), ern[kk] / (Z0.arad * Tb^4), frn[kk+1] / (Z0.c * ern[kk]), fe[kk], tt[kk])
            end
        end
        if haskey(ENV, "PROFAT") && k in parse.(Int, split(ENV["PROFAT"], ","))
            rr = TF.var(:rn); uu = TF.var(:un); dd = TF.var(:dn); tt = TF.var(:tn)
            pr("  profile at step %d, t = %.4e yr:\n", k, TF.getscalar(:timen) / 3.156e7)
            for j in 1:8:nsh
                kk = a + j - 1
                pr("    j %3d m %.3e r %.3e pc  n %.3e  T %7.1f  u %7.3f km/s  f_H2 %.2e\n", j, 4π * xmi(kk+1) / Msun, rr[kk+1] / 3.0857e18,
                   dd[kk] / (1.22 * PC.mH), tt[kk], uu[kk+1] / 1e5, 2ycomp(kk)[3])
            end
        end
        ρc_ = TF.var(:dn)[a]; Tc_ = TF.var(:tn)[a]; ln = log10(ρc_ / PC.mH)
        if ln >= nextln[] || k % PRINTEVERY == 0
            umin = minimum(TF.var(:un)[a:b+1]) / 1e5
            pr("%7d %8.3f %10.4e %9.1f %9.1f %10.3e %10.3f %9.2e %9.2e%s\n", k, ln, TF.getscalar(:timen) / 3.156e7, Tc_, on98T(ln),
               2ycomp(a)[3], umin, dT, TF.getscalar(:dtime), cool == "lines" ? @sprintf("  L_lines %.3e  bins %d, thick lines %d", LLINE[], H2Lines.LAST[].nbins, H2Lines.LAST[].nthick) : "")
            ln >= nextln[] && (nextln[] += 0.25)
            if HOTON
                g_ = SG[]; tn_ = TF.var(:tn); un_ = TF.var(:un); jc_ = findfirst(<(-2e5), un_[a:b+1]); ms_ = jc_ === nothing ? NaN : 4π * xmi(a+jc_-1) / Msun; rs_ = jc_ === nothing ? NaN : TF.var(:rn)[a+jc_-1] / 6.957e10
                REC.mstar[] = ms_; REC.rstar[] = rs_
                pr("        hot: t %.6e yr, m* %.5g Msun, R* %.4g Rsun, L_nuc %.3e erg/s, X_c %.5f, D_c/D_0 %.3e, CNO_c %.2e, hot zones %d, conv %d, split %.1e/%d/%.1e, corr %.1e\n", TF.getscalar(:timen) / 3.156e7, ms_, rs_, LNUC[],
                   g_.Y[SM.IH, a], g_.Y[NN.IDX[:d], a] / (DH * HX), SM.xyz(view(g_.Y, :, a))[3], count(k -> tn_[k] > PM.HOT[].T1, a:b), HOTM[].nconv, LASTSPLIT[]..., HOTR[].corr)
            end
            if get(ENV, "UDIAG", "0") == "1"            # where the fastest infall is
                uu = TF.var(:un)[a:b+1]; ku = argmin(uu); rr = TF.var(:rn)[a:b+1]
                pr("        u_min %.3f km/s at interface %d of %d, r = %.3e pc, m = %.3e Msun, T(zone) = %.1f, n = %.3e\n", uu[ku] / 1e5, ku, nsh + 1,
                   rr[ku] / 3.0857e18, 4π * xmi(a+ku-1) / Msun, TF.var(:tn)[a+min(ku, nsh)-1], TF.var(:dn)[a+min(ku, nsh)-1] / (1.22 * PC.mH))
            end
            if get(ENV, "DTDIAG", "0") == "1"           # which variable/zone limits Δt (timstp: kx, dx per equation)
                kx = TF.var(:kx); dx = TF.var(:dx); names = ("r", "m", "rho", "u", "T", "E_r", "F_r"); ne_ = RAD ? 6 : 5   # (F_r is not in smax)
                j = argmax(abs.(dx[[1, 2, 3, 5, 6][1:ne_-1]])); j = (1, 2, 3, 5, 6)[j]; xs = TF.getscalar(:stol) / TF.getscalar(:smax)
                REC.dtlim[] = string(names[j], "@", kx[j])
                pr("        dt limit: %s in zone %d (rel. change %.2e, stol/smax = %.2f)%s\n", names[j], kx[j], dx[j], xs,
                   (RAD ? @sprintf("; E_r: zone %d %.2e, rho: zone %d %.2e", kx[6], dx[6], kx[3], dx[3]) : ""))
            end
        end
        if ENVEUL !== nothing
            envsink!()
            if k % PRINTEVERY == 0
                if TX.SLABFLUX[] && !isempty(TX.SFRATIO)   # SLABFLUX: the interfaces whose optical depth changed by more than x1.5 (zone numbers from the inner boundary)
                    jf = [j for j in a+1:min(b, length(TX.SFRATIO)) if !(1 / 1.5 < TX.SFRATIO[j] < 1.5)]
                    pr("        SLABFLUX: %d interface(s) with dtau/dtau_T outside [2/3, 3/2]%s\n", length(jf), join((@sprintf(" | %d: %.3g", j - a + 1, TX.SFRATIO[j]) for j in jf[1:min(end, 8)]), ""))
                end
                if TX.SUBZONE[] && !isempty(TX.SZOMEGA)    # SUBZONE: the interfaces with a front weight > 0.05 (zone numbers from the inner boundary)
                    jf = [j for j in a+1:min(b, length(TX.SZOMEGA)) if TX.SZOMEGA[j] > 0.05]
                    pr("        SUBZONE: %d front interface(s)%s\n", length(jf), join((@sprintf(" | %d: w %.3f dtau/dtau_T %.3f", j - a + 1, TX.SZOMEGA[j], TX.SZRATIO[j]) for j in jf), ""))
                end
                j_ = envsurface(); uu = TF.var(:un)[a:b+1]; ku = argmin(uu); rr = TF.var(:rn)[a:b+1]; mm = 4π .* envmenc() ./ Msun
                pr("        env: chem fallbacks %d, t %.6e yr, Mc+grid(<surf) %.4f Msun, sink %.2e Msun/yr, surface r %.1f R_sun, u_min %.1f km/s at r %.1f (v_esc %.1f), dt %.2e s\n", PM.NCHEMFB[], TF.getscalar(:timen) / 3.156e7, mm[j_],
                   4π * ENVMDOT[] / Msun * 3.156e7, rr[j_] / Rsun, uu[ku] / 1e5, rr[ku] / Rsun, sqrt(2Z0.G * mm[ku] * Msun / rr[ku]) / 1e5, TF.getscalar(:dtime))
                # radiation force on the infall: Γ = χF/c / (G m/r²) per zone outside 2 × the surface; Ṁ = −4π r² ρ u at the surface and at 1e5 R_sun
                let fr_ = TF.var(:frn)[a:b+1], ch_ = TF.var(:chifn)[a:b], dn_ = TF.var(:dn)[a:b], Γ = zeros(nsh)
                    for j in 1:nsh
                        rr[j] > 2rr[j_] || continue
                        rc = 0.5(rr[j] + rr[j+1]); Γ[j] = ch_[j] * 0.5(fr_[j] + fr_[j+1]) / ChemNet.c / (Z0.G * 0.5(mm[j] + mm[j+1]) * Msun / rc^2)
                    end
                    jg = argmax(Γ); w5 = findfirst(>=(1e5 * Rsun), rr); w4 = findfirst(>=(1e4 * Rsun), rr)
                    mdot(j) = j === nothing || j > nsh ? NaN : -4π * rr[j]^2 * 0.5(dn_[max(j-1, 1)] + dn_[min(j, nsh)]) * uu[j] / Msun * 3.156e7
                    pr("        feedback: Gamma_max %.3f at r %.0f R_sun; Gamma(1e4) %.3f, Gamma(1e5) %.3f; Mdot surface %.3e, 1e4 %.3e, 1e5 %.3e Msun/yr; xadv1 %d; to table %d, from table %d\n",
                       Γ[jg], 0.5(rr[jg] + rr[jg+1]) / Rsun, w4 === nothing || w4 > nsh ? NaN : Γ[w4], w5 === nothing || w5 > nsh ? NaN : Γ[w5],
                       mdot(j_), mdot(w4), mdot(w5), NXADV1[], PM.NSWIN[], PM.NSWOUT[])
                end
                GRIDEQP === nothing || pr("%s", grideq_line())
                # causality of the radiation field outside the surface: F/cE at the interfaces (FLUXTIE or not)
                let fz = TF.var(:frn), ez = TF.var(:ern), rz = TF.var(:rn), q = 0.0, kq = 0, k1 = 0
                    for k in a+1:b
                        rz[k] > 2rr[j_] || continue
                        v = fz[k] / (ChemNet.c * 0.5(ez[k-1] + ez[k])); v > q && (q = v; kq = k)
                    end
                    if TD.FLUXTIE[][2] > 0 && length(TD.TIEW) >= b
                        k1 = something(findfirst(k -> TD.TIEW[k] > 0, a+1:b), 0); k1 > 0 && (k1 += a)
                    end
                    pr("        causality: max F/cE %.3f at zone %d (r %.0f R_sun)%s\n", q, kq, kq > 0 ? rz[kq] / Rsun : NaN,
                       k1 > 0 ? @sprintf("; flux tie from zone %d (r %.0f R_sun, h %.3f)", k1, rz[k1] / Rsun, TD.HEDD[k1]) : "")
                end
                if PM.SLDIAG[] !== nothing
                    sd_ = PM.SLDIAG[]; rz_ = TF.var(:rn); xs_ = String[]
                    for lr in (3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
                        kk = findfirst(k -> rz_[k] / Rsun >= 10^lr, sd_.kph+1:length(sd_.kLW))
                        kk === nothing && continue; kk += sd_.kph
                        push!(xs_, @sprintf("r 1e%.1f: kHm %.1e kH %.1e kLW %.1e fsh %.1e", lr, sd_.kHm[kk], sd_.kH[kk], sd_.kLW[kk], sd_.fsh[kk]))
                    end
                    pr("        light: photosphere zone %d r %.0f R_sun, T_eff %.0f K; %s\n", sd_.kph, TF.var(:rn)[sd_.kph+1] / Rsun, sd_.Te, join(xs_, "; "))
                end
            end
        end
        (Tc_ > Tstop || ρc_ / PC.mH > nstop) && (pr("stopped: T_c %.4e K (TSTOP %.2e), n_c %.4e (NSTOP %.2e) at step %d\n", Tc_, Tstop, ρc_ / PC.mH, nstop, k); TF.setscalar!(:tmax, TF.getscalar(:timen)))
        TF.getscalar(:dtime) < DTSTOP && (pr("stopped: dt = %.2e s below DTSTOP at step %d\n", TF.getscalar(:dtime), k); TF.setscalar!(:tmax, TF.getscalar(:timen)))
        if TIMING && k % PRINTEVERY == 0
            pr("        timing since last report: step (hydro+radiation) %.1f s, chemistry %.1f s, H2 lines %.1f s, other output %.1f s\n",
               TSTEP[], TCHEM[], TLINE[], time() - tnow - 0.0 + TOTH[] - TCHEM[] - TLINE[])
            REC.on && recperf!(k, tnow)
            TSTEP[] = TCHEM[] = TLINE[] = TOTH[] = 0.0
        else
            TOTH[] += time() - tnow
        end
        TMARK[] = time()
        # CKDUMP=k1:k2[,k3:k4…] (diagnostic, default off): also write a full checkpoint <SAVE>_k<step> at each listed step (restart A/B tests)
        haskey(ENV, "CKDUMP") && haskey(ENV, "SAVE") && any(r_ -> (p_ = parse.(Int, split(r_, ":")); p_[1] <= k <= p_[end]), split(ENV["CKDUMP"], ",")) && checkpoint(ENV["SAVE"] * "_k" * string(k), k; jinc = 1)
        haskey(ENV, "SAVE") && CKEVERY > 0 && k % CKEVERY == 0 && checkpoint(ENV["SAVE"] * "_auto" * string((k ÷ CKEVERY) % 2), k; jinc = 1)
        haskey(ENV, "SAVE") && (LASTOK[] = TD.snapshot())     # the accepted state, for the checkpoint of a failed step

end

"The entry point: `args` are the command-line arguments (`cool [nsh]`); evaluates the run's statements (popiii_run_body.jl) in this module."
function main(args = ARGS)
    args === Base.ARGS || (empty!(Base.ARGS); append!(Base.ARGS, args))
    Base.include(@__MODULE__, joinpath(@__DIR__, "popiii_run_body.jl"))
    nothing
end

end # module
