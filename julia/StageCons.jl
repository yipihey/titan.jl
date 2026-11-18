"""
    StageCons

STAGECONS=n (popiii_m2c, diagnostic, default off = bit-identical): on every n-th accepted step (Trbdf2.ESTATS[:steps] % n == 0) print `STAGECONS` lines that test whether the
two TR-BDF2 stages and the FSAL derivative see the same right-hand side. Everything the diagnostic touches (NewtonKrylov.Snap: all COMMON arrays and scalars, closure arrays,
species iterate; plus the correction `rhs`) is restored, so the step proceeds exactly as without it.

All numbers are normalized rates [1/s]: ρ_row[k] = −s·R_row[k] / C[k], R the TRUE residual (TitanResidFast.residual! at the state, zero correction), C the matching conserved
combination of the start-of-step state yₙ (Trbdf2.combos), s the factor that turns the row into d C/dt. Every row is "(C_new − C_old)/dt + (fluxes, sources)" so that ρ = (dC/dt)/C the row
implies; compare with ε = (z2 − z1)/(yₙ Δ) of DTERRDIAG (z = Δ × derivative). Row (unknown it belongs to), C, s:
  r   (1)  C = r (:r)                 s = 1     (lgrid 3 row (rn − ro)/dt − u: ρ = u_θ − (rn − ro)/dt)
  m   (3)  continuity,  C = m = ρΔV (:m)           s = 1
  P   (4)  momentum,    C = u(m_{k−1} + m_k) (:P)  s = 2    (row is Δ[u(m_{k−1}+m_k)]/(2dt))
  E   (5)  energy row, s = 1, of a type that depends on the zone (printed as `Etype` at the start): `etot` (ETOTAL and not on the dual-energy switch) Δ[m(e_g + E_r/ρ) + m(u_k²+u_{k+1}²)/4]/dt, C = C_mer + K,
           K the kinetic energy of yₙ; `gas` (GASROW) Δ[m e_g]/dt, C = m e_g (:me); `mer` (neither) Δ[m(e_g + E_r/ρ)]/dt, C = :mer
  EV  (6)  E_r, C = E_r ΔV (:EV)      s = 1
  QF  (7)  F_r, C = F·2(m_{k−1}+m_k)/(ρ_{k−1}+ρ_k) (:QF)  s = 2c²   (row is ΔQF/(2c² dt); in FLUXTIE zones (TIEW > 0) the row is blended with the tie condition: not a derivative)
  x1, x2 (8, 9) species, C = m x (:x1, :x2)    s = 1
  xme (2)  mass-coordinate row x_me(k−1) − x_me(k) − m_{k−1}: no time derivative; printed R/m_{k−1} (dimensionless, NOT a rate)
Item 2 evaluates the energy row with etotal = false (all zones gas-only or mer type) and compares :me and :mer = :me + :EV with the FSAL D.

Items per printed step (jstep J): (1) at the start-of-step state with new := old (all primary unknowns and the species; update! with a zero correction then refreshes ghost
cells and derived arrays; no predict!) the residual at (θ, dt) = (½, γΔ), (1, γΔ), (1, C2Δ), (½, 10γΔ): with new = old all time differences vanish, so every row must be independent of
θ and dt; (2) −R(new = old) against the FSAL derivative DSAL of the previous step, (f − D)/C (energy from the row with etotal = false); (1b: θ-independence there, where the old-level face values enter; 2b) the same with r and x_me displaced by Δt × their FSAL velocities (the rows contain the mesh velocity (rn − ro)/dt and the mass flux dmdt = −(xmen − xmeo)/dt, which are 0 at new = old) and the time-difference part (C_new − C_old)/dt of each row removed; (3) the true
residual at each converged stage (TR: θ ½, dt γΔ, old = yₙ; BDF2: θ 1, dt C2Δ, old = y*), with Newton iterations and final dmax; (4) every old-level array the residual reads against its
new-level counterpart at the start state (max relative difference with zone, values in zones 31–33), and the closure arrays against their values at the end of the previous step; (6) the old state of the step against the state at the end of the last converged BDF2 stage (anything done to the state between the steps).
Zones are numbered from the inner boundary (zone Z = array index ngrs + Z − 1, as in DTERRDIAG). Skipped on a step without FSAL (first step after a restart: DSAL === nothing).
"""
module StageCons

using Printf
import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..NewtonKrylov, ..Trbdf2
const TC = TitanCore; const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const NK = NewtonKrylov; const TR = Trbdf2
vf(s) = TF.var(s)::Vector{Float64}
const GATE = Dict{Symbol,Tuple{Float64,Symbol,Int}}()      # per printed step: worst value, row, zone of item 1 / item 1dt / item 2b (summary line `gate`)
gate!(k, v, row, z) = (!haskey(GATE, k) || v > GATE[k][1]) && (GATE[k] = (v, row, z))
const EVERY = Ref(0); const ZS = 28:34; const ZV = 31:33
const ENDCL = Ref{Any}(nothing)       # closure arrays at the end of the last converged BDF2 stage
const ENDST = Ref{Any}(nothing)       # primary unknowns and species at the end of the last converged BDF2 stage

const PRIM = (:rn => :ro, :xmen => :xmeo, :dn => :do, :un => :uo, :tn => :to, :ern => :ero, :frn => :fro)
endstate() = (; prim = Dict(n_ => copy(vf(n_)) for (n_, _) in PRIM), xs = [copy(TX.XSN[][i]) for i in 1:2])
closure() = (; fedd = copy(vf(:fedd)), hedd = copy(TD.HEDD), hxedd = copy(TD.HXEDD), fxedd = copy(TD.FXEDD), tiew = copy(TD.TIEW))

function evalR(θ, dt; etotal = TX.ETOTAL[])
    TD.s!(:thet, θ); TD.s!(:dtime, dt)
    c = TX.context(; etotal); ws = TX.getws(Float64, c.mgr); R = zeros(c.neqn, c.nb)
    TX.residual!(R, ws, c, zeros(c.neqn * c.nb)); (R, c, ws.a)
end

"Normalized rates ρ_row[k] (NaN where undefined) of a residual, see the module docstring; :xme is R/m_{k−1}."
function rates(R, c, Cn; Cnew = nothing, dt = 1.0)
    ir, im, id, iu, it, ie, iF = c.ix; a, b, k0 = c.ngrs, c.ngre, c.k0; cc = Float64(TD.g(:cc)); n = length(Cn[:m]); out = Dict{Symbol,Vector{Float64}}()
    # with Cnew (the combos of the displaced state of item 2b): ρ = (−s R + (C_new − C_old)/dt)/C, i.e. the flux part of the row only
    function put(name, row, s, C, ks, Cnw = nothing)
        v = fill(NaN, n); for k in ks; C[k] != 0 && (v[k] = (-s * R[row, k-k0+1] + (Cnw === nothing ? 0.0 : (Cnw[k] - C[k]) / dt)) / C[k]); end; out[name] = v
    end
    nw(key) = Cnew === nothing ? nothing : Cnew[key]
    put(:r, ir, 1.0, Cn[:r], a:b+1); put(:m, id, 1.0, Cn[:m], a:b, nw(:m)); put(:P, iu, 2.0, Cn[:P], a+1:b, nw(:P))
    ty = etypes(c); uo = TR.SAVE[:uo]
    CE = zeros(n); for k in a:b; CE[k] = ty[k] === :etot ? (c.lrad ? Cn[:mer][k] : Cn[:me][k]) + 0.25Cn[:m][k] * (uo[k]^2 + uo[k+1]^2) : ty[k] === :gas ? Cn[:me][k] : Cn[:mer][k]; end
    CEn = Cnew === nothing ? nothing : [ty[k] === :gas ? Cnew[:me][k] : ty[k] === :mer ? Cnew[:mer][k] : NaN for k in 1:n]
    put(:E, it, 1.0, CE, a:b, CEn)
    c.lrad && (put(:EV, ie, 1.0, Cn[:EV], a:b, nw(:EV)); put(:QF, iF, 2cc^2, Cn[:QF], a+1:b, nw(:QF)))
    if TX.COUPLED[]
        zs = [k for k in a:b if TX.cplz(k)]; put(:x1, 8, 1.0, Cn[:x1], zs, nw(:x1)); put(:x2, 9, 1.0, Cn[:x2], zs, nw(:x2))
    end
    v = fill(NaN, n); for k in a+1:b+1; v[k] = R[im, k-k0+1] / Cn[:m][k-1]; end; out[:xme] = v
    out
end
const ROWS = (:r, :m, :P, :E, :EV, :QF, :x1, :x2, :xme)
"Row type of the energy row per zone (see the module docstring)."
etypes(c) = [(c.etotal && !c.useint[k]) ? :etot : (TX.ECOMPAT[] || !c.lrad || TX.GASROW[]) ? :gas : :mer for k in 1:length(c.useint)]

function maxabs(v, a)
    best = 0.0; kb = 0
    for k in eachindex(v); x = v[k]; (isfinite(x) && abs(x) > best) && (best = abs(x); kb = k); end
    (best, kb == 0 ? 0 : kb - a + 1)
end
zl(v, a, Zs = ZS) = "[" * join((@sprintf("%.2e", v[a+Z-1]) for Z in Zs), " ") * "]"

function item1(J, Δ, Cn, a)
    γ = TR.γ; cases = ((0.5, γ * Δ), (1.0, γ * Δ), (1.0, TR.C2 * Δ), (0.5, 10γ * Δ)); ρ = Dict{Int,Any}(); xkn = Float64[]; dn_ = Float64[]
    for (i, (θ, dt)) in enumerate(cases)
        R, c, W = evalR(θ, dt); ρ[i] = rates(R, c, Cn)
        i == 1 && item4(J, c, W, Cn, a)
        i == 2 && (xkn = copy(W.xken); dn_ = copy(W.dn))
    end
    for s in ROWS
        haskey(ρ[1], s) || continue
        h, o, o2, h10 = ρ[1][s], ρ[2][s], ρ[3][s], ρ[4][s]
        m1, z1 = maxabs(h .- o, a); m2, z2 = maxabs(h .- h10, a); m3, z3 = maxabs(o .- o2, a)
        gate!(:item1, m1, s, z1); gate!(:item1dt, max(m2, m3), s, m2 >= m3 ? z2 : z3)
        @printf("STAGECONS jstep %d item1 %-4s max|R(½,γΔ)-R(1,γΔ)| %.3e zone %d | max|R(½,γΔ)-R(½,10γΔ)| %.3e zone %d | max|R(1,γΔ)-R(1,C2Δ)| %.3e zone %d | zones %d..%d R½ %s R1 %s diff %s\n",
                J, s, m1, z1, m2, z2, m3, z3, ZS[1], ZS[end], zl(h, a), zl(o, a), zl(h .- o, a))
    end
    if TX.EXCAP[] > 0       # item 7: the explicit dt of the exchange factor `excap` (TitanResidFast.excap) at dt = γΔ and C2Δ, and the E_r-row difference it makes there
        b = Int(TD.g(:ngre)); cc = Float64(TD.g(:cc)); ex(dt) = [TX.EXCAP[] > 0 ? min(1.0, TX.EXCAP[] / (cc * xkn[k] * dn_[k] * dt)) : 1.0 for k in a:b]
        eg_, ec_ = ex(γ * Δ), ex(TR.C2 * Δ); Zs = [Z for Z in eachindex(eg_) if eg_[Z] < 1 || ec_[Z] < 1]; dEV = ρ[2][:EV] .- ρ[3][:EV]
        @printf("STAGECONS jstep %d item7 excap<1 in %d zones (dt=gammaDelta) / %d zones (dt=C2Delta); zones %s; excap(gammaDelta) %s; excap(C2Delta) %s; E_r row rho(1,gammaDelta)-rho(1,C2Delta) %s\n", J, count(<(1), eg_), count(<(1), ec_),
                join(Zs, ","), join((@sprintf("%.3g", eg_[Z]) for Z in Zs), ","), join((@sprintf("%.3g", ec_[Z]) for Z in Zs), ","), join((@sprintf("%.2e", dEV[a+Z-1]) for Z in Zs), ","))
    end
    ρ[2]
end

"Chemistry rate per unit species (S/x, 1/s) of the coupled species in the residual's last evaluation (W.cpS1,2 / W.xs1,2)."
chemrates(W) = Dict(:x1 => [W.xs1[k] != 0 ? W.cpS1[k] / W.xs1[k] : NaN for k in eachindex(W.xs1)], :x2 => [W.xs2[k] != 0 ? W.cpS2[k] / W.xs2[k] : NaN for k in eachindex(W.xs2)])
const LASTCHEM = Ref{Any}(nothing)
const FIELDS = (:r, :m, :P, :me, :mer, :EV, :QF, :x1, :x2)
"The rates of ρi per conserved combination: row E turned into :me (gas) and :mer (gas + radiation) according to the zone's row type; :r only if `withr`."
function fieldset(ρi, ty, Cn, withr)
    F = Dict{Symbol,Vector{Float64}}(nm => ρi[nm] for nm in (:m, :P, :EV, :QF, :x1, :x2) if haskey(ρi, nm)); withr && (F[:r] = ρi[:r])
    if haskey(ρi, :EV)
        F[:me] = [ty[k] === :gas ? ρi[:E][k] : (ρi[:E][k] * Cn[:mer][k] - ρi[:EV][k] * Cn[:EV][k]) / Cn[:me][k] for k in eachindex(ρi[:E])]
        F[:mer] = [ty[k] === :gas ? (ρi[:E][k] * Cn[:me][k] + ρi[:EV][k] * Cn[:EV][k]) / Cn[:mer][k] : ρi[:E][k] for k in eachindex(ρi[:E])]
    end
    F
end
function item2(J, Cn, a, Δ; mesh = :none, disp = nothing)
    D = TR.DSAL[]; tag = mesh === :fsal ? "item2b" : mesh === :tr ? "item2c" : "item2"; b = Int(TD.g(:ngre)); dt = TR.γ * Δ
    if D === nothing
        @printf("STAGECONS jstep %d %s skipped: DSAL === nothing (first step after a restart: no FSAL derivative, the step runs without an error estimate)\n", J, tag); return nothing
    end
    Cnw = nothing
    if mesh !== :none      # items 2b/2c: as item 2, but r and x_me displaced by dt × the FSAL velocities (2b: the mesh velocity of the previous step's BDF2 stage) or by the TR stage's own displacement (2c, called at the TR hook), everything else new = old; the time-difference part of every row is subtracted
        for (n_, key, o_, i_) in ((:rn, :r, :ro, 1), (:xmen, :xme, :xmeo, 2)); x = vf(n_); xo = vf(o_); for k in a:b+1; x[k] = mesh === :tr ? disp[i_][k] : xo[k] + dt * D[key][k]; end; end
        TF.var(:rhs) .= 0.0; TD.r(:update); Cnw = TR.combos_new()
        @printf("STAGECONS jstep %d %s: mesh displaced by %s; face velocity (r_displaced - r_old)/dt at zone-32 inner face %.6e, FSAL D[:r] there %.6e, u_old there %.6e\n", J, tag,
                mesh === :tr ? "the TR stage's (r_gamma - r_n), (xme_gamma - xme_n)" : "dt*D[:r], dt*D[:xme]", (vf(:rn)[a+31] - vf(:ro)[a+31]) / dt, D[:r][a+31], vf(:uo)[a+31])
    end
    Ri, ci, Wi = evalR(1.0, dt; etotal = false); ρi = rates(Ri, ci, Cn; Cnew = Cnw, dt); ty = etypes(ci); TX.COUPLED[] && (LASTCHEM[] = chemrates(Wi))
    if mesh === :fsal       # item 1b: θ-independence at the displaced state (dmdt, u_rel ≠ 0 there, so the old-level face values enter): ρ(½) − ρ(1), same dt, same state
        Rh, ch, _ = evalR(0.5, dt; etotal = false); ρh = rates(Rh, ch, Cn); ρ1 = rates(Ri, ci, Cn)
        for nm in (:m, :P, :E, :EV, :QF, :x1, :x2)
            haskey(ρh, nm) || continue
            d = ρh[nm] .- ρ1[nm]; mx, z = maxabs(d, a)
            @printf("STAGECONS jstep %d item1b %-4s max|R(½)-R(1)| %.3e zone %d | zones %d..%d R½ %s R1 %s diff %s\n", J, nm, mx, z, ZS[1], ZS[end], zl(ρh[nm], a), zl(ρ1[nm], a), zl(d, a))
        end
    end
    Fs = fieldset(ρi, ty, Cn, mesh === :none); pairs = [(nm, nm, Fs[nm]) for nm in FIELDS if haskey(Fs, nm) && haskey(D, nm)]
    for (nm, key, f) in pairs
        d = [ (k <= length(D[key]) && Cn[key][k] != 0) ? f[k] - D[key][k] / Cn[key][k] : NaN for k in eachindex(f)]
        mx, z = maxabs(d, a); mesh === :fsal && gate!(:item2b, mx, nm, z)
        @printf("STAGECONS jstep %d %s %-4s max|f-D|/C %.3e zone %d | zones %d..%d f/C %s D/C %s (f-D)/C %s\n", J, tag, nm, mx, z, ZS[1], ZS[end], zl(f, a),
                zl([k <= length(D[key]) && Cn[key][k] != 0 ? D[key][k] / Cn[key][k] : NaN for k in eachindex(f)], a), zl(d, a))
    end
    Fs
end

"Item 5 (at the TR hook): the TR stage's mean rate G = (C_γ − C_n)/(γΔ C_n), the flux part F_n of the row at y_n (with the TR stage's mesh motion, item 2c), F_γ at the converged TR state (θ = 1, old = y_n, flux part), and the FSAL D. ε_pred = 2(G − D/C) = (z2 − z1)/(C Δ) of DTERRDIAG."
function item5(J, Cn, a, Δ, Cg, Fn, Fg, Sg)
    D = TR.DSAL[]; (D === nothing || Fn === nothing) && return
    for nm in FIELDS
        (nm === :r || !haskey(Fn, nm) || !haskey(Fg, nm) || !haskey(D, nm)) && continue
        G = [Cn[nm][k] != 0 ? (Cg[nm][k] - Cn[nm][k]) / (TR.γ * Δ * Cn[nm][k]) : NaN for k in eachindex(Fn[nm])]
        Dn = [Cn[nm][k] != 0 ? D[nm][k] / Cn[nm][k] : NaN for k in eachindex(Fn[nm])]
        hs = 0.5 .* (Fn[nm] .+ Fg[nm]); ep = 2 .* (G .- Dn); dl = G .- hs
        mx, z = maxabs(ep, a)
        if haskey(Sg, nm) && LASTCHEM[] !== nothing      # species: chemistry part S/x and the rest (advection) of the flux part, at y_n (TR mesh) and at y_γ
            Sn = LASTCHEM[][nm]; Sγ = Sg[nm]
            @printf("STAGECONS jstep %d item5chem %-3s zones %d..%d S_n/x %s S_gamma/x %s A_n=F_n-S_n/x %s A_gamma=F_gamma-S_gamma/x %s\n", J, nm, ZS[1], ZS[end], zl(Sn, a), zl(Sγ, a), zl(Fn[nm] .- Sn, a), zl(Fg[nm] .- Sγ, a))
        end
        @printf("STAGECONS jstep %d item5 %-4s max|2(G-D/C)| %.3e zone %d | zones %d..%d G %s F_n(TRmesh) %s F_gamma %s G-(F_n+F_gamma)/2 %s D/C %s eps_pred=2(G-D/C) %s\n", J, nm, mx, z, ZS[1], ZS[end], zl(G, a), zl(Fn[nm], a), zl(Fg[nm], a), zl(dl, a), zl(Dn, a), zl(ep, a))
    end
end
"Flux-part rates at the current (converged TR) state: θ = 1, dt γΔ, old = yₙ, row energy of the etotal = false form."
function fgamma(Cn, Δ)
    Cg = TR.combos_new(); R, c, W = evalR(1.0, TR.γ * Δ; etotal = false); ρ = rates(R, c, Cn; Cnew = Cg, dt = TR.γ * Δ)
    Cg, fieldset(ρ, etypes(c), Cn, false), TX.COUPLED[] ? chemrates(W) : Dict{Symbol,Vector{Float64}}()
end

function item3(J, tag, Cn, a)
    R, c, _ = evalR(TD.g(:thet), TD.g(:dtime)); ρ = rates(R, c, Cn)
    let dt = TD.g(:dtime), D = TR.DSAL[], rn = vf(:rn), ro = vf(:ro), xn = vf(:xmen), xo = vf(:xmeo), un = vf(:un), uo = vf(:uo), ks = [a + Z - 1 for Z in 31:34]
        f(v) = join((@sprintf("%.4e", x) for x in v), " ")
        @printf("STAGECONS jstep %d item3 %s mesh at inner faces of zones 31..34: v=(rn-ro)/dt [%s] u_new [%s] u_old [%s] (xmeo-xmen)/dt [%s] FSAL D[:r] [%s]\n", J, tag, f((rn[k] - ro[k]) / dt for k in ks), f(un[k] for k in ks),
                f(uo[k] for k in ks), f((xo[k] - xn[k]) / dt for k in ks), D === nothing ? "n/a" : f(D[:r][k] for k in ks))
    end
    @printf("STAGECONS jstep %d item3 stage %s theta %.2f dt %.4e try %d: Newton iterations %d final dmax %.3e\n", J, tag, TD.g(:thet), TD.g(:dtime), Int(TD.g(:jtry)), Int(TD.g(:iter)), TD.g(:dmax))
    for s in ROWS
        haskey(ρ, s) || continue
        mx, z = maxabs(ρ[s], a); @printf("STAGECONS jstep %d item3 %s %-4s max|rho| %.3e zone %d | zones %d..%d rho %s\n", J, tag, s, mx, z, ZS[1], ZS[end], zl(ρ[s], a))
    end
end

function rel(old, new, rng)
    best = 0.0; kb = first(rng)
    for k in rng; x, y = old[k], new[k]; s = max(abs(x), abs(y)); r = s > 0 ? abs(x - y) / s : 0.0; (isfinite(r) && r > best) && (best = r; kb = k); end
    (best, kb)
end
"Item 4: the old-level arrays of the residual's context `c` (and the closure) against the new-level ones the residual computes (workspace W) at the start state."
function item4(J, c, W, Cn, a)
    b = c.ngre; r1 = a:b; r2 = a:b+1; P = c.p; lrad = c.lrad
    L = Tuple{String,Any,Any,Any}[("pgo/pgn", c.pgo, W.pgn, r1), ("ego/egn", c.ego, W.egn, r1), ("aso/asn", c.aso, W.asn, r1), ("dvolo/dvoln", c.dvolo, W.dvoln, r1),
        ("do/dn", c.do_, W.dn, r1), ("ro/rn", c.ro, W.rn, r2), ("xmeo/xmen", c.xmeo, W.xmen, r2), ("uo/un", c.uo, W.un, r2),
        ("xmo/(xmtot-xmen)", c.xmo, P.xmtot .- W.xmen, r2), ("oldslope d/dsn", W.dso, W.dsn, r2), ("oldslope u/usn", W.uso, W.usn, r2)]
    if lrad
        append!(L, [("chifo/chifn", c.chifo, W.chifn, r1), ("xkeo/xken", c.xkeo, W.xken, r1), ("xkpo/xkpn", c.xkpo, W.xkpn, r1), ("plfo/plfn", c.plfo, W.plfn, r1),
            ("egro/egrn", c.egro, W.egrn, r1), ("egro/(ego+ero/do) [old only]", c.egro, c.ego .+ c.ero ./ c.do_, r1), ("ero/ern", c.ero, W.ern, r1), ("fro/frn", c.fro, W.frn, r2),
            ("avchio/harm(chifn)", vf(:avchio), [k > a ? 2 / (1 / W.chifn[k-1] + 1 / W.chifn[k]) : 0.0 for k in eachindex(W.chifn)], a+1:b),
            ("avchio/avchin(COMMON)", vf(:avchio), vf(:avchin), a+1:b), ("oldslope egr/egsn", W.egrso, W.egsn, r2), ("oldslope er/ersn", W.erso, W.ersn, r2), ("oldslope fr/frsn", W.frso, W.frsn, r2)])
        TX.GASROW[] && push!(L, ("oldslope eg/egsng", W.egso, W.egsng, r2))
    else
        push!(L, ("oldslope eg/egsn", W.egso, W.egsn, r2))
    end
    if TX.COUPLED[]
        for i in 1:2
            push!(L, ("XSO$i/xs$i", TX.XSO[][i], i == 1 ? W.xs1 : W.xs2, r1), ("oldslope x$i/xsn$i", i == 1 ? W.xso1 : W.xso2, i == 1 ? W.xsn1 : W.xsn2, r2),
                     ("XSRCO$i/cpS$i", TX.XSRCO[][i], i == 1 ? W.cpS1 : W.cpS2, r1))
        end
        push!(L, ("XLAMO/cpL", TX.XLAMO[], W.cpL, r1))
    end
    TX.MATERIAL[] !== nothing && TX.has_source(TX.MATERIAL[]) && push!(L, ("SRCO/snuc", TX.SRCO[], W.snuc, r1))
    gh = [a-1, b+1, b+2]       # ghost cells: the old-level slope at the outer faces reads the old ghost values
    append!(L, [("ghost do/dn", c.do_, W.dn, gh), ("ghost uo/un", c.uo, W.un, gh), ("ghost ego/egn", c.ego, W.egn, gh), ("ghost ro/rn", c.ro, W.rn, gh)])
    lrad && append!(L, [("ghost egro/egrn", c.egro, W.egrn, gh), ("ghost ero/ern", c.ero, W.ern, gh), ("ghost fro/frn", c.fro, W.frn, gh)])
    for (nm, o, n, rng) in L
        m, k = rel(o, n, rng)
        @printf("STAGECONS jstep %d item4 %-30s max rel diff %.3e zone %d | zones %d,%d,%d old %.6e %.6e %.6e new %.6e %.6e %.6e\n", J, nm, m, k - a + 1, ZV[1], ZV[2], ZV[3],
                o[a+ZV[1]-1], o[a+ZV[2]-1], o[a+ZV[3]-1], n[a+ZV[1]-1], n[a+ZV[2]-1], n[a+ZV[3]-1])
    end
    let k = b       # slope debug at the last zone: old/new gas-energy slope and the values they are made from
        gq = TX.GASROW[] ? W.egsng : W.egsn
        @printf("STAGECONS jstep %d item4 slope@k=b(zone %d): ego[k-1,k,k+1,k+2] %.9e %.9e %.9e %.9e | egn %.9e %.9e %.9e %.9e | a_old %.3e b_old %.3e a_new %.3e b_new %.3e | oldslope %.6e newslope %.6e\n", J, k - a + 1,
                c.ego[k-1], c.ego[k], c.ego[k+1], c.ego[k+2], W.egn[k-1], W.egn[k], W.egn[k+1], W.egn[k+2], c.ego[k] - c.ego[k-1], c.ego[k+1] - c.ego[k], W.egn[k] - W.egn[k-1], W.egn[k+1] - W.egn[k], W.egso[k], gq[k])
    end
    # closure arrays: the start state (new = old, after update!) against the end of the previous step
    E = ENDCL[]; now = closure()
    if E !== nothing
        for (nm, o, n) in (("fedd", E.fedd, now.fedd), ("hedd", E.hedd, now.hedd), ("hxedd", E.hxedd, now.hxedd), ("fxedd", E.fxedd, now.fxedd), ("tiew", E.tiew, now.tiew))
            (length(o) == length(n) && length(o) > b + 1) || continue
            m, k = rel(o, n, a:min(b, length(o)))
            @printf("STAGECONS jstep %d item4 closure %-6s end-of-last-step/start-state: max rel diff %.3e zone %d | zones %d,%d,%d end %.6e %.6e %.6e start %.6e %.6e %.6e\n", J, nm, m, k - a + 1, ZV[1], ZV[2], ZV[3],
                    o[a+ZV[1]-1], o[a+ZV[2]-1], o[a+ZV[3]-1], n[a+ZV[1]-1], n[a+ZV[2]-1], n[a+ZV[3]-1])
        end
    end
end

"Item 6: the old state of this step against the state at the end of the last converged BDF2 stage (anything done to the state between the steps shows here): max relative difference with zone, values in zones 31-33."
function item6(J, a, b)
    E = ENDST[]; E === nothing && return
    for (n_, o_) in PRIM
        m, k = rel(E.prim[n_], vf(o_), a:b+1)
        @printf("STAGECONS jstep %d item6 %-5s end-of-last-step/old: max rel diff %.3e zone %d | zones 31,32,33 end %.9e %.9e %.9e old %.9e %.9e %.9e\n", J, string(n_), m, k - a + 1, [E.prim[n_][a+Z-1] for Z in ZV]..., [vf(o_)[a+Z-1] for Z in ZV]...)
    end
    for i in 1:2
        m, k = rel(E.xs[i], TX.XSO[][i], a:b)
        @printf("STAGECONS jstep %d item6 x%d    end-of-last-step/old: max rel diff %.3e zone %d | zones 31,32,33 end %.9e %.9e %.9e old %.9e %.9e %.9e\n", J, i, m, k - a + 1, [E.xs[i][a+Z-1] for Z in ZV]..., [TX.XSO[][i][a+Z-1] for Z in ZV]...)
    end
end
"new := old for the seven primary unknowns (a-1..b+2) and the species, then update! with a zero correction (ghost cells, EOS, opacities, closure)."
function new_eq_old!(a, b)
    for (n_, o_) in (:rn => :ro, :xmen => :xmeo, :dn => :do, :un => :uo, :tn => :to, :ern => :ero, :frn => :fro); copyto!(view(vf(n_), a-1:b+2), view(vf(o_), a-1:b+2)); end
    for i in 1:2; copyto!(TX.XSN[][i], TX.XSO[][i]); end
    TF.var(:rhs) .= 0.0; TD.r(:update)
end
"STAGEDUMP=zlo,zhi (with STAGECONS): the old state of zones zlo..zhi at every printed step start (T, rho, E_r, x_H+, x_H2, chi_f, kappa_E, f); STAGEONLY=1: nothing else is evaluated."
function stagedump(J, a)
    zl_, zh_ = parse.(Int, split(ENV["STAGEDUMP"], ","))
    @printf("STAGEDUMP jstep %d dt %.4e: Z  T  rho  E_r  x_H+  x_H2  chif  kappaE  fedd (old level)\n", J, TD.g(:dtime))
    for Z in zl_:zh_
        k = a + Z - 1
        @printf("STAGEDUMP jstep %d Z %d T %.6e rho %.6e Er %.6e xHp %.6e xH2 %.6e chif %.6e kE %.6e f %.6e\n", J, Z, vf(:to)[k], vf(:do)[k], vf(:ero)[k], TX.XSO[][1][k], TX.XSO[][2][k], vf(:chifo)[k], vf(:xkeo)[k], vf(:fedd)[k])
    end
end
const ACTIVE = Ref(false)
function hook(what, Δ, Cn)
    what === :bdf && (ENDCL[] = closure(); ENDST[] = endstate())
    what === :start && (ACTIVE[] = TR.ESTATS[:steps] % EVERY[] == 0)
    ACTIVE[] || return nothing
    (what !== :start && get(ENV, "STAGEONLY", "0") == "1") && return nothing
    J = Int(TD.g(:jstep)); a = Int(TD.g(:ngrs)); b = Int(TD.g(:ngre))
    if what === :start && haskey(ENV, "STAGEDUMP")
        stagedump(Int(TD.g(:jstep)), a); get(ENV, "STAGEONLY", "0") == "1" && return nothing
    end
    S = NK.Snap(); NK.take!(S); rhs0 = copy(TF.var(:rhs)); θ0 = TD.g(:thet); dt0 = TD.g(:dtime)
    try
        if what === :start
            ds = PRIM; item6(J, a, b)
            new_eq_old!(a, b)              # ghost cells and derived arrays of the new = old state (no predict!)
            mx = [maximum(abs(vf(n_)[k] - vf(o_)[k]) for k in a:b+1) for (n_, o_) in ds]
            mx = vcat(mx, [maximum(abs(TX.XSN[][i][k] - TX.XSO[][i][k]) for k in a:b) for i in 1:2])
            @printf("STAGECONS jstep %d start Delta %.4e (try %d): new := old, max|new-old| over zones a..b+1 for r m rho u T E_r F_r x1 x2: %s\n", J, Δ, Int(TD.g(:jtry)), join((@sprintf("%.1e", x) for x in mx), " "))
            ce = TX.context(); @printf("STAGECONS jstep %d Etype zones %d..%d: %s; FLUXTIE weights TIEW: %s\n", J, ZS[1], ZS[end], join((etypes(ce)[a+Z-1] for Z in ZS), " "), length(TD.TIEW) > b ? zl(TD.TIEW, a) : "n/a")
            empty!(GATE); item1(J, Δ, Cn, a); item2(J, Cn, a, Δ); item2(J, Cn, a, Δ; mesh = :fsal)
            gt(k) = get(GATE, k, (NaN, :none, 0))
            @printf("STAGECONS jstep %d gate: item1 max %.3e (%s zone %d) | item1dt max %.3e (%s zone %d) | item2b max %.3e (%s zone %d)\n", J, gt(:item1)..., gt(:item1dt)..., gt(:item2b)...)
        else
            item3(J, what === :tr ? "TR" : "BDF2", Cn, a)
            if what === :tr                # item 2c: the start state with the TR stage's own mesh displacement
                rg = copy(vf(:rn)); xg = copy(vf(:xmen)); Cg, Fg, Sg = fgamma(Cn, Δ); TD.s!(:thet, θ0); TD.s!(:dtime, dt0)
                new_eq_old!(a, b); Fn = item2(J, Cn, a, Δ; mesh = :tr, disp = (rg, xg)); item5(J, Cn, a, Δ, Cg, Fn, Fg, Sg)
            end
        end
    finally
        NK.restore!(S); TF.var(:rhs) .= rhs0; TD.s!(:thet, θ0); TD.s!(:dtime, dt0)
    end
    flush(stdout); nothing
end

"""
Item 8 (STAGEPEQ=1, with STAGECONS): at the step start, inside `peq!` after new → old and before `predict!` moves the new level, the derived arrays are recomputed at the unmoved state
with the step's held inputs (TitanCore._ghost_and_derived!(:peq), EOS, opacity, Eddington closure) and compared with what was just copied to the old level (pgo, ego, aso, chifo, xkeo,
xkpo, plfo, egro, avchio, xneo) and with the stored closure (fedd, hedd, hxedd, fxedd, tiew). Everything is restored afterwards.
"""
function peqcheck()
    TR.ESTATS[:steps] % EVERY[] == 0 || return nothing
    a = Int(TD.g(:ngrs)); b = Int(TD.g(:ngre)); J = Int(TD.g(:jstep))
    S = NK.Snap(); NK.take!(S); rhs0 = copy(TF.var(:rhs)); cl0 = closure()
    try
        L = ((:pgo, :pgn), (:ego, :egn), (:aso, :asn), (:chifo, :chifn), (:xkeo, :xken), (:xkpo, :xkpn), (:plfo, :plfn), (:egro, :egrn), (:avchio, :avchin), (:xneo, :xnen))
        old = Dict(o => copy(vf(o)) for (o, _) in L)
        for i in 1:2; copyto!(TX.XSN[][i], TX.XSO[][i]); end      # (cpl_newstep! has just replaced XSN by the local-chemistry guess of the new step; the old state's species are XSO)
        TF.var(:rhs) .= 0.0; TC._ghost_and_derived!(:peq, () -> TD.EOS[](), () -> TD.OPAC[]())
        cl1 = closure()
        for (o, n) in L
            m, k = rel(old[o], vf(n), a:b)
            @printf("STAGECONS jstep %d item8 %-8s stored/recomputed: max rel diff %.3e zone %d | zones 31,32,33 stored %.6e %.6e %.6e recomputed %.6e %.6e %.6e\n", J, string(o), m, k - a + 1, [old[o][a+Z-1] for Z in ZV]..., [vf(n)[a+Z-1] for Z in ZV]...)
        end
        for nm in (:fedd, :hedd, :hxedd, :fxedd, :tiew)
            o_, n_ = getfield(cl0, nm), getfield(cl1, nm); (length(o_) > b + 1 && length(n_) == length(o_)) || continue
            m, k = rel(o_, n_, a:b)
            @printf("STAGECONS jstep %d item8 %-8s stored/recomputed: max rel diff %.3e zone %d | zones 31,32,33 stored %.6e %.6e %.6e recomputed %.6e %.6e %.6e\n", J, string(nm), m, k - a + 1, [o_[a+Z-1] for Z in ZV]..., [n_[a+Z-1] for Z in ZV]...)
        end
    finally
        NK.restore!(S); TF.var(:rhs) .= rhs0
    end
    nothing
end

"STAGECONS=n: install the hook (every n-th accepted step)."
function configure!(n::Int)
    EVERY[] = n; TR.STAGECONS_HOOK[] = hook
    if get(ENV, "STAGEPEQ", "0") == "1"; TC.PEQ_HOOK[] === nothing || error("STAGEPEQ: TitanCore.PEQ_HOOK already in use"); TC.PEQ_HOOK[] = peqcheck; end
    nothing
end

end # module
