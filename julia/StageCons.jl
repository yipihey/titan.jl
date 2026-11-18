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
θ and dt; (2) −R(new = old) against the FSAL derivative DSAL of the previous step, (f − D)/C (energy from the row with etotal = false); (2b) the same with r and x_me displaced by Δt × their FSAL velocities (the rows contain the mesh velocity (rn − ro)/dt and the mass flux dmdt = −(xmen − xmeo)/dt, which are 0 at new = old) and the time-difference part (C_new − C_old)/dt of each row removed; (3) the true
residual at each converged stage (TR: θ ½, dt γΔ, old = yₙ; BDF2: θ 1, dt C2Δ, old = y*), with Newton iterations and final dmax; (4) every old-level array the residual reads against its
new-level counterpart at the start state (max relative difference with zone, values in zones 31–33), and the closure arrays against their values at the end of the previous step.
Zones are numbered from the inner boundary (zone Z = array index ngrs + Z − 1, as in DTERRDIAG). Skipped on a step without FSAL (first step after a restart: DSAL === nothing).
"""
module StageCons

using Printf
import ..TitanDriver, ..TitanResidFast, ..NewtonKrylov, ..Trbdf2
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const NK = NewtonKrylov; const TR = Trbdf2
vf(s) = TF.var(s)::Vector{Float64}
const EVERY = Ref(0); const ZS = 28:34; const ZV = 31:33
const ENDCL = Ref{Any}(nothing)       # closure arrays at the end of the last converged BDF2 stage

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
    γ = TR.γ; cases = ((0.5, γ * Δ), (1.0, γ * Δ), (1.0, TR.C2 * Δ), (0.5, 10γ * Δ)); ρ = Dict{Int,Any}()
    for (i, (θ, dt)) in enumerate(cases)
        R, c, W = evalR(θ, dt); ρ[i] = rates(R, c, Cn)
        i == 1 && item4(J, c, W, Cn, a)
    end
    for s in ROWS
        haskey(ρ[1], s) || continue
        h, o, o2, h10 = ρ[1][s], ρ[2][s], ρ[3][s], ρ[4][s]
        m1, z1 = maxabs(h .- o, a); m2, z2 = maxabs(h .- h10, a); m3, z3 = maxabs(o .- o2, a)
        @printf("STAGECONS jstep %d item1 %-4s max|R(½,γΔ)-R(1,γΔ)| %.3e zone %d | max|R(½,γΔ)-R(½,10γΔ)| %.3e zone %d | max|R(1,γΔ)-R(1,C2Δ)| %.3e zone %d | zones %d..%d R½ %s R1 %s diff %s\n",
                J, s, m1, z1, m2, z2, m3, z3, ZS[1], ZS[end], zl(h, a), zl(o, a), zl(h .- o, a))
    end
    ρ[2]
end

function item2(J, Cn, a, Δ; mesh = :none, disp = nothing)
    D = TR.DSAL[]; tag = mesh === :fsal ? "item2b" : mesh === :tr ? "item2c" : "item2"; b = Int(TD.g(:ngre)); dt = TR.γ * Δ
    if D === nothing
        @printf("STAGECONS jstep %d %s skipped: DSAL === nothing (first step after a restart: no FSAL derivative, the step runs without an error estimate)\n", J, tag); return
    end
    Cnw = nothing
    if mesh !== :none      # items 2b/2c: as item 2, but r and x_me displaced by dt × the FSAL velocities (2b: the mesh velocity of the previous step's BDF2 stage) or by the TR stage's own displacement (2c, called at the TR hook), everything else new = old; the time-difference part of every row is subtracted
        for (n_, key, o_, i_) in ((:rn, :r, :ro, 1), (:xmen, :xme, :xmeo, 2)); x = vf(n_); xo = vf(o_); for k in a:b+1; x[k] = mesh === :tr ? disp[i_][k] : xo[k] + dt * D[key][k]; end; end
        TF.var(:rhs) .= 0.0; TD.r(:update); Cnw = TR.combos_new()
        @printf("STAGECONS jstep %d %s: mesh displaced by %s; face velocity (r_displaced - r_old)/dt at zone-32 inner face %.6e, FSAL D[:r] there %.6e, u_old there %.6e\n", J, tag,
                mesh === :tr ? "the TR stage's (r_gamma - r_n), (xme_gamma - xme_n)" : "dt*D[:r], dt*D[:xme]", (vf(:rn)[a+31] - vf(:ro)[a+31]) / dt, D[:r][a+31], vf(:uo)[a+31])
    end
    Ri, ci, _ = evalR(1.0, dt; etotal = false); ρi = rates(Ri, ci, Cn; Cnew = Cnw, dt); ty = etypes(ci)
    pairs = Tuple{Symbol,Symbol,Vector{Float64}}[]
    for nm in (mesh !== :none ? (:m, :P, :EV, :QF, :x1, :x2) : (:r, :m, :P, :EV, :QF, :x1, :x2)); haskey(ρi, nm) && haskey(D, nm) && push!(pairs, (nm, nm, ρi[nm])); end
    if haskey(D, :me) && haskey(D, :mer) && haskey(ρi, :EV)         # E row: gas (m e_g) or mer (m e_g + E_r ΔV) form
        fme = [ty[k] === :gas ? ρi[:E][k] : (ρi[:E][k] * Cn[:mer][k] - ρi[:EV][k] * Cn[:EV][k]) / Cn[:me][k] for k in eachindex(ρi[:E])]
        fmer = [ty[k] === :gas ? (ρi[:E][k] * Cn[:me][k] + ρi[:EV][k] * Cn[:EV][k]) / Cn[:mer][k] : ρi[:E][k] for k in eachindex(ρi[:E])]
        push!(pairs, (:me, :me, fme), (:mer, :mer, fmer))
    end
    for (nm, key, f) in pairs
        d = [ (k <= length(D[key]) && Cn[key][k] != 0) ? f[k] - D[key][k] / Cn[key][k] : NaN for k in eachindex(f)]
        mx, z = maxabs(d, a)
        @printf("STAGECONS jstep %d %s %-4s max|f-D|/C %.3e zone %d | zones %d..%d f/C %s D/C %s (f-D)/C %s\n", J, tag, nm, mx, z, ZS[1], ZS[end], zl(f, a),
                zl([k <= length(D[key]) && Cn[key][k] != 0 ? D[key][k] / Cn[key][k] : NaN for k in eachindex(f)], a), zl(d, a))
    end
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
    L = Tuple{String,Any,Any,UnitRange{Int}}[("pgo/pgn", c.pgo, W.pgn, r1), ("ego/egn", c.ego, W.egn, r1), ("aso/asn", c.aso, W.asn, r1), ("dvolo/dvoln", c.dvolo, W.dvoln, r1),
        ("do/dn", c.do_, W.dn, r1), ("ro/rn", c.ro, W.rn, r2), ("xmeo/xmen", c.xmeo, W.xmen, r2), ("uo/un", c.uo, W.un, r2),
        ("xmo/(xmtot-xmen)", c.xmo, P.xmtot .- W.xmen, r2), ("dso/dsn", c.dso, W.dsn, r2), ("uso/usn", c.uso, W.usn, r2)]
    if lrad
        append!(L, [("chifo/chifn", c.chifo, W.chifn, r1), ("xkeo/xken", c.xkeo, W.xken, r1), ("xkpo/xkpn", c.xkpo, W.xkpn, r1), ("plfo/plfn", c.plfo, W.plfn, r1),
            ("egro/egrn", c.egro, W.egrn, r1), ("egro/(ego+ero/do) [old only]", c.egro, c.ego .+ c.ero ./ c.do_, r1), ("ero/ern", c.ero, W.ern, r1), ("fro/frn", c.fro, W.frn, r2),
            ("avchio/harm(chifn)", vf(:avchio), [k > a ? 2 / (1 / W.chifn[k-1] + 1 / W.chifn[k]) : 0.0 for k in eachindex(W.chifn)], a+1:b),
            ("avchio/avchin(COMMON)", vf(:avchio), vf(:avchin), a+1:b), ("egrso/egsn", c.egrso, W.egsn, r2), ("erso/ersn", c.erso, W.ersn, r2), ("frso/frsn", c.frso, W.frsn, r2)])
        TX.GASROW[] && push!(L, ("egso/egsng", c.egso, W.egsng, r2))
    end
    if TX.COUPLED[]
        for i in 1:2
            push!(L, ("XSO$i/xs$i", TX.XSO[][i], i == 1 ? W.xs1 : W.xs2, r1), ("XSSO$i/xsn$i", TX.XSSO[][i], i == 1 ? W.xsn1 : W.xsn2, r2),
                     ("XSRCO$i/cpS$i", TX.XSRCO[][i], i == 1 ? W.cpS1 : W.cpS2, r1))
        end
        push!(L, ("XLAMO/cpL", TX.XLAMO[], W.cpL, r1))
    end
    TX.MATERIAL[] !== nothing && TX.has_source(TX.MATERIAL[]) && push!(L, ("SRCO/snuc", TX.SRCO[], W.snuc, r1))
    for (nm, o, n, rng) in L
        m, k = rel(o, n, rng)
        @printf("STAGECONS jstep %d item4 %-30s max rel diff %.3e zone %d | zones %d,%d,%d old %.6e %.6e %.6e new %.6e %.6e %.6e\n", J, nm, m, k - a + 1, ZV[1], ZV[2], ZV[3],
                o[a+ZV[1]-1], o[a+ZV[2]-1], o[a+ZV[3]-1], n[a+ZV[1]-1], n[a+ZV[2]-1], n[a+ZV[3]-1])
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

"new := old for the seven primary unknowns (a-1..b+2) and the species, then update! with a zero correction (ghost cells, EOS, opacities, closure)."
function new_eq_old!(a, b)
    for (n_, o_) in (:rn => :ro, :xmen => :xmeo, :dn => :do, :un => :uo, :tn => :to, :ern => :ero, :frn => :fro); copyto!(view(vf(n_), a-1:b+2), view(vf(o_), a-1:b+2)); end
    for i in 1:2; copyto!(TX.XSN[][i], TX.XSO[][i]); end
    TF.var(:rhs) .= 0.0; TD.r(:update)
end
const ACTIVE = Ref(false)
function hook(what, Δ, Cn)
    what === :bdf && (ENDCL[] = closure())
    what === :start && (ACTIVE[] = TR.ESTATS[:steps] % EVERY[] == 0)
    ACTIVE[] || return nothing
    J = Int(TD.g(:jstep)); a = Int(TD.g(:ngrs)); b = Int(TD.g(:ngre))
    S = NK.Snap(); NK.take!(S); rhs0 = copy(TF.var(:rhs)); θ0 = TD.g(:thet); dt0 = TD.g(:dtime)
    try
        if what === :start
            ds = (:rn => :ro, :xmen => :xmeo, :dn => :do, :un => :uo, :tn => :to, :ern => :ero, :frn => :fro)
            new_eq_old!(a, b)              # ghost cells and derived arrays of the new = old state (no predict!)
            mx = [maximum(abs(vf(n_)[k] - vf(o_)[k]) for k in a:b+1) for (n_, o_) in ds]
            mx = vcat(mx, [maximum(abs(TX.XSN[][i][k] - TX.XSO[][i][k]) for k in a:b) for i in 1:2])
            @printf("STAGECONS jstep %d start Delta %.4e (try %d): new := old, max|new-old| over zones a..b+1 for r m rho u T E_r F_r x1 x2: %s\n", J, Δ, Int(TD.g(:jtry)), join((@sprintf("%.1e", x) for x in mx), " "))
            ce = TX.context(); @printf("STAGECONS jstep %d Etype zones %d..%d: %s; FLUXTIE weights TIEW: %s\n", J, ZS[1], ZS[end], join((etypes(ce)[a+Z-1] for Z in ZS), " "), length(TD.TIEW) > b ? zl(TD.TIEW, a) : "n/a")
            item1(J, Δ, Cn, a); item2(J, Cn, a, Δ); item2(J, Cn, a, Δ; mesh = :fsal)
        else
            item3(J, what === :tr ? "TR" : "BDF2", Cn, a)
            if what === :tr                # item 2c: the start state with the TR stage's own mesh displacement
                rg = copy(vf(:rn)); xg = copy(vf(:xmen)); new_eq_old!(a, b); item2(J, Cn, a, Δ; mesh = :tr, disp = (rg, xg))
            end
        end
    finally
        NK.restore!(S); TF.var(:rhs) .= rhs0; TD.s!(:thet, θ0); TD.s!(:dtime, dt0)
    end
    flush(stdout); nothing
end

"STAGECONS=n: install the hook (every n-th accepted step)."
function configure!(n::Int)
    EVERY[] = n; TR.STAGECONS_HOOK[] = hook; nothing
end

end # module
