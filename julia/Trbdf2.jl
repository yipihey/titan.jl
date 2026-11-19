"""
    Trbdf2

Second-order, L-stable time stepping for the Julia TITAN (Eulerian grid, lgrid = 2; Lagrangian grid, lgrid = 3): TR-BDF2 (Bank et al. 1985) built
from TITAN's own θ-scheme, so the residual is unchanged.

One step of size Δ from y_n:
1. **TR stage:** TITAN's θ-scheme with θ = ½ (implicit midpoint) over γΔ, γ = 2 − √2, giving y_γ.
2. **BDF2 stage:** y_{n+1} − [A y_γ − B y_n] = (1−γ)/(2−γ) Δ · f(y_{n+1}), A = 1/(γ(2−γ)), B = (1−γ)²/(γ(2−γ)).
   That is a backward-Euler step (θ = 1) of size (1−γ)/(2−γ)Δ from the "old" state y* = A y_γ − B y_n. TITAN's
   residual reads the old state only through conserved combinations (mass m = ρ dV, m e_g, m(e_g + E/ρ), E dV,
   mass-weighted momentum, the flux-row inertia F·2M/(ρ₋+ρ₊), the exterior mass coordinate), so y* is formed on those
   and written back into the old arrays (do, ego, egro, ero, uo, fro, xmeo). Terms with weight 1 − θ vanish at θ = 1.
Newton's initial guess of every solve is the grid extrapolated with the grid velocities over the stage (TitanCore.predict!,
the second half of peq: r = r_old + (dr/dt) dt, x_me likewise, ghost zones, derived quantities), also after a reset and for
the backward-Euler fallback. Without it the BDF stage and the fallback start from y_γ resp. y_n: with zones moving
u Δ ~ their width (the Δ-limiting infall envelope), Newton's first step is garbage (relative T changes of 1e2–1e6 in
distant zones through a block system with condition number 1e19), damped to 10 % per iteration, and the stage fails;
the θ-scheme never sees this because step.f starts every step with peq. peq itself is called with dtime = γΔ so that the
TR stage is extrapolated over γΔ, not Δ.
After the step the old arrays are restored to y_n, so TITAN's step-size control (timstp) and the energy ledger see the
whole step Δ. If either stage fails, the same Δ is retried with backward Euler (FALLBACK_BE); only if that fails too is
the step restarted from y_n with Δ × TFAC[]. Steps taken that way are first order (counted in STATS[:be_fallback]).

Lagrangian grid (lgrid = 3, Population III collapse, G2): the radii are combined linearly, ro* = A r_γ − B r_n, the old zone
volumes dvolo* are recomputed from ro*, and the old densities are ρ* = m*/dvolo* (m* = A m_γ − B m_n, which is m for a
Lagrangian zone), so that every time-derivative term of the residual, (ρ dV)_{n+1} − (ρ* dV*), reads the combined
conserved mass. Also saved/restored: ro, dvolo. With TitanCore.FREE_OUTER[] (PEXT) the outermost interface is dynamic and
its old velocity is combined with the mass of the outer half zone. REGULAR_CENTRE needs nothing (r = 0, u = 0 at the
centre), XMINT keeps x_me frozen (A − B = 1), EXTMASS and DT_DYN do not involve old-state arrays.
Limits: the energy ledger (totnrg) books the step as one backward-Euler step, so it is approximate.
Use: `TitanDriver.STEPFN[] = Trbdf2.step!`.
"""
module Trbdf2

import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..BlockSolver
using Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast
g(s) = TD.g(s); s!(s, x) = TD.s!(s, x); r(s) = TD.r(s)
vf(s) = TF.var(s)::Vector{Float64}

const γ = 2 - sqrt(2)
const A = 1 / (γ * (2 - γ)); const B = (1 - γ)^2 / (γ * (2 - γ)); const C2 = (1 - γ) / (2 - γ)
# A y_γ − B y_n written as y_γ + B (y_γ − y_n) (A = 1 + B): exact when the two states coincide (a frozen x_me of 1e47 would
# otherwise acquire round-off 1e31 and feed dm/dt ≠ 0 into the advective terms)
@inline cb(yγ, yn) = yγ + B * (yγ - yn)
const TFAC = Ref(0.05)
const STATS = Dict(:steps => 0, :fail_tr => 0, :fail_bdf => 0, :be_fallback => 0)
"If either stage fails, retry the same Δ with backward Euler (θ = 1) before cutting Δ (default true). The implicit
midpoint stage is not dissipative and fails at under-resolved shocks where backward Euler converges."
const FALLBACK_BE = Ref(true)
"After a fallback, the next STICKY[] steps go straight to backward Euler (saves the failing TR stage)."
const STICKY = Ref(10); const BE_LEFT = Ref(0)
const OLD = (:do, :ego, :egro, :ero, :fro, :uo, :xmeo, :ro, :dvolo)
const SAVE = Dict{Symbol,Vector{Float64}}()
# COUPLED: the old species (TitanResidFast.XSO) are saved, restored and combined with the other old arrays
const SAVEX = [Float64[], Float64[]]
save!() = (for s in OLD; copyto!(get!(() -> similar(vf(s)), SAVE, s), vf(s)); end;
           TX.COUPLED[] && for i in 1:2; length(SAVEX[i]) == length(TX.XSO[][i]) || (SAVEX[i] = similar(TX.XSO[][i])); copyto!(SAVEX[i], TX.XSO[][i]); end;
           TX.MG[] && mg_save!())
restore!() = (for s in OLD; copyto!(vf(s), SAVE[s]); end; TX.COUPLED[] && for i in 1:2; copyto!(TX.XSO[][i], SAVEX[i]); end;
              TX.MG[] && mg_restore!(); isempty(TX.KSTAR[]) || (TX.KSTAR[] = Float64[]))
# MGROUPS: the old group fields (TitanResidFast MGEO, MGFO) saved, restored and combined like E_r (on E ΔV) and F_r (on the flux-row inertia)
const SAVEMG = [Vector{Vector{Float64}}(), Vector{Vector{Float64}}()]
mg_save!() = (SAVEMG[1] = deepcopy(TX.MGEO[]); SAVEMG[2] = deepcopy(TX.MGFO[]); nothing)
mg_restore!() = (for g in 1:TX.MGN[]; copyto!(TX.MGEO[][g], SAVEMG[1][g]); copyto!(TX.MGFO[][g], SAVEMG[2][g]); end; nothing)
mgkeys(g) = (Symbol(:EV, g), Symbol(:QF, g))

"""
ERESTOL[] (provisional, ENV ERESTOL [L_sun], default 0 = off): residual-based stopping test of the Newton iteration, in addition to
the correction test (docs/POPIII_ENERGY.md §17). After the correction test has passed (and TITAN's closing full iteration), the energy
the solve's own rows leave unbalanced, |Σ_k (gas energy row + radiation energy row)_k + Σ_j u_j (gas momentum row)_j| × 4π, is evaluated
at the state reached (`eres`); while it exceeds ERESTOL, further full Newton iterations are taken, up to niter in all. The iteration
ends without meeting the tolerance when an iteration does not halve the measure (round-off of the rows; counted in ERES[:unmet], the
stage is accepted on the correction test as before). ERES counts the stages tested, the extra iterations and the stages left above the
tolerance; ERESMAX the largest measure accepted.
"""
const ERESTOL = Ref(0.0); const ERES = Dict(:checks => 0, :extra => 0, :unmet => 0, :niter => 0); const ERESMAX = Ref(0.0)
function eres()
    c = TX.context(); R = zeros(c.neqn, c.nb); ws = TX.getws(Float64, c.mgr)
    TX.residual!(R, ws, c, zeros(c.neqn * c.nb)); W = ws.a
    iu, it, ie = c.ix[4], c.ix[5], c.ix[6]; kb(k) = k - c.k0 + 1; s = 0.0
    for k in c.ngrs:c.ngre; s += R[it, kb(k)] + (c.lrad ? R[ie, kb(k)] : 0.0); end
    for j in c.ngrs+1:c.ngre; s += W.u[j] * R[iu, kb(j)]; end
    4π * abs(s) / 3.828e33
end
# NWTRACE=<file> (diagnostic, default off = bit-identical): one line per Newton iteration of every attempt (iteration, dmax, the base measure
# max(|dr|,|dm|,|drho|,|dT|,|dE|), the species' weighted measure with its equation and zone, and for each of the equations the largest
# raw correction with its zone), a header per attempt and a result line.
const NWDUMP = Ref{Any}(nothing)       # (diagnostic) state dump hook called at the start of the NWPROBE attempts
const NWTRACE = Ref{Any}(nothing); const NWATT = Ref(0); const NWPROBE = Ref(Float64[])
"STAGECONS=n (diagnostic, default off = bit-identical; StageCons.jl): `(what, Δ, Cn)` called with what = :start (after peq/save!, before the first attempt), :tr, :bdf (after a converged stage); it restores everything it touches."
const STAGECONS_HOOK = Ref{Any}(nothing)
# ---- stalldiag diagnostics (all default off; SDZONES=lo:hi turns the per-zone prints on) -----------------------------------------------
# SDZONES=lo:hi  zones (numbered from the inner boundary) printed: SDSTATE (state of each accepted step), SDEST (the weighted DTERR estimate of
#                every attempt, first filter f1 and second filter f2), SDFLAG (regime flags via SD_FLAGS[]), SDDEC (decide() flips, with SDDEC=1)
# SDFIXDT=<s>   fixed step size, never rejected for the error estimate (a Newton failure still cuts it): convergence/Richardson tests of the true error
const SDZ = let s = get(ENV, "SDZONES", ""); isempty(s) ? Int[] : collect(parse(Int, split(s, ":")[1]):parse(Int, split(s, ":")[2])) end
const SDON = !isempty(SDZ)
const SDFIXDT = parse(Float64, get(ENV, "SDFIXDT", "0"))
const SDCLO = SDON && get(ENV, "SDCLO", "0") == "1"   # per Newton iterate: closure arrays (h, f, tie weight, chi) of the zones in SDZONES
function sd_clo()
    a = Int(g(:ngrs)); fe, ch, er, tn = vf(:fedd), vf(:chifn), vf(:ern), vf(:tn); hE = TD.HEDD; tw = TD.TIEW
    io = IOBuffer()
    for zk in SDZ
        k = a + zk - 1
        ds = vf(:dsf); sfv = vf(:sf); bs = vf(:bsf)
        @printf(io, " z%d[S %.5e dsf- %.3e dsf+ %.3e bsf %.3e]", zk, sfv[k], ds[k-1], ds[k], bs[k])
        @printf(io, " z%d[h %.6f f %.6f w %.4f chi %.5e E %.5f T %.2f]", zk, k <= length(hE) ? hE[k] : NaN, fe[k], k <= length(tw) ? tw[k] : NaN, ch[k], er[k], tn[k])
    end
    println("SDCLO jstep ", Int(g(:jstep)), " thet ", g(:thet), " dt ", g(:dtime), " it ", Int(g(:iter)), " dmax ", g(:dmax), String(take!(io)))
end
const SDESTNK = get(ENV, "SDESTNK", "0") == "1"   # the estimate's linear solves also get the Krylov refinement (exact J v), as the Newton solves do; default: band J_b only
const SD_FLAGS = Ref{Any}(nothing)          # (popiii_m2c) k -> String: regime flags of zone k (table/network, coupled)
const SDACC = Dict{Tuple{Symbol,Int,Int},Float64}(); const SDTAG = Ref(:e0)
function nwtrace_it()
    io = NWTRACE[]; n = Int(g(:neqn)); kx = TF.var(:kx); dx = TF.var(:dx)
    base = max(abs(dx[Int(g(:ir))]), abs(dx[Int(g(:im))]), abs(dx[Int(g(:id))]), abs(dx[Int(g(:it))]), abs(dx[Int(g(:ie))])) + 1.0e-30
    cw = TX.CPLDM[]
    println(io, "NWT att ", NWATT[], " it ", Int(g(:iter)), " dmax ", g(:dmax), " base ", base, " species ", cw[1], " eq ", cw[2], " zone ", cw[3], " | ",
            join(("$(i):$(round(dx[i], sigdigits = 4))@$(kx[i])" for i in 1:n), " "))
    flush(io)
end
function newton!()
    NWTRACE[] === nothing && return newton0!()
    NWATT[] += 1
    println(NWTRACE[], "NWT head att ", NWATT[], " jstep ", Int(g(:jstep)), " thet ", g(:thet), " dtime ", g(:dtime), " niter ", Int(g(:niter)))
    ok = newton0!()
    println(NWTRACE[], "NWT end att ", NWATT[], " ok ", ok, " iter ", Int(g(:iter))); flush(NWTRACE[])
    ok
end
"""
Newton difficulty in the step-size control (all schemes with DTERR): `NQ[]` is the smallest, over the Newton solves since it was last reset to 1, of
f = clamp(ϑ_req/ϑ, 0.1, 1) with ϑ = (d_last/d_first)^(1/(n−1)) the geometric-mean contraction of the convergence measure `dmax` (of the FULL Newton step, also when the
globalisation cut it) over the n iterations of the solve, and ϑ_req = (conv/d_first)^(2/niter) the contraction that would converge in niter/2 iterations. A solve that
converges at least that fast has f = 1. `step_err!` bounds the next Δt by NQ·Δt after an accepted step and a failed attempt's retry by min(DTERRFAIL, max(0.1, NQ))·Δt.
"""
const NQ = Ref(1.0)
function newton_quality!(dm)
    n = length(dm); n < 2 && return
    d1, dl = dm[1], dm[end]; cv = g(:conv)
    d1 > cv || return
    ϑ = (dl / d1)^(1 / (n - 1)); ϑreq = (cv / d1)^(2 / g(:niter))
    f = isfinite(ϑ) && ϑ > 0 ? clamp(ϑreq / ϑ, 0.1, 1.0) : (isfinite(ϑ) ? 1.0 : 0.1)
    NQ[] = min(NQ[], f)
end
function newton0!()
    s!(:iter, 1); dm = Float64[]
    try
        while true
            TD.newton_linear!(); r(:maxdel); NWTRACE[] === nothing || nwtrace_it(); r(:update); SDCLO && sd_clo(); push!(dm, g(:dmax))
            g(:dmax) <= g(:conv) && break
            s!(:iter, g(:iter) + 1)
            g(:iter) > g(:niter) && (newton_quality!(dm); return false)
        end
        s!(:lboos, 1); TD.newton_linear!(); r(:maxdel); r(:update); s!(:lboos, 0)
        if ERESTOL[] > 0
            e = eres(); ERES[:checks] += 1
            while e > ERESTOL[]
                g(:iter) >= g(:niter) && (ERES[:niter] += 1; break)
                s!(:iter, g(:iter) + 1); s!(:lboos, 1); TD.newton_linear!(); r(:maxdel); r(:update); s!(:lboos, 0); ERES[:extra] += 1
                e1 = eres(); stall = !(e1 < 0.5e); e = e1
                stall && break
            end
            e > ERESTOL[] && (ERES[:unmet] += 1)
            ERESMAX[] = max(ERESMAX[], e)
        end
    catch e
        e isa InterruptException && rethrow()
        NWTRACE[] === nothing || println(NWTRACE[], "NWT exception ", first(sprint(showerror, e, catch_backtrace()), 700))
        s!(:lboos, 0); newton_quality!(dm); return false
    end
    newton_quality!(dm)
    all(isfinite, vf(:dn)) && all(isfinite, vf(:tn))
end

"Old arrays := y* = A·y_γ − B·y_n on conserved combinations (y_n in SAVE, y_γ in the new arrays). false if m* ≤ 0."
function bdf_old!()
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1
    dn, dvn, egn, un, xmen = vf(:dn), vf(:dvoln), vf(:egn), vf(:un), vf(:xmen)
    lag = Int(g(:lgrid)) == 3 || TitanCore.GRIDMOVE[]     # (fixed interfaces: cb(r, r) = r exactly)
    d0, dv0, eg0, u0, xm0 = SAVE[:do], SAVE[:dvolo], SAVE[:ego], SAVE[:uo], SAVE[:xmeo]
    do_, ego, uo, xmeo = vf(:do), vf(:ego), vf(:uo), vf(:xmeo)
    dvo = vf(:dvolo)                   # old zone volume read by the residual (y_n's; recomputed from r* when lagrangian)
    rng = a:b                          # zones read by the residual's time-derivative terms (ghosts keep y_n)
    mγ(k) = dn[k] * dvn[k]; mn(k) = d0[k] * dv0[k]
    ms = zeros(length(dn))
    @inbounds for k in rng
        ms[k] = cb(mγ(k), mn(k))
        ms[k] > 0 || return false
    end
    if lag                             # combined radii and the volumes they enclose
        rn_, ro, r0 = vf(:rn), vf(:ro), SAVE[:ro]; mu1 = Int(g(:mu)) + 1; xmup1 = g(:xmup1)
        @inbounds for k in a-1:b+2; ro[k] = cb(rn_[k], r0[k]); end
        @inbounds for k in a:b
            dvo[k] = (ro[k+1]^mu1 - ro[k]^mu1) / xmup1
            dvo[k] > 0 || return false
        end
    end
    @inbounds for k in rng
        ego[k] = cb(mγ(k) * egn[k], mn(k) * eg0[k]) / ms[k]
        do_[k] = ms[k] / dvo[k]
    end
    @inbounds for k in a-1:b+2; xmeo[k] = cb(xmen[k], xm0[k]); end
    if TX.COUPLED[]                    # species on ρ x ΔV (the residual reads x_old only as ρ_old x_old ΔV_old)
        for i in 1:2
            xn_, xo_ = TX.XSN[][i], TX.XSO[][i]
            @inbounds for k in rng; TX.cplz(k) && (xo_[k] = cb(mγ(k) * xn_[k], mn(k) * SAVEX[i][k]) / ms[k]); end
        end
    end
    if lrad
        egrn, ern, frn = vf(:egrn), vf(:ern), vf(:frn)
        egr0, er0, fr0 = SAVE[:egro], SAVE[:ero], SAVE[:fro]
        egro, ero, fro = vf(:egro), vf(:ero), vf(:fro)
        @inbounds for k in rng
            egro[k] = cb(mγ(k) * egrn[k], mn(k) * egr0[k]) / ms[k]
            ero[k] = cb(ern[k] * dvn[k], er0[k] * dv0[k]) / dvo[k]
        end
        @inbounds for k in a+1:b      # flux-row inertia Q = F·2M/(ρ₋ + ρ₊), M = m₋ + m₊
            Qγ = frn[k] * 2 * (mγ(k-1) + mγ(k)) / (dn[k-1] + dn[k])
            Qn = fr0[k] * 2 * (mn(k-1) + mn(k)) / (d0[k-1] + d0[k])
            fro[k] = cb(Qγ, Qn) * (do_[k-1] + do_[k]) / (2 * (ms[k-1] + ms[k]))
        end
        if TX.MG[]                     # MGROUPS: each group as E_r and F_r
            for g in 1:TX.MGN[]
                en_, fn_, eo_, fo_ = TX.MGEN[][g], TX.MGFN[][g], TX.MGEO[][g], TX.MGFO[][g]; e0_, f0_ = SAVEMG[1][g], SAVEMG[2][g]
                @inbounds for k in rng; eo_[k] = cb(en_[k] * dvn[k], e0_[k] * dv0[k]) / dvo[k]; end
                @inbounds for k in a+1:b
                    Qγ = fn_[k] * 2 * (mγ(k-1) + mγ(k)) / (dn[k-1] + dn[k]); Qn = f0_[k] * 2 * (mn(k-1) + mn(k)) / (d0[k-1] + d0[k])
                    fo_[k] = cb(Qγ, Qn) * (do_[k-1] + do_[k]) / (2 * (ms[k-1] + ms[k]))
                end
            end
        end
    end
    @inbounds for k in a+1:b          # momentum: u·(m₋ + m₊) (boundary velocities are fixed)
        uo[k] = cb(un[k] * (mγ(k-1) + mγ(k)), u0[k] * (mn(k-1) + mn(k))) / (ms[k-1] + ms[k])
    end
    if TitanCore.FREE_OUTER[]          # dynamic outer interface: momentum of the outer half zone
        k = b + 1
        uo[k] = cb(un[k] * mγ(k-1), u0[k] * mn(k-1)) / ms[k-1]
    end
    if TX.ECOMPAT[]                    # kinetic energy of y* minus the BDF2 combination of the stage kinetic energies (TX.KSTAR)
        D = zeros(length(dn))
        @inbounds for k in rng
            D[k] = 0.25 * (ms[k] * (uo[k]^2 + uo[k+1]^2) - cb(mγ(k) * (un[k]^2 + un[k+1]^2), mn(k) * (u0[k]^2 + u0[k+1]^2)))
        end
        TX.KSTAR[] = D
    end
    true
end

"Drop-in for TitanDriver.step! (TR-BDF2). With DTERR set, `step_err!` (error-controlled step size) instead."
function step!()
    ERRCTL[] && return step_err!()
    Int(g(:lgrid)) in (2, 3) || error("Trbdf2: lgrid = 2 or 3 only")
    # peq extrapolates the grid over dtime: that is the TR stage (γΔ), not the whole step (BE_LEFT > 0: backward Euler over Δ)
    Δ0 = g(:dtime); BE_LEFT[] > 0 || s!(:dtime, γ * Δ0)
    r(:peq)
    s!(:dtime, BE_LEFT[] > 0 ? g(:dtime) : g(:dtime) / γ); s!(:timen, g(:timeo) + g(:dtime))   # peq may have halved dtime
    θ0 = g(:thet); t0 = g(:timeo)
    save!()
    s!(:jtry, 1); s!(:jback, 0)
    while true
        Δ = g(:dtime); TD.FULL_DT[] = Δ
        if BE_LEFT[] > 0                  # sticky backward-Euler phase after a fallback
            BE_LEFT[] -= 1
            s!(:thet, 1.0); s!(:timen, t0 + Δ)
            ok = newton!()
            ok && TD.STAGE_HOOK[] !== nothing && TD.STAGE_HOOK[](:be, Δ)
            s!(:thet, θ0)
            if ok; STATS[:be_sticky] = get(STATS, :be_sticky, 0) + 1; break; end
            restore!(); TD.msg(1)
            s!(:tfac, TFAC[]); s!(:dtime, TFAC[] * Δ); s!(:jtry, g(:jtry) + 1); g(:jtry) > g(:ntry) && TD.msg(2)
            r(:reset); s!(:timen, t0 + g(:dtime)); TitanCore.predict!(; halve = false); continue
        end
        s!(:thet, 0.5); s!(:dtime, γ * Δ); s!(:timen, t0 + γ * Δ)
        ok = newton!()
        if !ok
            STATS[:fail_tr] += 1
        else
            TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](:tr, γ * Δ)
            if bdf_old!()
                s!(:thet, 1.0); s!(:dtime, C2 * Δ); s!(:timen, t0 + Δ)
                TitanCore.predict!(; halve = false)       # initial guess y* + C2 Δ f from the grid velocities
                ok = newton!()
                ok ? (TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](:bdf, C2 * Δ)) : (STATS[:fail_bdf] += 1)
            else
                ok = false; STATS[:fail_bdf] += 1
            end
        end
        restore!(); s!(:thet, θ0); s!(:dtime, Δ); s!(:timen, t0 + Δ)
        ok && break
        if FALLBACK_BE[]          # same Δ with backward Euler from y_n before cutting the step
            r(:reset); s!(:thet, 1.0); s!(:timen, t0 + Δ); TitanCore.predict!(; halve = false)
            ok = newton!()
            ok && TD.STAGE_HOOK[] !== nothing && TD.STAGE_HOOK[](:be, Δ)
            s!(:thet, θ0)
            ok && (STATS[:be_fallback] += 1; BE_LEFT[] = STICKY[]; break)
        end
        TD.msg(1)
        s!(:tfac, TFAC[]); s!(:dtime, TFAC[] * Δ)
        s!(:jtry, g(:jtry) + 1); g(:jtry) > g(:ntry) && TD.msg(2)
        r(:reset); s!(:timen, t0 + g(:dtime)); TitanCore.predict!(; halve = false, dtscale = γ)   # next: the TR stage
    end
    STATS[:steps] += 1
    r(:totnrg)
    r(:timstp)
    nothing
end


# ================================================================== DTERR: step-size control from the embedded error estimate
"""
DTERR=rtol[,xfloor[,ffloor[,efloor]]] (ENV, default off; off = bit-identical to the scheme above): TR-BDF2's embedded local-error
estimate is the ONLY step-size control (`step_err!`). TITAN's heuristic (timstp: stol on relative changes, CDYN, DT_HOOKS — DTENV,
the COUPLED x_H⁺ cap, DTRADMAX —, the line-cooling cap and STOLCHEM in popiii_m2c) no longer sets Δt; timstp still runs, and the
Δt it would have chosen is kept in `HEUR[]` for comparison.

**Estimate** (Hosea & Shampine 1996, Appl. Numer. Math. 20, 21–37, §2 and §6; the pair of Bank et al. 1985, IEEE Trans. CAD 4, 436).
TR-BDF2 is the FSAL DIRK pair (d = γ/2, w = √2/4, γ = 2 − √2)

    0 |  0  0  0           b  = (w, w, d)                        (order 2, the solution)
    γ |  d  d  0           b̂  = ((1−w)/3, (3w+1)/3, d/3)          (order 3, the companion)
    1 |  w  w  d

with stage derivatives z₁ = h f(yₙ), z₂ = h f(y_γ), z₃ = h f(yₙ₊₁), and est = ŷₙ₊₁ − yₙ₊₁ = Σ (b̂ᵢ − bᵢ) zᵢ (checked against the paper's
Butcher array, p. 23, and SUNDIALS' ARKODE table TRBDF2-3-3-2). Since d = (1−γ)/(2−γ) = C2, the BDF2 stage above *is* the third DIRK
stage, and the stages give z₂ = (y_γ − yₙ)/d − z₁, z₃ = (yₙ₊₁ − y*)/d (y* = A y_γ − B yₙ, the "old" state of the BDF stage), so

    est = −(4w/3) z₁ + (y_γ − yₙ)/(3d) − (2/3)(yₙ₊₁ − y*).

z₁ is the FSAL "smoothed scaled derivative" of Hosea & Shampine §6 (their z̃ = (yₙ − ψ)/d): z₁ = h·Dₙ with Dₙ = (yₙ − y*ₙ₋₁)/(d hₙ₋₁)
from the last accepted step (for an accepted backward-Euler step, Dₙ = (yₙ − yₙ₋₁)/hₙ₋₁). After a restart there is no Dₙ: that first
step is accepted without an estimate and Δt is kept. All y here are TITAN's conserved per-zone combinations (m, m e_g, m e_gr, E ΔV,
momentum, the flux-row inertia, x_me, r, ρ x ΔV of the coupled species), the variables in which the residual is linear in the old state.
**Filter** (Hosea & Shampine §6, "(I − h d J) Est = est"): the raw est grossly over-estimates stiff components (∝ hλ), so the estimate
used is Est = (I − h d J)⁻¹ est. That matrix is the BDF stage's Newton matrix: the old arrays are set to y* + est (`old_from!`) and a
matgen/matslv at the converged yₙ₊₁ gives Est, less the solve with the old arrays at y* (the converged state's own remaining
correction: F_r is not in TITAN's Newton convergence measure, so without this the estimate has an h-independent part), directly in TITAN's correction variables (relative r, x_me, ρ, T, E_r; u/unom, F/frnom;
δ ln x of the coupled species).

**Second filter** (DTERRFILT=2, default; 1 = Est only): EST = (I − h d J)⁻¹ Est, which Hosea & Shampine §6 note gives the correct
behaviour at infinity (Est still over-estimates a stiff component by the constant factor √2/(3d) instead of ∝ 1/(hλ)); see refilter!.
On the 100-zone grid the stiff F_r/E_r of the thin zones outside the photosphere made Est nearly h-independent (err 1.0–1.07 over
dt 2.8e3 → 1.5e3 s, every attempt rejected).
**Norm** (DTERRNORM=rms (default) | max): err = ‖Est / w‖ with per-variable weights
    ρ, T: rtol·|y|;  E_r: rtol·(E_r + efloor·aT⁴);  u: rtol·max(|u|, c_s) (c_s the larger adjacent sound speed);
    F_r: rtol·max(|F|, ffloor·c E_r);  x_me: rtol·max(x_me, adjacent zone masses) (TITAN's own measure, relative change of the mass coordinate; the zone masses
    themselves are controlled through ρ — weighting δx_me by the zone mass, a Courant-type test, held the production state, where the
    accretion shock crosses the fixed interfaces, at Δt ~ 1e5 s against the heuristic's 1.3e8 s);  r (moving grids): rtol·(smaller adjacent zone width);
    x_H⁺, x_H₂ (COUPLED zones): rtol·(x + xfloor) on |δx| = |δ ln x|·x.
Defaults: xfloor 1e-6, ffloor 1 (F_r measured on its causal scale cE: with 1e-2 the thin zones just outside the photosphere, where F crosses through small values, held Δt at ~400 s), efloor 0.

**Controller** (DTERRCTL=safety,facmin,facmax, default 0.9,0.2,5; DTERRPI=1): accepted when err ≤ 1, then
h_new = h·clamp(safety·err^(−0.7/3)·err_prev^(0.4/3), facmin, facmax) (PI control, Gustafsson 1991 / Hairer & Wanner II §IV.8; DTERRPI=0:
safety·err^(−1/3)); no growth on a step that followed a rejection. err > 1: the step is rejected and retried from yₙ with
h·max(facmin, safety·err^(−1/3)). A Newton failure in either stage is a rejection with h·DTERRFAIL (default 0.25). No backward-Euler
fallback and no sticky BE phase (DTERRBE=1 restores a non-sticky BE retry of the same h after a Newton failure, for testing).
COUPLED species and cooling are backward Euler also in the TR stage (TitanResidFast.CPL_THETA = 1): for those components the scheme is
first order and the estimate is O(h²) instead of O(h³) — it still measures their error, just more conservatively.
DTERRMAX=<s> (test only, default none): an upper bound on Δt, to separate step-size from history effects.
DTERRTRACE=n (default 20): a line every n steps (Δt, err, the limiting variable and zone, the heuristic Δt), and every rejection.
"""
const ERRCTL = Ref(haskey(ENV, "DTERR"))
"""
DT_BOUNDS: hard bounds on the next step that the error control respects (functions () -> Δt_max [s], evaluated after an accepted
step, with the step's Δt in `dtime` and its old/new state in the arrays): the next Δt is min(error-controlled Δt, bounds). Not an
accuracy criterion but a validity limit of a discretisation, e.g. GRIDEQ's mesh motion (no interface moves more than fmove × its
smaller neighbouring zone per step, popiii_m2c). TITAN's DT_HOOKS stay diagnostic under DTERR; ESTATS[:bound] counts the steps a
bound set.
"""
const DT_BOUNDS = Function[]
const EOPT = let v = haskey(ENV, "DTERR") ? parse.(Float64, split(ENV["DTERR"], ",")) : Float64[], c = parse.(Float64, split(get(ENV, "DTERRCTL", "0.9,0.2,5"), ","))
    (rtol = isempty(v) ? 1e-2 : v[1], xfl = length(v) >= 2 ? v[2] : 1e-6, ffl = length(v) >= 3 ? v[3] : 1.0, efl = length(v) >= 4 ? v[4] : 0.0,
     safety = c[1], facmin = c[2], facmax = c[3], pi = get(ENV, "DTERRPI", "1") == "1", fail = parse(Float64, get(ENV, "DTERRFAIL", "0.25")),
     hmax = parse(Float64, get(ENV, "DTERRMAX", "Inf")), be = get(ENV, "DTERRBE", "0") == "1", filt = parse(Int, get(ENV, "DTERRFILT", "2")), maxnorm = get(ENV, "DTERRNORM", "rms") == "max", trace = parse(Int, get(ENV, "DTERRTRACE", "20")))
end
const dW = γ / 2; const wW = sqrt(2) / 4
const ENTRY = Ref((-1, NaN)); const CSAVE = Ref{Any}(nothing); const LASTACC = Ref(0.0)
const HNEXT = Ref(0.0); const ERRPREV = Ref(0.0); const HEUR = Ref(0.0); const DSAL = Ref{Any}(nothing)
const ESTATS = Dict(:steps => 0, :rej_err => 0, :fail_newton => 0, :nofsal => 0, :be => 0, :neg_old => 0, :nonfinite => 0, :nonfinite_est => 0)
const LIMHIST = Dict{String,Int}(); const LASTERR = Ref((NaN, "", 0)); const ERR0 = Ref(NaN); const ERR1 = Ref(NaN); const DIAGZ = Ref(isempty(get(ENV, "DTERRDIAGZ", "")) ? Int[] : parse.(Int, split(ENV["DTERRDIAGZ"], ",")))   # (DTERRDIAGZ=z1,z2,..: zones printed by DTERRDIAG besides the limiting one)
const DIAG = Ref(parse(Int, get(ENV, "DTERRDIAG", "0")))   # (DTERRDIAG=n: stage derivatives of the limiting zone every n steps)   # (ERR0: norm of the converged state's remaining correction)
"Extra text for the DTERR trace line (popiii_m2c: coverage of the COUPLED zones)."
const TRACE_EXTRA = Ref{Any}(nothing)
const VNAMES = Dict(1 => "r", 2 => "m", 3 => "rho", 4 => "u", 5 => "T", 6 => "E_r", 7 => "F_r", 8 => "x_H+", 9 => "x_H2")

"Conserved combinations of a state (arrays of the new or the old time level; xs = the two coupled species)."
function combos(d, dv, eg, egr, er, u, fr, xme, rr, xs, mg = nothing)
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1; n = length(d)
    C = Dict{Symbol,Vector{Float64}}(s => zeros(n) for s in (:m, :me, :xme, :r, :P))
    m = C[:m]
    @inbounds for k in a:b; m[k] = d[k] * dv[k]; C[:me][k] = m[k] * eg[k]; end
    @inbounds for k in a-1:b+2; C[:xme][k] = xme[k]; C[:r][k] = rr[k]; end
    @inbounds for k in a+1:b; C[:P][k] = u[k] * (m[k-1] + m[k]); end
    TitanCore.FREE_OUTER[] && (C[:P][b+1] = u[b+1] * m[b])
    if lrad
        C[:mer] = zeros(n); C[:EV] = zeros(n); C[:QF] = zeros(n)
        @inbounds for k in a:b; C[:mer][k] = m[k] * egr[k]; C[:EV][k] = er[k] * dv[k]; end
        @inbounds for k in a+1:b; C[:QF][k] = fr[k] * 2 * (m[k-1] + m[k]) / (d[k-1] + d[k]); end
    end
    if TX.COUPLED[]
        for i in 1:2
            x = zeros(n); C[i == 1 ? :x1 : :x2] = x
            @inbounds for k in a:b; TX.cplz(k) && (x[k] = m[k] * xs[i][k]); end
        end
    end
    if mg !== nothing                  # MGROUPS: E_g ΔV and the flux-row inertia of F_g
        for g in 1:TX.MGN[]
            kE, kF = mgkeys(g); EV = zeros(n); QF = zeros(n); C[kE] = EV; C[kF] = QF; Eg = mg[1][g]; Fg = mg[2][g]
            @inbounds for k in a:b; EV[k] = Eg[k] * dv[k]; end
            @inbounds for k in a+1:b; QF[k] = Fg[k] * 2 * (m[k-1] + m[k]) / (d[k-1] + d[k]); end
        end
    end
    C
end
combos_new() = combos(vf(:dn), vf(:dvoln), vf(:egn), vf(:egrn), vf(:ern), vf(:un), vf(:frn), vf(:xmen), vf(:rn), TX.XSN[], TX.MG[] ? (TX.MGEN[], TX.MGFN[]) : nothing)
combos_old() = combos(vf(:do), vf(:dvolo), vf(:ego), vf(:egro), vf(:ero), vf(:uo), vf(:fro), vf(:xmeo), vf(:ro), TX.XSO[], TX.MG[] ? (TX.MGEO[], TX.MGFO[]) : nothing)
combos_saved() = combos(SAVE[:do], SAVE[:dvolo], SAVE[:ego], SAVE[:egro], SAVE[:ero], SAVE[:uo], SAVE[:fro], SAVE[:xmeo], SAVE[:ro], SAVEX, TX.MG[] ? (SAVEMG[1], SAVEMG[2]) : nothing)

"Old arrays := the state with conserved combinations C (as bdf_old!, for any combination). false if a mass or volume is ≤ 0."
function old_from!(C)
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1
    lag = Int(g(:lgrid)) == 3 || TitanCore.GRIDMOVE[]
    do_, ego, uo, xmeo, dvo = vf(:do), vf(:ego), vf(:uo), vf(:xmeo), vf(:dvolo)
    ms = C[:m]
    @inbounds for k in a:b; ms[k] > 0 || return false; end
    if lag
        ro = vf(:ro); mu1 = Int(g(:mu)) + 1; xmup1 = g(:xmup1)
        @inbounds for k in a-1:b+2; ro[k] = C[:r][k]; end
        @inbounds for k in a:b; dvo[k] = (ro[k+1]^mu1 - ro[k]^mu1) / xmup1; dvo[k] > 0 || return false; end
    end
    @inbounds for k in a:b; ego[k] = C[:me][k] / ms[k]; do_[k] = ms[k] / dvo[k]; end
    @inbounds for k in a-1:b+2; xmeo[k] = C[:xme][k]; end
    if TX.COUPLED[]
        for i in 1:2
            xo_ = TX.XSO[][i]; x = C[i == 1 ? :x1 : :x2]
            @inbounds for k in a:b; TX.cplz(k) && (xo_[k] = x[k] / ms[k]); end
        end
    end
    if lrad
        egro, ero, fro = vf(:egro), vf(:ero), vf(:fro)
        @inbounds for k in a:b; egro[k] = C[:mer][k] / ms[k]; ero[k] = C[:EV][k] / dvo[k]; end
        @inbounds for k in a+1:b; fro[k] = C[:QF][k] * (do_[k-1] + do_[k]) / (2 * (ms[k-1] + ms[k])); end
        if TX.MG[]
            for g in 1:TX.MGN[]
                kE, kF = mgkeys(g); eo_ = TX.MGEO[][g]; fo_ = TX.MGFO[][g]
                @inbounds for k in a:b; eo_[k] = C[kE][k] / dvo[k]; end
                @inbounds for k in a+1:b; fo_[k] = C[kF][k] * (do_[k-1] + do_[k]) / (2 * (ms[k-1] + ms[k])); end
            end
        end
    end
    @inbounds for k in a+1:b; uo[k] = C[:P][k] / (ms[k-1] + ms[k]); end
    TitanCore.FREE_OUTER[] && (uo[b+1] = C[:P][b+1] / ms[b])
    true
end


"(SDZONES) per zone: weighted estimate entries (units of the tolerance: 1 = rtol) of E_r, x_H+, x_H2, T, rho, F_r for the first (f1) and second (f2) filter, and the three stage increments z1,z2,z3 of the zone's E_r*dV, x_H+ mass."
function sd_est_print(Δ, C1, Cs, Cg, Cn, D)
    a = Int(g(:ngrs)); ie, it, id, iff = Int(g(:ie)), Int(g(:it)), Int(g(:id)), Int(g(:if))
    for zk in SDZ
        k = a + zk - 1
        e(t, i) = get(SDACC, (t, zk, i), NaN)
        zs = ""
        for s_ in (:EV, :x1)
            haskey(C1, s_) || continue
            z1 = D === nothing ? NaN : Δ * get(D, s_, zeros(k))[k]; z3 = (C1[s_][k] - Cs[s_][k]) / dW; z2 = (Cg[s_][k] - Cn[s_][k]) / dW - z1
            zs *= @sprintf(" %s[y_n %.4e z1 %.3e z2 %.3e z3 %.3e]", s_, Cn[s_][k], z1, z2, z3)
        end
        @printf("SDEST jstep %d dt %.4e zone %d | f2: Er %.3e xHp %.3e xH2 %.3e T %.3e rho %.3e F %.3e | f1: Er %.3e xHp %.3e xH2 %.3e T %.3e rho %.3e |%s\n", Int(g(:jstep)), Δ, zk,
                e(:f2, ie), e(:f2, 8), e(:f2, 9), e(:f2, it), e(:f2, id), e(:f2, iff), e(:f1, ie), e(:f1, 8), e(:f1, 9), e(:f1, it), e(:f1, id), zs)
    end
end
"(SDZONES) state of the accepted step's zones: E_r, T, x_H+, x_H2, rho, u, F_r, and the regime flags (SD_FLAGS)."
function sd_state_print(Δ)
    a = Int(g(:ngrs))
    TitanResidFast.SDDEC[] && println("SDDEC jstep ", Int(g(:jstep)), " flips since last print (zone x count):", TitanResidFast.sd_dec_text(a-1 .+ SDZ), "  CPL_BAD=", sort!(collect(TitanResidFast.CPL_BAD)))
     er, tn, dn, un, fr = vf(:ern), vf(:tn), vf(:dn), vf(:un), vf(:frn)
    rn, chi, asn, pg, b_ = vf(:rn), vf(:chifn), vf(:asn), vf(:pgn), Int(g(:ngre))
    tau = zeros(length(dn)); tt = 0.0; for k in b_:-1:a; tt += chi[k] * dn[k] * (rn[k+1] - rn[k]); tau[k] = tt; end     # optical depth from the outer boundary to the zone's inner face
    for zk in SDZ
        k = a + zk - 1
        x1 = TX.COUPLED[] ? TX.XSN[][1][k] : NaN; x2 = TX.COUPLED[] ? TX.XSN[][2][k] : NaN
        TitanResidFast.SDHD && haskey(TitanResidFast.SDHRAW, k) && (hr = TitanResidFast.SDHRAW[k]; @printf("SDHD jstep %d zone %d interface-donor raw_h %.6e hFS %.6e tau_i %.5e a %.5e\n", Int(g(:jstep)), zk, hr...))
        TitanResidFast.SDROW && haskey(TitanResidFast.SDROWV, k) && (rv = TitanResidFast.SDROWV[k]; @printf("SDROW jstep %d zone %d R_mom %.4e B %.4e w %.4f | dF: abs %.4e inert %.4e tie %.4e | dE_(k-1): mom %.4e tie %.4e\n", Int(g(:jstep)), zk, rv[1], rv[2], rv[3], rv[4], rv[8], rv[5], rv[6], rv[7]))
        @printf("SDSTATE jstep %d t %.9e yr dt %.4e zone %d Er %.9e T %.9e xHp %.9e xH2 %.9e rho %.9e u %.6e F %.6e r %.6e dr %.4e chi %.4e tau %.4e cs %.4e p %.4e %s\n", Int(g(:jstep)), g(:timen) / 3.156e7, Δ, zk, er[k], tn[k], x1, x2, dn[k], un[k], fr[k], rn[k], rn[k+1] - rn[k], chi[k], tau[k], asn[k], pg[k],
                SD_FLAGS[] === nothing ? "" : SD_FLAGS[](k))
    end
end
"""
Filtered error estimate of the step just solved (BDF stage converged: new = yₙ₊₁, old = y*, θ = 1, dtime = dΔ). Returns
(err, variable name, zone from the inner boundary); the old arrays are left at y* + est (the caller restores yₙ).
"""
function estimate!(Cn, Cg, C1, Cs, Δ)
    D = DSAL[]
    Co = Dict{Symbol,Vector{Float64}}()
    for (s, v1) in C1
        vs = Cs[s]; vg = Cg[s]; vn = Cn[s]; Dv = get(D, s, nothing)
        o = similar(v1)
        @inbounds for k in eachindex(v1)
            z1 = (Dv === nothing || k > length(Dv)) ? 0.0 : Δ * Dv[k]
            o[k] = vs[k] + (-(4wW / 3) * z1 + (vg[k] - vn[k]) / (3dW) - (2 / 3) * (v1[k] - vs[k]))
        end
        Co[s] = o
    end
    # the converged state's own remaining Newton correction (at old = y*; F_r is not in TITAN's convergence measure, so it need not
    # be small) is subtracted: Est = δ(y* + est) − δ(y*), the same matrix (linear in the old state)
    rhs = TF.var(:rhs)
    TitanCore.CSENS_SKIP[] = TD.CSENS.valid
    try
        s!(:iter, 1); r(:matgen); r(:matslv); SDESTNK && TD.KRYLOV[] !== nothing && TD.KRYLOV[](); δ0 = copy(rhs)
        old_from!(Co) || (ESTATS[:neg_old] += 1; return (Inf, "m<=0", 0, Inf))
        TX.RESID_ONLY[] = true; BlockSolver.Penta.REUSE_LU[] = true
        r(:matgen); r(:matslv); SDESTNK && TD.KRYLOV[] !== nothing && TD.KRYLOV[](); rhs .-= δ0
        SDTAG[] = :e0; e0 = wnorm(δ0); ERR0[] = e0[1]
        SDTAG[] = :f1
        EOPT.filt >= 2 || return wnorm(rhs)
        ERR1[] = wnorm(rhs)[1]
        refilter!(copy(rhs), Cs, δ0) || (ESTATS[:neg_old] += 1; return (Inf, "m<=0", 0, Inf))
    finally
        TitanCore.CSENS_SKIP[] = false; TX.RESID_ONLY[] = false; BlockSolver.Penta.REUSE_LU[] = false
    end
    SDTAG[] = :f2
    res = wnorm(rhs)
    SDON && sd_est_print(Δ, C1, Cs, Cg, Cn, D)
    if DIAG[] > 0 && ESTATS[:steps] % DIAG[] == 0 && res[3] > 0       # stage derivatives of the limiting zone: smooth flow has z1 ≈ z2 ≈ z3
        for zd in unique(vcat(res[3], DIAGZ[]))      # (DTERRDIAGZ=z1,z2,..: these zones too, in addition to the limiting one)
          k = Int(g(:ngrs)) + zd - 1
          for s_ in (:EV, :me, :m, :P, :QF, :x1, :x2)
            haskey(C1, s_) || continue
            z1 = D === nothing ? NaN : Δ * get(D, s_, zeros(k))[k]; z3 = (C1[s_][k] - Cs[s_][k]) / dW; z2 = (Cg[s_][k] - Cn[s_][k]) / dW - z1
            @printf("DTERRDIAG jstep %d zone %d %s: y_n %.4e  z1 %.3e z2 %.3e z3 %.3e  est %.3e\n", Int(g(:jstep)), zd, s_, Cn[s_][k], z1, z2, z3, Co[s_][k] - Cs[s_][k])
          end
        end
    end
    res
end

"""
Second filter (DTERRFILT=2, default; Hosea & Shampine 1996 §6, after their eq. for Est: "(I − h d J) EST = Est … to get the correct
behavior at infinity"): rhs := (I − h d J)⁻¹ Est. Est is in TITAN's correction variables; its conserved-combination image ΔQ = Q_y Est is
taken from the residual's own state evaluation (TitanResidFast.residual! at the converged state + ε Est and at + 0, difference / ε) and
injected through the old arrays as for the first filter.
"""
function refilter!(E1, Cs, δ0)
    c = TX.context(); Wws = TX.getws(Float64, c.mgr); W = Wws.a; neqn, k0, nb = c.neqn, c.k0, c.nb; R = zeros(neqn, nb)
    x = zeros(neqn * nb); mx = 0.0
    @inbounds for bb in 1:nb, i in 1:neqn; v = E1[i, k0+bb-1]; isfinite(v) || return false; x[(bb-1)*neqn+i] = v; mx = max(mx, abs(v)); end
    ε = mx > 0 ? min(1.0, 1e-3 / mx) : 1.0
    Wm = TX.MG[] ? TX.mgws(Float64, c.mgr) : nothing
    stw() = combos(W.dn, W.dvoln, W.egn, W.egrn, W.ern, W.un, W.frn, W.xmen, W.rn, (W.xs1, W.xs2), Wm === nothing ? nothing : (Wm.egn, Wm.fgn))
    TX.residual!(R, Wws, c, zeros(neqn * nb)); C0 = Dict(s_ => copy(v) for (s_, v) in stw())
    x .*= ε; TX.residual!(R, Wws, c, x); Cp = stw()
    Co = Dict{Symbol,Vector{Float64}}()
    for (s_, v0) in C0
        vs = Cs[s_]; vp = Cp[s_]; o = similar(v0)
        @inbounds for k in eachindex(v0); o[k] = vs[k] + (vp[k] - v0[k]) / ε; end
        Co[s_] = o
    end
    old_from!(Co) || return false
    rhs = TF.var(:rhs); r(:matgen); r(:matslv); SDESTNK && TD.KRYLOV[] !== nothing && TD.KRYLOV[](); rhs .-= δ0
    true
end

"Weighted norm of a correction matrix (TITAN's correction variables, see ERRCTL): (norm, variable, zone, largest weighted entry)."
function wnorm(rhs)
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1; lag = Int(g(:lgrid)) == 3 || TitanCore.GRIDMOVE[]
    ir, im, id, iu, it, ie, iff = (Int(g(s)) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    rn, xmen, dn, un, tn, unom, asn = vf(:rn), vf(:xmen), vf(:dn), vf(:un), vf(:tn), vf(:unom), vf(:asn)
    dvn = vf(:dvoln)
    o = EOPT; rt = o.rtol; ssq = 0.0; nn = 0; worst = 0.0; wv = 0; wk = 0
    acc(v, i, k) = (SDON && (k - a + 1) in SDZ && (SDACC[(SDTAG[], k - a + 1, i)] = v); isfinite(v) || (v = Inf); a2 = abs(v); ssq += a2^2; nn += 1; a2 > worst && (worst = a2; wv = i; wk = k - a + 1))
    @inbounds for k in a:b
        acc(rhs[id, k] / rt, id, k)
        acc(rhs[it, k] / rt, it, k)
        if lrad
            E = vf(:ern)[k]; acc(rhs[ie, k] * E / (rt * (E + o.efl * 7.565723e-15 * tn[k]^4)), ie, k)
        end
        if TX.COUPLED[] && TX.cplz(k)
            for (i, q) in ((8, TX.XSN[][1]), (9, TX.XSN[][2]))
                acc(rhs[i, k] * q[k] / (rt * (q[k] + o.xfl)), i, k)
            end
        end
        if TX.MG[]                     # MGROUPS: δE_g against the group's own E_g with the grey floor plus MGEFL·E_r (minor groups do not set Δt)
            E = vf(:ern)[k]
            for g in 1:TX.MGN[]
                Eg = TX.MGEN[][g][k]; acc(rhs[TX.mgiE(g), k] * Eg / (rt * (Eg + o.efl * 7.565723e-15 * tn[k]^4 + TX.MGEFL[] * E)), TX.mgiE(g), k)
            end
        end
    end
    @inbounds for k in a+1:b
        cs = max(asn[k-1], asn[k]); acc(rhs[iu, k] * unom[k] / (rt * max(abs(un[k]), cs)), iu, k)
        mz = max(abs(xmen[k]), dn[k-1] * dvn[k-1], dn[k] * dvn[k]); acc(rhs[im, k] * xmen[k] / (rt * mz), im, k)
        lag && acc(rhs[ir, k] * rn[k] / (rt * min(rn[k] - rn[k-1], rn[k+1] - rn[k])), ir, k)
        if lrad
            frn, ern, frnom = vf(:frn), vf(:ern), vf(:frnom)
            acc(rhs[iff, k] * frnom[k] / (rt * max(abs(frn[k]), o.ffl * 2.99792458e10 * 0.5 * (ern[k-1] + ern[k]))), iff, k)
            if TX.MG[]
                for g in 1:TX.MGN[]
                    Fg = TX.MGFN[][g][k]
                    acc(rhs[TX.mgiF(g), k] * TX.MGFNOM[][g][k] / (rt * max(abs(Fg), o.ffl * 2.99792458e10 * 0.5 * (ern[k-1] + ern[k]), TX.MGEFL[] * abs(frn[k]))), TX.mgiF(g), k)
                end
            end
        end
    end
    err = o.maxnorm ? worst : sqrt(ssq / max(nn, 1))
    (err, get(VNAMES, wv, string(wv)), wk, worst)
end

"All unknowns of the new time level finite (r, x_me, u, ρ, T, E_r, F_r and the coupled species)."
function finite_new()
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1
    for s_ in (lrad ? (:rn, :xmen, :un, :dn, :tn, :ern, :frn) : (:rn, :xmen, :un, :dn, :tn))
        v = vf(s_); @inbounds for k in a:b; isfinite(v[k]) || return false; end
    end
    if TX.COUPLED[]
        for i in 1:2; q = TX.XSN[][i]; @inbounds for k in a:b; TX.cplz(k) && !isfinite(q[k]) && return false; end; end
    end
    if TX.MG[]
        for g in 1:TX.MGN[], q in (TX.MGEN[][g], TX.MGFN[][g]); @inbounds for k in a:b; isfinite(q[k]) || return false; end; end
    end
    true
end

function etrace(kind, Δ, err, info)
    TRACE_EXTRA[] === nothing ? "" : TRACE_EXTRA[]()
end

# ================================================================== DTSCHEME: alternative one-step schemes under the DTERR controller
"""
DTSCHEME=trbdf2 (default; the scheme above, unchanged) | sdirk2 | be (ENV, only with DTERR; a test switch for the photospheric
fluctuation study, branch photfluct). Both alternatives are L-stable and have no non-dissipative (θ = ½) stage:

* `sdirk2`: Alexander's two-stage, second-order, stiffly accurate SDIRK (Alexander 1977, SIAM J. Numer. Anal. 14, 1006, γ = 1 − 1/√2;
  Hairer & Wanner II §IV.6). Stage 1: backward Euler over γΔ from yₙ (Y₁). Stage 2: y_{n+1} = yₙ + (1−γ)Δ f(Y₁) + γΔ f(y_{n+1}),
  i.e. backward Euler over γΔ from the "old" state y** = yₙ + (1−γ)/γ (Y₁ − yₙ) = Y₁ + √2 (Y₁ − yₙ), formed on the same conserved
  combinations as TR-BDF2's y* (`old_from!`). Embedded first-order solution ŷ = yₙ + Δ f(Y₁): est = y_{n+1} − ŷ = γΔ (f₂ − f₁)
  = (y_{n+1} − y**) − (Y₁ − yₙ) (no FSAL derivative needed), filtered with the stage-2 Newton matrix (I − γΔJ) as for TR-BDF2.
* `be`: backward Euler, est = ½[(y_{n+1} − yₙ) − Δ f(yₙ)] (f(yₙ) = (yₙ − yₙ₋₁)/Δₙ₋₁ of the previous BE step, FSAL), filtered with (I − ΔJ).

Both estimates are O(Δ²): the controller uses exponents 1/2 (PI: 0.7/2, 0.4/2) instead of 1/3. The stage hook books stage 1 as `:be`
and stage 2 of sdirk2 as `:sd2` (weights 1, linear burn path): the nuclear-energy ledger is approximate for sdirk2 (no burning at the
M* ≲ 1 M_sun states it is meant for).
"""
const SCHEME = Ref(Symbol(get(ENV, "DTSCHEME", "trbdf2")))
const γS = 1 - 1 / sqrt(2); const BS = sqrt(2.0)
combine(C1, C2_, f) = Dict(s_ => f(C1[s_], C2_[s_]) for s_ in keys(C1))

"Filtered error estimate for DTSCHEME ≠ trbdf2: old arrays := Cs + est at the converged new state (θ = 1, the last stage's dtime)."
function estimate_alt!(est, Cs)
    Co = combine(Cs, est, (a_, e_) -> a_ .+ e_)
    rhs = TF.var(:rhs)
    TitanCore.CSENS_SKIP[] = TD.CSENS.valid
    try
        s!(:iter, 1); r(:matgen); r(:matslv); δ0 = copy(rhs)
        old_from!(Co) || (ESTATS[:neg_old] += 1; return (Inf, "m<=0", 0, Inf))
        TX.RESID_ONLY[] = true; BlockSolver.Penta.REUSE_LU[] = true
        r(:matgen); r(:matslv); rhs .-= δ0
        ERR0[] = wnorm(δ0)[1]
        EOPT.filt >= 2 || return wnorm(rhs)
        ERR1[] = wnorm(rhs)[1]
        refilter!(copy(rhs), Cs, δ0) || (ESTATS[:neg_old] += 1; return (Inf, "m<=0", 0, Inf))
        return wnorm(rhs)
    finally
        TitanCore.CSENS_SKIP[] = false; TX.RESID_ONLY[] = false; BlockSolver.Penta.REUSE_LU[] = false
    end
end

"One attempt of the step Δ from yₙ (in SAVE; combos Cn) with DTSCHEME sdirk2 or be: (ok, err, info, Dnew)."
function attempt_alt!(Δ, t0, Cn)
    err = NaN; info = ("", 0, NaN); Dnew = nothing
    hook(k_, dt_) = TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](k_, dt_)
    fin(ok_) = (ok_ && !finite_new()) ? (ESTATS[:nonfinite] += 1; false) : ok_
    function est_of(est, Cs)
        e_, v_, z_, w_ = try estimate_alt!(est, Cs) catch ex
            ex isa InterruptException && rethrow()
            (ESTATS[:nonfinite_est] += 1) <= 3 && println("DTERR: estimate failed: ", first(sprint(showerror, ex, catch_backtrace()), 1500))
            (Inf, "estimate failed", 0, Inf)
        end
        isfinite(e_) || (e_ = Inf; ESTATS[:nonfinite_est] += 1)
        (e_, (v_, z_, w_))
    end
    if SCHEME[] === :be
        s!(:thet, 1.0); s!(:dtime, Δ); s!(:timen, t0 + Δ)
        ok = fin(newton!())
        ok || (STATS[:fail_tr] += 1; return (false, err, info, Dnew))
        hook(:be, Δ)
        C1 = combos_new(); Cs = combos_old()
        Dnew = combine(C1, Cn, (a_, b_) -> (a_ .- b_) ./ Δ)
        D = DSAL[]
        if D === nothing
            ESTATS[:nofsal] += 1
        else
            est = Dict(s_ => (C1[s_] .- Cn[s_] .- Δ .* (haskey(D, s_) && length(D[s_]) == length(C1[s_]) ? D[s_] : zero(C1[s_]))) ./ 2 for s_ in keys(C1))
            err, info = est_of(est, Cs)
        end
        return (true, err, info, Dnew)
    end
    # sdirk2
    s!(:thet, 1.0); s!(:dtime, γS * Δ); s!(:timen, t0 + γS * Δ)
    ok = fin(newton!())
    ok || (STATS[:fail_tr] += 1; return (false, err, info, Dnew))
    hook(:be, γS * Δ)
    Cg = combos_new()
    Cs = combine(Cg, Cn, (y1, yn) -> y1 .+ BS .* (y1 .- yn))
    old_from!(Cs) || (STATS[:fail_bdf] += 1; ESTATS[:neg_old] += 1; return (false, err, info, Dnew))
    s!(:dtime, γS * Δ); s!(:timen, t0 + Δ)
    TitanCore.predict!(; halve = false)
    ok = fin(newton!())
    ok || (STATS[:fail_bdf] += 1; return (false, err, info, Dnew))
    hook(:sd2, γS * Δ)
    C1 = combos_new()
    Dnew = combine(C1, Cs, (a_, b_) -> (a_ .- b_) ./ (γS * Δ))
    est = Dict(s_ => (C1[s_] .- Cs[s_]) .- (Cg[s_] .- Cn[s_]) for s_ in keys(C1))
    err, info = est_of(est, Cs)
    (true, err, info, Dnew)
end

"TR-BDF2 step with error control (DTERR): see the docstring of ERRCTL."
function step_err!()
    Int(g(:lgrid)) in (2, 3) || error("Trbdf2: lgrid = 2 or 3 only")
    o = EOPT; BE_LEFT[] = 0
    alt = SCHEME[] !== :trbdf2; γ1 = SCHEME[] === :be ? 1.0 : SCHEME[] === :sdirk2 ? γS : γ; qe = alt ? 2 : 3
    # A repeat of the same step (StarMaterial.solve_step!, HOT: the burn/hydro source iteration restores the snapshot and calls the
    # step again): same controller state as at the first call, and the Δt that call accepted
    key = (Int(g(:jstep)), g(:timen))
    if key == ENTRY[] && CSAVE[] !== nothing
        HNEXT[], ERRPREV[], DSAL[] = CSAVE[]; s!(:dtime, LASTACC[]); ESTATS[:repeat] = get(ESTATS, :repeat, 0) + 1
    else
        ENTRY[] = key; CSAVE[] = (HNEXT[], ERRPREV[], DSAL[])
        HNEXT[] > 0 && s!(:dtime, HNEXT[])
    end
    SDFIXDT > 0 && s!(:dtime, SDFIXDT)
    Δ0 = g(:dtime); s!(:dtime, γ1 * Δ0)
    r(:peq)
    s!(:dtime, g(:dtime) / γ1); s!(:timen, g(:timeo) + g(:dtime))
    θ0 = g(:thet); t0 = g(:timeo)
    save!(); Cn = combos_saved()
    s!(:jtry, 1); s!(:jback, 0)
    rejected = false; err = NaN; info = ("", 0, NaN); Dnew = nothing; nrej = 0
    if !isempty(NWPROBE[])        # (NWPROBE, diagnostic, default off) one trapezoidal-stage attempt at each given step size from this state, as a retry after a failure would start it
        NWDUMP[] === nothing || NWDUMP[]()
        for dtp in NWPROBE[]
            restore!(); s!(:thet, θ0); s!(:dtime, dtp); s!(:tfac, dtp / Δ0); r(:reset); s!(:timen, t0 + dtp)
            TitanCore.predict!(; halve = false, dtscale = γ1)
            s!(:thet, 0.5); s!(:dtime, γ * dtp); s!(:timen, t0 + γ * dtp)
            tw = time(); okp = newton!()
            @printf("NWPROBE step dt %.4e (TR stage %.4e): converged %s, iterations %d, final dmax %.3e, %.1f s\n", dtp, γ * dtp, okp, Int(g(:iter)), g(:dmax), time() - tw)
        end
        NWTRACE[] === nothing || close(NWTRACE[])
        flush(stdout); exit(0)
    end
    STAGECONS_HOOK[] === nothing || alt || STAGECONS_HOOK[](:start, g(:dtime), Cn)
    while true
        Δ = g(:dtime); TD.FULL_DT[] = Δ; err = NaN; Dnew = nothing; NQ[] = 1.0
        if alt
            ok, err, info, Dnew = attempt_alt!(Δ, t0, Cn)
        else
        s!(:thet, 0.5); s!(:dtime, γ * Δ); s!(:timen, t0 + γ * Δ)
        ok = newton!()
        ok && !finite_new() && (ok = false; ESTATS[:nonfinite] += 1)        # a non-finite stage is a failed attempt
        if !ok
            STATS[:fail_tr] += 1
        else
            TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](:tr, γ * Δ)
            STAGECONS_HOOK[] === nothing || STAGECONS_HOOK[](:tr, Δ, Cn)
            Cg = combos_new()
            if bdf_old!()
                s!(:thet, 1.0); s!(:dtime, C2 * Δ); s!(:timen, t0 + Δ)
                TitanCore.predict!(; halve = false)
                ok = newton!()
                ok && !finite_new() && (ok = false; ESTATS[:nonfinite] += 1)
                if ok
                    TD.STAGE_HOOK[] === nothing || TD.STAGE_HOOK[](:bdf, C2 * Δ)
                    STAGECONS_HOOK[] === nothing || STAGECONS_HOOK[](:bdf, Δ, Cn)
                    C1 = combos_new(); Cs = combos_old()
                    Dnew = Dict(s_ => (C1[s_] .- Cs[s_]) ./ (dW * Δ) for s_ in keys(C1))
                    if DSAL[] === nothing
                        ESTATS[:nofsal] += 1
                    else
                        e_, v_, z_, w_ = try estimate!(Cn, Cg, C1, Cs, Δ) catch ex
                            ex isa InterruptException && rethrow()
                            (ESTATS[:nonfinite_est] += 1) <= 3 && println("DTERR: estimate failed: ", first(sprint(showerror, ex, catch_backtrace()), 1500))
                            (Inf, "estimate failed", 0, Inf)
                        end
                        isfinite(e_) || (e_ = Inf; ESTATS[:nonfinite_est] += 1)  # a non-finite estimate rejects the step
                        err = e_; info = (v_, z_, w_)
                    end
                else
                    STATS[:fail_bdf] += 1
                end
            else
                ok = false; STATS[:fail_bdf] += 1
            end
        end
        end   # (alt)
        restore!(); s!(:thet, θ0); s!(:dtime, Δ); s!(:timen, t0 + Δ)
        ok && (SDFIXDT > 0 || !(err > 1)) && break
        if !ok && o.be                    # (testing) the same Δ with backward Euler, accepted without an estimate, not sticky
            r(:reset); s!(:thet, 1.0); s!(:timen, t0 + Δ); TitanCore.predict!(; halve = false)
            okb = newton!() && finite_new()
            okb && TD.STAGE_HOOK[] !== nothing && TD.STAGE_HOOK[](:be, Δ)
            s!(:thet, θ0)
            if okb
                ESTATS[:be] += 1; err = NaN; C1 = combos_new(); Dnew = Dict(s_ => (C1[s_] .- Cn[s_]) ./ Δ for s_ in keys(C1)); break
            end
        end
        if ok
            ESTATS[:rej_err] += 1; fac = max(o.facmin, o.safety * err^(-1 / qe))
        else
            ESTATS[:fail_newton] += 1; fac = min(o.fail, max(0.1, NQ[])); TD.msg(1)
        end
        nrej += 1; rejected = true
        nr = ESTATS[:rej_err] + ESTATS[:fail_newton]
        (nr <= 200 || nr % 100 == 0) && @printf("DTERR reject jstep %d: dt %.3e %s -> dt %.3e%s\n", Int(g(:jstep)), Δ,
            ok ? @sprintf("err %.3g (max %.3g: %s zone %d; resid %.2g, 1st filter %.3g)", err, info[3], info[1], info[2], ERR0[], ERR1[]) : "Newton failure", fac * Δ, etrace(:rej, Δ, err, info))
        s!(:tfac, fac); s!(:dtime, fac * Δ)
        s!(:jtry, g(:jtry) + 1); g(:jtry) > g(:ntry) && TD.msg(2)
        r(:reset); s!(:timen, t0 + g(:dtime)); TitanCore.predict!(; halve = false, dtscale = γ1)
    end
    Δ = g(:dtime); LASTACC[] = Δ
    DSAL[] = Dnew
    if isnan(err)
        hn = Δ
    else
        e = max(err, 1e-10)
        fac = o.pi ? o.safety * e^(-0.7 / qe) * (ERRPREV[] > 0 ? ERRPREV[]^(0.4 / qe) : 1.0) : o.safety * e^(-1 / qe)
        fac = clamp(fac, o.facmin, rejected ? 1.0 : o.facmax)
        ERRPREV[] = e; hn = fac * Δ
        LIMHIST[info[1]] = get(LIMHIST, info[1], 0) + 1
    end
    LASTERR[] = (err, info[1], info[2])
    NQ[] < 1 && (hn = min(hn, NQ[] * Δ))    # Newton difficulty of the accepted attempt (NQ): no growth beyond NQ·Δ
    hn = min(hn, o.hmax)                  # (DTERRMAX, a test bound only: default none)
    for f_ in DT_BOUNDS                   # hard bounds (mesh motion): see DT_BOUNDS
        hb = f_(); hb < hn && (hn = hb; ESTATS[:bound] = get(ESTATS, :bound, 0) + 1)
    end
    SDFIXDT > 0 && (hn = SDFIXDT)
    HNEXT[] = hn
    STATS[:steps] += 1; ESTATS[:steps] += 1
    r(:totnrg)
    r(:timstp)                            # (TITAN's rule and caps: diagnostic only)
    HEUR[] = g(:dtime); s!(:dtime, hn); s!(:tfac, hn / Δ)
    SDON && sd_state_print(Δ)
    if o.trace > 0 && ESTATS[:steps] % o.trace == 0
        @printf("DTERR jstep %d t %.6e yr: dt %.3e err %.3g (max %.3g: %s zone %d; resid %.2g, 1st filter %.3g) next %.3e | heuristic %.3e | steps %d rej %d fail %d%s\n",
            Int(g(:jstep)), g(:timen) / 3.156e7, Δ, err, info[3], info[1], info[2], ERR0[], ERR1[], hn, HEUR[], ESTATS[:steps], ESTATS[:rej_err], ESTATS[:fail_newton], etrace(:acc, Δ, err, info))
    end
    nothing
end

end # module
