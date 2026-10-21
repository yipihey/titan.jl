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

import ..TitanDriver, ..TitanCore
const TD = TitanDriver; const TF = TD.TF
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
save!() = (for s in OLD; copyto!(get!(() -> similar(vf(s)), SAVE, s), vf(s)); end)
restore!() = (for s in OLD; copyto!(vf(s), SAVE[s]); end)

function newton!()
    s!(:iter, 1)
    try
        while true
            r(:matgen); r(:matslv); r(:maxdel); r(:update)
            g(:dmax) <= g(:conv) && break
            s!(:iter, g(:iter) + 1)
            g(:iter) > g(:niter) && return false
        end
        s!(:lboos, 1); r(:matgen); r(:matslv); r(:maxdel); r(:update); s!(:lboos, 0)
    catch e
        e isa InterruptException && rethrow()
        s!(:lboos, 0); return false
    end
    all(isfinite, vf(:dn)) && all(isfinite, vf(:tn))
end

"Old arrays := y* = A·y_γ − B·y_n on conserved combinations (y_n in SAVE, y_γ in the new arrays). false if m* ≤ 0."
function bdf_old!()
    a, b = Int(g(:ngrs)), Int(g(:ngre)); lrad = Int(g(:lrad)) == 1
    dn, dvn, egn, un, xmen = vf(:dn), vf(:dvoln), vf(:egn), vf(:un), vf(:xmen)
    lag = Int(g(:lgrid)) == 3
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
    end
    @inbounds for k in a+1:b          # momentum: u·(m₋ + m₊) (boundary velocities are fixed)
        uo[k] = cb(un[k] * (mγ(k-1) + mγ(k)), u0[k] * (mn(k-1) + mn(k))) / (ms[k-1] + ms[k])
    end
    if TitanCore.FREE_OUTER[]          # dynamic outer interface: momentum of the outer half zone
        k = b + 1
        uo[k] = cb(un[k] * mγ(k-1), u0[k] * mn(k-1)) / ms[k-1]
    end
    true
end

"Drop-in for TitanDriver.step! (TR-BDF2)."
function step!()
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

end # module
