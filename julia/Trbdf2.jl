"""
    Trbdf2

Second-order, L-stable time stepping for the Julia TITAN (Eulerian grid, lgrid = 2): TR-BDF2 (Bank et al. 1985) built
from TITAN's own θ-scheme, so the residual is unchanged.

One step of size Δ from y_n:
1. **TR stage:** TITAN's θ-scheme with θ = ½ (implicit midpoint) over γΔ, γ = 2 − √2, giving y_γ.
2. **BDF2 stage:** y_{n+1} − [A y_γ − B y_n] = (1−γ)/(2−γ) Δ · f(y_{n+1}), A = 1/(γ(2−γ)), B = (1−γ)²/(γ(2−γ)).
   That is a backward-Euler step (θ = 1) of size (1−γ)/(2−γ)Δ from the "old" state y* = A y_γ − B y_n. TITAN's
   residual reads the old state only through conserved combinations (mass m = ρ dV, m e_g, m(e_g + E/ρ), E dV,
   mass-weighted momentum, the flux-row inertia F·2M/(ρ₋+ρ₊), the exterior mass coordinate), so y* is formed on those
   and written back into the old arrays (do, ego, egro, ero, uo, fro, xmeo). Terms with weight 1 − θ vanish at θ = 1.
After the step the old arrays are restored to y_n, so TITAN's step-size control (timstp) and the energy ledger see the
whole step Δ. A failed Newton iteration in either stage restarts the step from y_n with Δ × TFAC[].

Limits: lgrid = 2 only (on a moving grid the volumes in y* would also have to be combined); the energy ledger (totnrg)
books the step as one backward-Euler step, so it is approximate.
Use: `TitanDriver.STEPFN[] = Trbdf2.step!`.
"""
module Trbdf2

import ..TitanDriver, ..TitanCore
const TD = TitanDriver; const TF = TD.TF
g(s) = TD.g(s); s!(s, x) = TD.s!(s, x); r(s) = TD.r(s)
vf(s) = TF.var(s)::Vector{Float64}

const γ = 2 - sqrt(2)
const A = 1 / (γ * (2 - γ)); const B = (1 - γ)^2 / (γ * (2 - γ)); const C2 = (1 - γ) / (2 - γ)
const TFAC = Ref(0.05)
const STATS = Dict(:steps => 0, :fail_tr => 0, :fail_bdf => 0)
const OLD = (:do, :ego, :egro, :ero, :fro, :uo, :xmeo)
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
    d0, dv0, eg0, u0, xm0 = SAVE[:do], vf(:dvolo), SAVE[:ego], SAVE[:uo], SAVE[:xmeo]
    do_, ego, uo, xmeo = vf(:do), vf(:ego), vf(:uo), vf(:xmeo)
    rng = a-1:b+2
    mγ(k) = dn[k] * dvn[k]; mn(k) = d0[k] * dv0[k]
    ms = zeros(length(dn))
    @inbounds for k in rng
        ms[k] = A * mγ(k) - B * mn(k)
        ms[k] > 0 || return false
    end
    @inbounds for k in rng
        ego[k] = (A * mγ(k) * egn[k] - B * mn(k) * eg0[k]) / ms[k]
        do_[k] = ms[k] / dv0[k]
        xmeo[k] = A * xmen[k] - B * xm0[k]
    end
    if lrad
        egrn, ern, frn = vf(:egrn), vf(:ern), vf(:frn)
        egr0, er0, fr0 = SAVE[:egro], SAVE[:ero], SAVE[:fro]
        egro, ero, fro = vf(:egro), vf(:ero), vf(:fro)
        @inbounds for k in rng
            egro[k] = (A * mγ(k) * egrn[k] - B * mn(k) * egr0[k]) / ms[k]
            ero[k] = (A * ern[k] * dvn[k] - B * er0[k] * dv0[k]) / dv0[k]
        end
        @inbounds for k in a:b+2      # flux-row inertia Q = F·2M/(ρ₋ + ρ₊), M = m₋ + m₊
            Qγ = frn[k] * 2 * (mγ(k-1) + mγ(k)) / (dn[k-1] + dn[k])
            Qn = fr0[k] * 2 * (mn(k-1) + mn(k)) / (d0[k-1] + d0[k])
            fro[k] = (A * Qγ - B * Qn) * (do_[k-1] + do_[k]) / (2 * (ms[k-1] + ms[k]))
        end
    end
    @inbounds for k in a:b+2          # momentum: u·(m₋ + m₊)
        uo[k] = (A * un[k] * (mγ(k-1) + mγ(k)) - B * u0[k] * (mn(k-1) + mn(k))) / (ms[k-1] + ms[k])
    end
    true
end

"Drop-in for TitanDriver.step! (TR-BDF2)."
function step!()
    Int(g(:lgrid)) == 2 || error("Trbdf2: Eulerian grid (lgrid = 2) only")
    r(:peq)
    θ0 = g(:thet); t0 = g(:timeo)
    save!()
    s!(:jtry, 1); s!(:jback, 0)
    while true
        Δ = g(:dtime)
        s!(:thet, 0.5); s!(:dtime, γ * Δ); s!(:timen, t0 + γ * Δ)
        ok = newton!()
        if !ok
            STATS[:fail_tr] += 1
        else
            if bdf_old!()
                s!(:thet, 1.0); s!(:dtime, C2 * Δ); s!(:timen, t0 + Δ)
                ok = newton!()
                ok || (STATS[:fail_bdf] += 1)
            else
                ok = false; STATS[:fail_bdf] += 1
            end
        end
        restore!(); s!(:thet, θ0); s!(:dtime, Δ); s!(:timen, t0 + Δ)
        ok && break
        TD.msg(1)
        s!(:tfac, TFAC[]); s!(:dtime, TFAC[] * Δ)
        s!(:jtry, g(:jtry) + 1); g(:jtry) > g(:ntry) && TD.msg(2)
        r(:reset); s!(:timen, t0 + g(:dtime))
    end
    STATS[:steps] += 1
    r(:totnrg)
    r(:timstp)
    nothing
end

end # module
