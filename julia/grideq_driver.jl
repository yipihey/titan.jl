# Included by popiii_m2c.jl (top level, after envmenc/gm_tau): the GRIDEQ moving mesh. Uses the driver's globals a, b, nsh.
include(joinpath(@__DIR__, "GridEq.jl"))
"""
GRIDEQ=mode (with ENVEUL=0,… (regular centre), lgrid = 2; default off): one global, scale-free equidistribution rule for the whole
mesh (GridEq.jl): every zone between the first interface above the centre and the outer edge (both fixed) covers the same
∫ M d ln r, M = sqrt(1 + Σ_q w_q (∂ln q/∂ln r)²) over q = ρ, P, T, E_r, m and τ (clamped to [τlo, τhi]), with a smooth gradient cap
gmax and Dorfi–Drury smoothing α (neighbours within (1 + α)/α). No markers, thresholds or zone budgets.

  mode 1 (preferred): TITAN's implicit grid equation over the whole mesh (TitanResidFast.GRIDEQ), time-smoothed per zone with
          γ_j = κ t_cross,j / Δt, t_cross = r_c/(|u| + c_s) the zone's own signal-crossing time; the τ term is frozen per step;
          Δt limited so that no interface moves more than fmove of its smaller neighbouring zone per step.
  mode 0 (fallback): the explicit target of the same rule (GridEq.target); every interface relaxes to it in max(3 Δt, κ t_cross),
          at most fmove of the zone it moves into per step and |ṙ| ≤ fv (|u| + c_s); no zone shrinks below half its narrower
          neighbour. Mesh motion never limits Δt in this mode.

GEQN=N (with ICEUL): the IC is built on N zones (even; default: envlr's count). GEQW=wρ,wP,wT,wE,wm,wτ,wu (default 1,1,1,1,1,1,0; wu: the term wu ((r ∂u/∂r)/(|u| + c_s))², mode 1 only); GEQC=gmax,alph,κ,fv,fmove,τlo,τhi,ctol,κf (default 100,1.5,30,0.5,0.5,0.1,1e4,0.3,0).

Feature tracking (κf > 0): the features of M move with the velocity v = −(∂M/∂t)(∂M/∂ln r)/((∂M/∂ln r)² + (2M)²) (in ln r per s;
the pattern speed of ∂M/∂t + v ∂M/∂ln r = 0, damped where M is flat: scale-free), estimated from M of this step and of the last
one; each zone then relaxes in max(min(κ t_cross, κf Δln r/|v|), Δln r t_cross/fv) (|v| the largest of the zone and its
neighbours; the floor keeps the grid speed below about fv (|u| + c_s)), so that the mesh keeps up with a moving photosphere or
shock without being fast everywhere.
With ICEUL the IC is built on the rule's mesh of the IC. No state beyond the mesh: a checkpoint restarts with the same setting.

GEQX=wx[,xfl] (default 0: off): a composition term wx Σ_s (∂ ln(x_s + xfl)/∂ ln r)² over every species of the network (abundances
per H nucleus, ChemNet.full; xfl default 1e-3), frozen per step like the τ term (mode 1) and in the explicit target (mode 0). A
species far below xfl contributes nothing, so only changes of composition that carry a share ≳ xfl of the nuclei (ionization and
H2 dissociation fronts) are resolved; no list of species, no thresholds per front.
GEQRAMP=κ0,n[,f] (default off): κ of the time smoothing (and with it the speed floor) rises geometrically from κ0 to GEQC's κ over
the first n steps after the start of the process, if the mesh is then far from the rule's, max_j(M_j Δln r_j) > f × median (default
f 10; a fixed-mesh checkpoint with its front in one zone: 26.5 for zdiagB's S5; an equidistributed mesh: ~3.4). A restart from a
fixed-mesh checkpoint needs it: κ 30 would take ~30 crossing times of the front to adapt the mesh, while Δt stays at the fixed
mesh's and the fixed-mesh run fails within ~30 steps. The next segment of the same run (mesh already adapted) does not ramp again.
"""
const GRIDEQP = haskey(ENV, "GRIDEQ") ? parse(Int, ENV["GRIDEQ"]) : nothing
geqv(key, defs) = (v = haskey(ENV, key) ? parse.(Float64, split(ENV[key], ",")) : Float64[]; [i <= length(v) ? v[i] : defs[i] for i in eachindex(defs)])
const GEQW = geqv("GEQW", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
const GEQC = geqv("GEQC", [100.0, 1.5, 30.0, 0.5, 0.5, 0.1, 1e4, 0.3, 0.0])
const GEQX = geqv("GEQX", [0.0, 1e-3])
const GEQRAMP = haskey(ENV, "GEQRAMP") ? parse.(Float64, split(ENV["GEQRAMP"], ",")) : nothing
const GEQJ0 = Ref(-1)                  # jstep at the first grideq_step! of this process (GEQRAMP; -2: no ramp, the mesh is close to the rule)
"max_j(M_j Δln r_j) / median over the zones above the centre: how far the mesh is from equidistribution (GEQRAMP, GEQREMAP)."
function geq_quality(s)
    D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.X, wx = GEQX[1])
    q = [D.M[i] * log(s.r[i+2] / s.r[i+1]) for i in eachindex(D.M)]; qs = sort(q)
    (; ratio = maximum(q) / qs[(length(qs) + 1) ÷ 2], k = argmax(q) + 1)
end
"""
GEQREMAP=f (with GRIDEQ and a restart from an lgrid = 2 checkpoint, regular centre; default 0 = off): if the checkpoint's mesh is far
from the rule's, max_j(M_j Δln r_j) > f × median (a fixed-mesh checkpoint whose front sits in one zone), the state is rebuilt once on
GRIDEQ's explicit target mesh of the same zone count (envbuild!: enclosed mass kept at every source interface with ENVCONS=1; T, E_r,
x interpolated in ln r between zone centres; u, L at the interfaces; lagged material state reset as at the hand-over), mass zero point
and fixed-mesh flags set again as ICEUL. A mesh already close to the rule (the next segment of the same run) is not remapped.
"""
const GEQREMAP = parse(Float64, get(ENV, "GEQREMAP", "0"))
const GEQPREV = Ref{Any}(nothing)      # (t, x, M) of the last step: the monitor function's own motion (feature tracking)
const GEQLAST = Ref{Any}(nothing)
"State on the present mesh (interfaces centre first; m per steradian)."
function geq_state()
    r = TF.var(:rn)[a:b+1]; u = TF.var(:un)[a:b+1]; m = envmenc()
    (; r, u, m, ρ = TF.var(:dn)[a:b], P = TF.var(:pgn)[a:b], T = TF.var(:tn)[a:b], E = TF.var(:ern)[a:b], cs = TF.var(:asn)[a:b], τ = gm_tau(a),
       X = GEQX[1] > 0 ? geq_lx() : nothing)
end
"GEQX: ln(x_s + xfl) of every network species (rows) in the zones a…b (columns)."
function geq_lx()
    mat isa PM.NetGas || return nothing
    d = [ChemNet.full(mat.x[k], mat.net) for k in a:b]; sp = sort(collect(keys(d[1])))
    Float64[log(max(Float64(ForwardDiff.value(d[j][q])), 0.0) + GEQX[2]) for q in sp, j in eachindex(d)]
end
geq_target(s; N = nsh - 1) = GridEq.target(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ, N; w = Tuple(GEQW), gmax = GEQC[1], alph = GEQC[2], τlo = GEQC[6], τhi = GEQC[7],
                                            X = s.X, wx = GEQX[1])
"GEQREMAP: rebuild a restarted fixed-mesh state on the rule's mesh (see GEQREMAP)."
function geq_remap!()
    ENVCENTRE || error("GEQREMAP needs a regular centre (ENVEUL=0,...)")
    s = geq_state(); Q = geq_quality(s)
    pr("GEQREMAP: max M dln r / median %.3g (threshold %.3g) at zone %d (r %.4g R_sun)\n", Q.ratio, GEQREMAP, Q.k, s.r[Q.k] / Rsun)
    Q.ratio > GEQREMAP || return nothing
    E0 = geq_energy(); rq = geq_target(s)
    envbuild!(; rIset = vcat(0.0, rq))
    get(ENV, "ICXMZ", "1") == "1" && xmzero!()
    envflags!(); ENVMS0[] = 0.0
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    E1 = geq_energy()
    pr("GEQREMAP: rebuilt on the rule's mesh (%d zones, dln r %.4f .. %.4f); mass %.10e -> %.10e Msun; E_int %.6e -> %.6e, E_kin %.4e -> %.4e, E_rad %.4e -> %.4e, E_grav %.6e -> %.6e erg (total change %.3e of |E_grav|)\n",
       nsh, minimum(diff(log.(rq))), maximum(diff(log.(rq))), E0.M / Msun, E1.M / Msun, E0.ei, E1.ei, E0.ek, E1.ek, E0.er, E1.er, E0.eg, E1.eg,
       ((E1.ei + E1.ek + E1.er + E1.eg) - (E0.ei + E0.ek + E0.er + E0.eg)) / abs(E0.eg))
    nothing
end
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
    (; nsk, rsk = s.r[ku], uk = s.u[ku], nph, nph5, dT, dlmin = minimum(dl), rat)
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
    tfe = fill(Inf, n); nfast = 0                                 # feature-crossing time of each zone (feature tracking)
    if κf > 0
        D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.X, wx = GEQX[1])
        t = TF.getscalar(:timen); pv = GEQPREV[]
        if pv !== nothing && t > pv[1]
            Mp = [GridEq.lin(pv[2], pv[3], x) for x in D.x]; nd = length(D.x); v = zeros(nd)
            for i in 1:nd
                il, ih = max(i - 1, 1), min(i + 1, nd); Md = (D.M[ih] - D.M[il]) / (D.x[ih] - D.x[il])
                v[i] = -(D.M[i] - Mp[i]) / (t - pv[1]) * Md / (Md^2 + (2D.M[i])^2)
            end
            for i in 1:nd
                vm = maximum(abs, v[max(i - 1, 1):min(i + 1, nd)]); j = i + 1      # (D covers zones 2…n)
                vm > 0 && (tfe[j] = log(s.r[j+1] / s.r[j]) / vm)
            end
        end
        GEQPREV[] = (t, D.x, D.M)
    end
    if GRIDEQP == 1
        tauv = zeros(length(TF.var(:rn))); sx = zeros(length(TF.var(:rn)))
        lτ = log.(clamp.(max.(s.τ, 1e-300), GEQC[6], GEQC[7])); lτc = [0.5 * (lτ[j] + lτ[j+1]) for j in 1:n]
        for j in 1:n
            # (feature tracking; floor: a zone moving its own width in t_relax goes at most fv (|u| + c_s))
            dxj = j >= 2 ? log(s.r[j+1] / s.r[j]) : 1.0
            # (κf = 0: no feature tracking; 0 × Inf would be NaN and min() propagates it into every grid row)
            tr_ = max(κf > 0 ? min(κ * tc[j], κf * tfe[j]) : κ * tc[j], dxj * tc[j] / GEQC[4] * κ / GEQC[3]); tr_ < κ * tc[j] && (nfast += 1); tauv[a+j-1] = tr_   # (GEQRAMP: floor ramped with κ)
            jm, jp = max(j - 1, 2), min(j + 1, n); h = jp - jm == 2 ? 0.5 : 1.0
            sx[a+j-1] = GEQW[6] * (h * (lτc[jp] - lτc[jm]))^2
            # velocity jump in units of the signal speed, (r ∂u/∂r)/(|u| + c_s): resolves shocks (the zone's own Δu, per zone index)
            GEQW[7] > 0 && (sx[a+j-1] += GEQW[7] * ((s.u[j+1] - s.u[j]) / (0.5 * (abs(s.u[j]) + abs(s.u[j+1])) + s.cs[j]))^2)
            # GEQX: composition, ln(x + xfl) of every species, centred zone differences as the τ term
            s.X === nothing || (sx[a+j-1] += GEQX[1] * sum((h * (s.X[q, jp] - s.X[q, jm]))^2 for q in axes(s.X, 1)))
        end
        TX.GRIDEQ[] = (; k1 = a + 1, k2 = b + 1, w = (GEQW[1], GEQW[2], GEQW[3], GEQW[4], GEQW[5]), alph = GEQC[2], gmax = GEQC[1], tauv, sx)
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
    @sprintf("        grid eq: shock at %.4g R_sun (u %.1f km/s) across %d zones, photosphere (0.1 <= tau <= 10) across %d interfaces (1 <= tau <= 1e5: %d), largest |dln T| per interface %.3f; dln r min %.3e, max neighbour ratio %.2f; grid speed max %.3f (|u| + c_s); %d zones relaxing at their feature speed; kappa %.3g\n",
             d.rsk / Rsun, d.uk / 1e5, d.nsk, d.nph, d.nph5, d.dT, d.dlmin, d.rat, vmax, d.nfast, d.κ)
end
