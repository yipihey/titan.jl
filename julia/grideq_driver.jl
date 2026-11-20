# Included by popiii_m2c.jl (top level, after envmenc/gm_tau): the GRIDEQ moving mesh. Uses the driver's globals a, b, nsh.
include(joinpath(@__DIR__, "GridEq.jl"))
"""
GRIDEQ=mode (with ENVEUL=0,… (regular centre), lgrid = 2; default off): one global, scale-free equidistribution rule for the whole
mesh (GridEq.jl): every zone between the first interface above the centre and the outer edge (both fixed in radius; with ICMLOG the first interface
is instead Lagrangian, its grid row holds its enclosed mass: TitanResidFast.GEQ_LAGFIRST, and it counts in the fmove limit) covers the same
∫ M d ln r, M = sqrt(1 + Σ_q w_q (∂ln q/∂ln r)²) over q = ρ, P, T, E_r, m and τ (clamped to [τlo, τhi]), with a smooth gradient cap
gmax and Dorfi–Drury smoothing α (neighbours within (1 + α)/α). No markers, thresholds or zone budgets.

  mode 1 (preferred): TITAN's implicit grid equation over the whole mesh (TitanResidFast.GRIDEQ), time-smoothed per zone with
          γ_j = κ t_cross,j / Δt, t_cross = r_c/(|u| + c_s) the zone's own signal-crossing time; the τ term (GEQW entry 6, wτ) is formed from the Newton iterate (τ = Σχρ Δr from the outer zone, ln τ smoothly confined to [ln τlo, ln τhi]) and differentiated; the speed floor of the relaxation time too;
          Δt limited so that no interface moves more than fmove of its smaller neighbouring zone per step.
  mode 0 (fallback): the explicit target of the same rule (GridEq.target); every interface relaxes to it in max(3 Δt, κ t_cross),
          at most fmove of the zone it moves into per step and |ṙ| ≤ fv (|u| + c_s); no zone shrinks below half its narrower
          neighbour. Mesh motion never limits Δt in this mode.

GEQN=N (with ICEUL): the IC is built on N zones (even; default: envlr's count). GEQW=wρ,wP,wT,wE,wm,wτ,wu,wχ (default 1,1,1,1,1,1,0,0; wτ: (Δ ln τ)² of the Newton iterate within τlo..τhi, mode 1; wχ: (∂ln(χρ)/∂ln r)² of the Newton iterate, mode 1 only; wu: the term wu ((r ∂u/∂r)/(|u| + c_s))², mode 1 only); GEQC=gmax,alph,κ,fv,fmove,τlo,τhi,ctol,κf,f_dyn (default 100,1.5,30,0.5,0.5,0.1,1e4,0.3,0,0.2; mode 1: κf is no longer used in τ; f_dyn: the grid relaxation time of every zone is smin(κ t_cross, f_dyn/√(Gρ)) of the Newton iterate, see TitanResidFast.grideq_rows!).

Feature tracking (κf > 0): the features of M move with the velocity v = −(∂M/∂t)(∂M/∂ln r)/((∂M/∂ln r)² + (2M)²) (in ln r per s;
the pattern speed of ∂M/∂t + v ∂M/∂ln r = 0, damped where M is flat: scale-free), estimated from M of this step and of the last
one; each zone then relaxes in max(min(κ t_cross, κf Δln r/|v|), Δln r t_cross/fv) (|v| the largest of the zone and its
neighbours; the floor keeps the grid speed below about fv (|u| + c_s)), so that the mesh keeps up with a moving photosphere or
shock without being fast everywhere.
With ICEUL the IC is built on the rule's mesh of the IC. No state beyond the mesh: a checkpoint restarts with the same setting.

GEQDEC=pr,pm[,C*] (default 6,6,1): resolution floor, at least pr zones per decade of r and pm per decade of enclosed mass. M² gets the constant
ω_r² and the term ω_m² (∂ln m/∂ln r)² (added to the monitor's G² before the gmax cap; mode 1 and mode 0), ω_x = p_x C*/ln 10. Equidistribution gives Δln r_j = C/M_j
with C = ∫M dln r/N; M ≥ sqrt(1 + ω_r²) > ω_r and, at fixed ∂ln m/∂ln r, Δln m_j ≤ C/ω_m, so Δlog10 r ≤ (C/C*)/pr and Δlog10 m ≤ (C/C*)/pm
(the floor is exact when C ≤ C*; C is printed in the grid-eq line as 'C'; it is a weight, not a clamp: if physics raises C above C*, the floor degrades as C/C*). The floor terms act on the shells from the first interface outward (with ICMLOG the first interface is Lagrangian and not a grid row; its enclosed mass is pinned, so ∂ln m/∂ln r of the first shells is one-sided from the pinned interface). The 6 (not 5) leaves ~25% slack for C > C* so that >= 5 zones per decade survive. The first
term of M (the 1) and wm still act on top.
GEQX=wx[,xfl] (default 0: off): a composition term wx Σ_s (∂ ln(x_s + xfl)/∂ ln r)² over every species of the network (abundances
per H nucleus, ChemNet.full; xfl default 1e-3), frozen per step like the τ term (mode 1) and in the explicit target (mode 0). A
species far below xfl contributes nothing, so only changes of composition that carry a share ≳ xfl of the nuclei (ionization and
H2 dissociation fronts) are resolved; no list of species, no thresholds per front.
GEQRAMP=κ0,n[,f[,fvr]] (default off; fvr: the ramped speed floor never allows more than fvr (|u| + c_s), default unlimited): κ of the time smoothing (and with it the speed floor) rises geometrically from κ0 to GEQC's κ over
the first n steps after the start of the process, if the mesh is then far from the rule's, max_j(M_j Δln r_j) > f × median (default
f 10; a fixed-mesh checkpoint with its front in one zone: 26.5 for zdiagB's S5; an equidistributed mesh: ~3.4). A restart from a
fixed-mesh checkpoint needs it: κ 30 would take ~30 crossing times of the front to adapt the mesh, while Δt stays at the fixed
mesh's and the fixed-mesh run fails within ~30 steps. The next segment of the same run (mesh already adapted) does not ramp again.
"""
const GRIDEQP = let v_ = get(ENV, "GRIDEQ", "off"); v_ in ("", "off") ? nothing : parse(Int, v_) end   # GRIDEQ=off: the fixed mesh (the stack default is 1)
geqv(key, defs) = (v = haskey(ENV, key) ? parse.(Float64, split(ENV[key], ",")) : Float64[]; [i <= length(v) ? v[i] : defs[i] for i in eachindex(defs)])
const GEQW = geqv("GEQW", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])   # entry 8 wχ: ln(χρ) of the Newton iterate (mode 1; with wτ = 0 and GEQX off the monitor has no lagged term)
const GEQC = geqv("GEQC", [100.0, 1.5, 30.0, 0.5, 0.5, 0.1, 1e4, 0.3, 0.0, 0.2])
const GEQX = geqv("GEQX", [0.0, 1e-3])
# GEQGP=w_g,L_ref [L_sun] (provisional, default 0 = off; implicit grid equation only): density jump x gravity power of the zone in the monitor
# (TitanResidFast.grideq_rows!, docs/POPIII_ENERGY.md §14.6)
const GEQGP = geqv("GEQGP", [0.0, 100.0])
const GEQDEC = geqv("GEQDEC", [6.0, 6.0, 1.0])   # zones per decade in r, in enclosed mass; C* (see above)
const GEQFL = ((GEQDEC[1] * GEQDEC[3] / log(10))^2, (GEQDEC[2] * GEQDEC[3] / log(10))^2)   # (ω_r², ω_m²)
const GEQRAMP = haskey(ENV, "GEQRAMP") ? parse.(Float64, split(ENV["GEQRAMP"], ",")) : nothing
const GEQJ0 = Ref(-1)                  # jstep at the first grideq_step! of this process (GEQRAMP; -2: no ramp, the mesh is close to the rule)
"max_j(M_j Δln r_j) / median over the zones above the centre: how far the mesh is from equidistribution (GEQRAMP)."
function geq_quality(s)
    D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL)
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
                                            X = s.Xd, wx = 1.0, fl = GEQFL)
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
    Dn = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL)
    Cq = sum(Dn.M .* dl) / length(s.ρ)                 # C = ∫M dln r / N of the explicit monitor (the floor holds when C <= C*)
    (; nsk, rsk = s.r[ku], uk = s.u[ku], nph, nph5, dT, dlmin = minimum(dl), rat, Cq)
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
        D = GridEq.density(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ; w = Tuple(GEQW), gmax = GEQC[1], τlo = GEQC[6], τhi = GEQC[7], X = s.Xd, wx = 1.0, fl = GEQFL)
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
        TX.GRIDEQ[] = (; k1 = a + 1, k2 = b + 1, w = (GEQW[1], GEQW[2], GEQW[3], GEQW[4], GEQW[5]), alph = GEQC[2], gmax = GEQC[1], tauv, tfl, flc = 1 / min(GEQC[4] * GEQC[3] / κ, fvr), wtau = GEQW[6], tlo = GEQC[6], thi = GEQC[7], kap = κ, fdyn = GEQC[10], sx, wχ = GEQW[8], gp = (GEQGP[1], GEQGP[2] * 3.828e33 / (4π)), fr = GEQFL[1], fm = GEQFL[2])
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
    @sprintf("        grid eq: shock at %.4g R_sun (u %.1f km/s) across %d zones, photosphere (0.1 <= tau <= 10) across %d interfaces (1 <= tau <= 1e5: %d), largest |dln T| per interface %.3f; dln r min %.3e, max neighbour ratio %.2f; grid speed max %.3f (|u| + c_s); %d zones relaxing at their feature speed; kappa %.3g; C %.3g (C* %.3g)\n",
             d.rsk / Rsun, d.uk / 1e5, d.nsk, d.nph, d.nph5, d.dT, d.dlmin, d.rat, vmax, d.nfast, d.κ, d.Cq, GEQDEC[3])
end
