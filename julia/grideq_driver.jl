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

GEQN=N (with ICEUL): the IC is built on N zones (even; default: envlr's count). GEQW=wρ,wP,wT,wE,wm,wτ (default 1,1,1,1,1,1); GEQC=gmax,alph,κ,fv,fmove,τlo,τhi,ctol (default 100,1.5,3,0.5,0.5,0.1,1e4,0.3).
With ICEUL the IC is built on the rule's mesh of the IC. No state beyond the mesh: a checkpoint restarts with the same setting.
"""
const GRIDEQP = haskey(ENV, "GRIDEQ") ? parse(Int, ENV["GRIDEQ"]) : nothing
geqv(key, defs) = (v = haskey(ENV, key) ? parse.(Float64, split(ENV[key], ",")) : Float64[]; [i <= length(v) ? v[i] : defs[i] for i in eachindex(defs)])
const GEQW = geqv("GEQW", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
const GEQC = geqv("GEQC", [100.0, 1.5, 3.0, 0.5, 0.5, 0.1, 1e4, 0.3])
const GEQLAST = Ref{Any}(nothing)
"State on the present mesh (interfaces centre first; m per steradian)."
function geq_state()
    r = TF.var(:rn)[a:b+1]; u = TF.var(:un)[a:b+1]; m = envmenc()
    (; r, u, m, ρ = TF.var(:dn)[a:b], P = TF.var(:pgn)[a:b], T = TF.var(:tn)[a:b], E = TF.var(:ern)[a:b], cs = TF.var(:asn)[a:b], τ = gm_tau(a))
end
geq_target(s; N = nsh - 1) = GridEq.target(s.r, s.m, s.ρ, s.P, s.T, s.E, s.τ, N; w = Tuple(GEQW), gmax = GEQC[1], alph = GEQC[2], τlo = GEQC[6], τhi = GEQC[7])
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
    nph = count(t -> 0.1 <= t <= 10, s.τ)
    dl = [log(s.r[j+1] / s.r[j]) for j in 2:length(s.ρ)]
    rat = maximum(max(dl[j] / dl[j-1], dl[j-1] / dl[j]) for j in 2:length(dl))
    (; nsk, rsk = s.r[ku], uk = s.u[ku], nph, dlmin = minimum(dl), rat)
end
function grideq_step!()
    s = geq_state(); n = nsh; dt = TF.getscalar(:dtime); tc = geq_tcross(s); κ = GEQC[3]
    if GRIDEQP == 1
        tauv = zeros(length(TF.var(:rn))); sx = zeros(length(TF.var(:rn)))
        lτ = log.(clamp.(max.(s.τ, 1e-300), GEQC[6], GEQC[7])); lτc = [0.5 * (lτ[j] + lτ[j+1]) for j in 1:n]
        for j in 1:n
            tauv[a+j-1] = κ * tc[j]
            jm, jp = max(j - 1, 2), min(j + 1, n); h = jp - jm == 2 ? 0.5 : 1.0
            sx[a+j-1] = GEQW[6] * (h * (lτc[jp] - lτc[jm]))^2
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
    GEQLAST[] = geq_diag(s)
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
    @sprintf("        grid eq: shock at %.4g R_sun (u %.1f km/s) across %d zones, photosphere (0.1 <= tau <= 10) across %d interfaces; dln r min %.3e, max neighbour ratio %.2f; grid speed max %.3f (|u| + c_s)\n",
             d.rsk / Rsun, d.uk / 1e5, d.nsk, d.nph, d.dlmin, d.rat, vmax)
end
