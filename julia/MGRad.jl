# ------------------------------------------------------------------ MGROUPS: implicit multigroup radiation (prototype)
# Included into TitanResidFast. Design and tests: docs/POPIII_MULTIGROUP.md. Everything here is behind MGROUPS (default off: MG[] is
# false, no code path below runs and no existing path changes; bit-identical off).
"""
MGROUPS (popiii_m2c `MGROUPS=n` (default edge set) or `MGROUPS=e1,e2,…` (interior group edges in eV, ascending); default off):
the radiation field as G groups, each with its own energy density E_g and flux F_g, all Newton unknowns of the hydro system
(E_g at MGI0 + 2g − 1, F_g at MGI0 + 2g; MGI0 = 9 with COUPLED, so neqn = 9 + 2G). The grey rows of E_r, F_r become the sums
E_r = Σ_g E_g, F_r = Σ_g F_g (algebraic rows), so every term of the hydro that reads E_r, F_r and the grey Eddington factors
(the total-energy and gas+radiation rows, the pressure f E_r, the formal solution, Δt control) sees the total field unchanged, and:

* per group the moment rows of TITAN's grey system (TitanResidFast RADIATION ENERGY / MOMENTUM, comoving frame, Mihalas & Mihalas
  1984; Gehmeyr & Mihalas 1994) with the group's exchange x_cap (c κ_E,g E_g − 4π κ_P,g B_g(T)) ρ ΔV, B_g = b_g(T) σT⁴/π, and its
  interface opacity χ̄_g (harmonic, or SLABFLUX's χ_n of the group's χ_g, SUBZONE on top), the FLUXTIE closure with the grey h
  (SLABTIE: the donor h at the group's slab depth), the outer boundary F_g = c g E_g,s − π(2g+1) I_cmb b_g(T_cmb);
* the radiation force on the gas is Σ_g χ̄_g F_g/c (the grey rows' χ̄ F/c is replaced; the energy rows need nothing: their work and
  flux terms read E_r = Σ E_g and F_r = Σ F_g, with the grey f common to all groups in this version);
* group opacities: κ_X,g = κ_X(iterate) × r_X,g, the grey means of the material at the Newton iterate (exact AD) times lagged group
  ratios r_X,g (per step, PopIIIMaterial.mg_ratios!: G2 continuum per-group means over the code's own grey means of the same call;
  1 in table zones, the diffusion-limit split), plus, in the ionizing group of lit network zones, the absolute H(1s) bound–free
  absorption MGK1 (lagged; case B: no ground-state recombination emission, the network's own convention);
* the chemistry's lagged energy moves (RAYQ, GASQ) come out of group MGQG (optical).
Left for later stages (docs): per-group formal solution (f_g, g_g, h_g), Doppler shifts between groups, exact group opacities at the
iterate, per-process star-light factors at the iterate.
"""
const MG = Ref(false)
const MGN = Ref(0)
const MGEDGE = Ref(Float64[])                   # interior edges [eV], ascending
const MGI0 = Ref(9)
const MGEVK = Ref(1.602176634e-12 / 1.380649e-16)   # eV / k_B [K]; set by the driver from ChemNet's constants
const MGQG = Ref(2)                             # the group that carries RAYQ/GASQ (optical)
const MGBFLOOR = 1e-200                         # floor of a group's Planck fraction (energy 1e-200 aT⁴: none; keeps ln E_g finite)
const MGCLAMP = Ref(3.0)                        # clamp of a correction of ln E_g per Newton iteration
const MGEFL = Ref(1e-3)                         # convergence / error weight floor of a group: changes below MGEFL·E_r do not count
const MGTCMB = Ref(0.0)
vvec(n, m) = [zeros(m) for _ in 1:n]
const MGEN = Ref(Vector{Vector{Float64}}()); const MGFN = Ref(Vector{Vector{Float64}}())     # new level (Newton iterate base)
const MGEO = Ref(Vector{Vector{Float64}}()); const MGFO = Ref(Vector{Vector{Float64}}())     # old level (start of the stage)
const MGESO = Ref(Vector{Vector{Float64}}()); const MGFSO = Ref(Vector{Vector{Float64}}())   # limiter slopes of the old E_g/ρ and face F_g
const MGFNOM = Ref(Vector{Vector{Float64}}())   # nominal F_g scale per interface (correction F_g += δ·nom)
const MGRE = Ref(Vector{Vector{Float64}}()); const MGRP = Ref(Vector{Vector{Float64}}()); const MGRR = Ref(Vector{Vector{Float64}}())
const MGK1 = Ref(Vector{Vector{Float64}}())     # absolute H(1s) absorption [cm² g⁻¹] per group (ionizing groups of lit network zones)
"(diagnostic) of the last residual: Σ_g E_g/E_r − 1 worst, the group luminosities through each interface are in MGEN/MGFN"
const MGNEWTON = Ref(0)
mgiE(g) = MGI0[] + 2g - 1; mgiF(g) = MGI0[] + 2g

"Allocate the group arrays for mgr zones (values kept where present)."
function mg_alloc!(mgr)
    G = MGN[]
    for R_ in (MGEN, MGFN, MGEO, MGFO, MGESO, MGFSO, MGFNOM, MGK1)
        (length(R_[]) == G && all(length(a) == mgr for a in R_[])) || (R_[] = vvec(G, mgr))
    end
    for R_ in (MGRE, MGRP, MGRR)
        (length(R_[]) == G && all(length(a) == mgr for a in R_[])) || (R_[] = [ones(mgr) for _ in 1:G])
    end
    nothing
end

# ---------------------------------------------------------------- Planck group fractions (generic in the number type)
"∫₀ˣ t³/(eᵗ − 1) dt for x < 1: Bernoulli series B₀…B₁₀ (as RadGroups of branch radgroups-3g; next term ≲ 1e-10 relative at x = 1)."
@inline function mg_plow(x)
    x2 = x * x
    x^3 * (1 / 3 - x / 8 + x2 / 60 - x2^2 / 5040 + x2^3 / 272160 - x2^4 / 13305600 + x2^5 * (5 / 66) / (13 * 3628800))
end
"∫ₓ^∞ t³/(eᵗ − 1) dt for x ≥ 1: Σₙ e^{−nx}(x³/n + 3x²/n² + 6x/n³ + 6/n⁴) (no cancellation for large x)."
@inline function mg_phigh(x)
    s = zero(x); xv = val(x)
    for n in 1:60
        e = exp(-n * x)
        s += e * (x^3 / n + 3x^2 / n^2 + 6x / n^3 + 6 / n^4)
        exp(-n * xv) < 1e-18 && break
    end
    s
end
"Fraction of σT⁴/π above x = hν/kT."
@inline mg_above(x) = val(x) >= 1 ? mg_phigh(x) * (15 / π^4) : 1 - mg_plow(x) * (15 / π^4)
"Fraction of σT⁴/π below x."
@inline mg_below(x) = val(x) >= 1 ? 1 - mg_phigh(x) * (15 / π^4) : mg_plow(x) * (15 / π^4)
"""
Planck fraction b_g(T) of group g (edges MGEDGE, eV), floored at MGBFLOOR. Interior groups whose lower edge has x ≥ 1 are the difference
of two tails (no cancellation); else of two lower integrals.
"""
@inline function mg_bfrac(g, T, edges = MGEDGE[])
    G = length(edges) + 1; s = MGEVK[] / T
    b = if g == 1
        mg_below(edges[1] * s)
    elseif g == G
        mg_above(edges[G-1] * s)
    else
        xl = edges[g-1] * s; xh = edges[g] * s
        val(xl) >= 1 ? mg_phigh(xl) * (15 / π^4) - mg_phigh(xh) * (15 / π^4) : mg_below(xh) - mg_below(xl)
    end
    val(b) < MGBFLOOR ? MGBFLOOR + zero(b) : b
end

# ---------------------------------------------------------------- Newton hooks (chained with the COUPLED ones by popiii_m2c)
"Old-level limiter slopes of E_g/ρ (zones) and of the face value F_g·2/(ρ₋+ρ₊), as the species slopes (cpl_newstep!)."
function mg_slopes!()
    ngrs = gi(:ngrs); ngre = gi(:ngre); cadv = gf(:cadv); eps = gf(:epsadv); dO = v(:do)
    q = get!(() -> zeros(TF.MGR), BUFS, :mgq)::Vector{Float64}
    for g in 1:MGN[]
        for (kind, src, sl) in ((1, MGEO[][g], MGESO[][g]), (2, MGFO[][g], MGFSO[][g]))
            fill!(sl, 0.0)
            @inbounds for k in ngrs-2:ngre+2
                q[k] = kind == 1 ? src[k] / dO[k] : (k > 1 ? src[k] * 2 / (dO[k] + dO[k-1]) : 0.0)
            end
            @inbounds for k in ngrs:ngre+1
                a = q[k] - q[k-1]; b = q[k+1] - q[k]; sdq = a + b; den = sdq + (sdq >= 0 ? eps : -eps)
                sl[k] = a * b >= 0 ? cadv * a * b / den : 0.0
            end
        end
    end
    nothing
end
"Start of a step: old := new, slopes, nominal flux scales."
function mg_newstep!()
    ngrs = gi(:ngrs); ngre = gi(:ngre)
    for g in 1:MGN[]
        copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g])
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]
        @inbounds for k in ngrs-2:ngre+2
            nom[k] = max(abs(F[k]), 1e-3 * gf(:cc) * 0.5 * (E[max(k-1, 1)] + E[k]), 1e-300)
        end
    end
    mg_slopes!()
    nothing
end
"Reset after a failed iterate: new := old."
mg_reset!() = (for g in 1:MGN[]; copyto!(MGEN[][g], MGEO[][g]); copyto!(MGFN[][g], MGFO[][g]); end; nothing)
"Apply the Newton correction: E_g ← E_g e^δ, F_g ← F_g + δ·nom (all blocks the solve returned)."
function mg_update!(rhs, a, b)
    @inbounds for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; iE = mgiE(g); iF = mgiF(g)
        for k in a:b; E[k] *= exp(rhs[iE, k]); end
        for k in a+1:b+1; F[k] += rhs[iF, k] * nom[k]; end
        mg_ghosts!(E, F, a, b)
    end
    nothing
end
"Ghost values of a group as its ghost rows: E copied, F antisymmetric at the centre (lribc = 2) and copied outside."
@inline function mg_ghosts!(E, F, a, b)
    E[a-1] = E[a]; E[b+1] = E[b]; E[b+2] = E[b+1]
    F[a] = 0.0; F[a-1] = gi(:lribc) == 2 ? -F[a+1] : F[a-1]; F[b+2] = F[b+1]
    nothing
end
"maxdel: clamp ln E_g corrections and return the convergence measure (|δ ln E_g| weighted by E_g/(E_g + MGEFL·E_r))."
function mg_dmax!(rhs, a, b)
    d = 0.0; cl = MGCLAMP[]; er = v(:ern)
    @inbounds for g in 1:MGN[]
        E = MGEN[][g]; iE = mgiE(g)
        for k in a:b
            δ = rhs[iE, k]
            isfinite(δ) || (rhs[iE, k] = 0.0; d = Inf; continue)
            rhs[iE, k] = clamp(δ, -cl, cl)
            w = E[k] / (E[k] + MGEFL[] * er[k])
            d = max(d, abs(rhs[iE, k]) * w)
        end
        for k in a:b+1; isfinite(rhs[mgiF(g), k]) || (rhs[mgiF(g), k] = 0.0; d = Inf); end
    end
    d
end

# ---------------------------------------------------------------- residual pieces
struct MGWS{T}
    egn::Vector{Vector{T}}; fgn::Vector{Vector{T}}; eg::Vector{Vector{T}}; fg::Vector{Vector{T}}
    chif::Vector{Vector{T}}; chifn::Vector{Vector{T}}; xke::Vector{Vector{T}}; emit::Vector{Vector{T}}; avchi::Vector{Vector{T}}
    erb::Vector{Vector{T}}; ersn::Vector{Vector{T}}; qe::Vector{Vector{T}}; qeo::Vector{Vector{T}}
    frb::Vector{Vector{T}}; frsn::Vector{Vector{T}}; qF::Vector{Vector{T}}; qFo::Vector{Vector{T}}
    frc::Vector{T}                                  # Σ_g χ̄_g F_g at interfaces (radiation force)
end
const MGWSC = Dict{DataType,Any}()
mgws(::Type{T}, mgr) where {T} = get!(MGWSC, T) do
    G = MGN[]; z() = [zeros(T, mgr) for _ in 1:G]
    MGWS{T}(z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), zeros(T, mgr))
end::MGWS{T}

"The group unknowns of the iterate (blocks k0 … ngre+2) into the workspace."
function mg_iterate!(M::MGWS{T}, x, neqn, nb, k0, ngre) where {T}
    @inbounds for g in 1:MGN[]
        en = M.egn[g]; fn = M.fgn[g]; E = MGEN[][g]; F = MGFN[][g]; nom = MGFNOM[][g]; iE = mgiE(g); iF = mgiF(g)
        for k in k0:ngre+2; en[k] = E[k]; fn[k] = F[k]; end
        for bb in 1:nb
            k = k0 + bb - 1
            en[k] = E[k] * exp(x[(bb-1)*neqn+iE]); fn[k] = F[k] + x[(bb-1)*neqn+iF] * nom[k]
        end
    end
    nothing
end

"""
θ-level group quantities after the grey material (chif, xke, xkp, plf at the θ level, tn): χ_g, κ_E,g, emission 4π κ_P,g B_g(T)
(θ-weighted B_g), interface opacity χ̄_g (harmonic, SLABFLUX, SUBZONE as the grey), and the advected face values.
"""
function mg_theta!(M::MGWS{T}, c, W, θ, ngrs, ngre, fle, flu) where {T}
    P = c.p; tn = W.tn; to = v(:to); dn = W.dn; dO = c.do_
    @inbounds for g in 1:MGN[]
        RE = MGRE[][g]; RP = MGRP[][g]; RR = MGRR[][g]; K1 = MGK1[][g]
        eg, fg, en, fn = M.eg[g], M.fg[g], M.egn[g], M.fgn[g]; EO = MGEO[][g]; FO = MGFO[][g]
        for k in ngrs-1:ngre+2
            eg[k] = θ * en[k] + (1 - θ) * EO[k]; fg[k] = θ * fn[k] + (1 - θ) * FO[k]
        end
        chi = M.chif[g]; chn = M.chifn[g]; ke = M.xke[g]; em = M.emit[g]
        for k in ngrs-1:ngre+1
            chi[k] = W.chif[k] * RR[k]; chn[k] = W.chifn[k] * RR[k]
            ke[k] = W.xke[k] * RE[k] + K1[k]
            # emission per unit mass: 4π κ_P,g B_g, B_g = b_g(T) σT⁴/π at each level, θ-weighted (κ_P at the θ level as the grey row)
            bn = mg_bfrac(g, tn[k]); bo = θ < 1 ? mg_bfrac(g, to[k]) : 0.0
            em[k] = 4P.cpi * W.xkp[k] * RP[k] * (θ * W.plfn[k] * bn + (1 - θ) * c.plfo[k] * bo)
        end
        av = M.avchi[g]
        for k in ngrs:ngre+1; av[k] = 2 / (1 / chi[k-1] + 1 / chi[k]); end
        SLABFLUX[] && sf_avchi!(av, chi, ngrs, ngre)
        SUBZONE[] && sz_avchi!(av, chi, W.d, W.dvol, W.rmu, ngrs, ngre)
        qe, qeo, qF, qFo = M.qe[g], M.qeo[g], M.qF[g], M.qFo[g]
        for k in ngrs-1:ngre+2
            qe[k] = en[k] / dn[k]; qeo[k] = EO[k] / dO[k]
            qF[k] = fn[k] * 2 / (dn[k] + dn[k-1]); qFo[k] = FO[k] * 2 / (dO[k] + dO[k-1])
        end
        advect!(M.erb[g], M.ersn[g], qe, qeo, MGESO[][g], fle, θ, P.cadv, P.eps, ngrs, ngre, true, 0)
        advect!(M.frb[g], M.frsn[g], qF, qFo, MGFSO[][g], flu, θ, P.cadv, P.eps, ngrs, ngre, false, 0)
    end
    fc = M.frc
    @inbounds for k in ngrs:ngre+1
        s = zero(T)
        for g in 1:MGN[]; s += M.avchi[g][k] * M.fg[g][k]; end
        fc[k] = s
    end
    nothing
end

"VEF_SPH outer-boundary energy density of a group (vef_esurf with the group's E, F, χ and background E_b,g)."
@inline function mg_esurf(rn, en, chn, dn, fn, fedd, N, xmu, cc, Eb)
    rc = 0.5 * (rn[N] + rn[N+1]); Ex = en[N] - Eb
    fx = ForwardDiff.value(Ex) != 0 ? clamp(ForwardDiff.value((fedd[N] * en[N] - Eb / 3) / Ex), 1 / 3, 1.0) : 1 / 3
    s = 0.5xmu * (3 - 1 / fx) * log(rn[N+1] / rc)
    Eb + exp(-s) * Ex - exp(-0.5s) * chn[N] * dn[N] * fn[N+1] * (rn[N+1] - rc) / (cc * fx)
end

"""
The radiation rows with MGROUPS: the grey rows ie/iF become E_r − Σ E_g, F_r − Σ F_g (interior zones / interfaces and the outer boundary
flux), and the group rows (energy per zone, momentum per interface, boundary and ghost rows as the grey ones).
"""
function mg_rows!(R, M::MGWS{T}, c, W, x, kb, P, dt, θ, xmu, cc, sph, tie) where {T}
    ngrs, ngre, neqn, k0, nb = c.ngrs, c.ngre, c.neqn, c.k0, c.nb
    ie, iF = c.ix[6], c.ix[7]
    X(i, k) = x[(kb(k) - 1) * neqn + i]
    ern, frn, rn, dn, dvoln, u, r, rmu, d, dvol, dmdt, urel = W.ern, W.frn, W.rn, W.dn, W.dvoln, W.u, W.r, W.rmu, W.d, W.dvol, W.dmdt, W.urel
    fedd = c.fedd; G = MGN[]
    # grey rows → sums
    @inbounds for k in ngrs:ngre
        s = zero(T); for g in 1:G; s += M.egn[g][k]; end
        R[ie, kb(k)] = ern[k] - s
    end
    @inbounds for k in ngrs+1:ngre+1
        s = zero(T); for g in 1:G; s += M.fgn[g][k]; end
        R[iF, kb(k)] = frn[k] - s
    end
    Tcmb = MGTCMB[]
    @inbounds for g in 1:G
        iE = mgiE(g); iFg = mgiF(g)
        en, fn, eg, fg = M.egn[g], M.fgn[g], M.eg[g], M.fg[g]; EO = MGEO[][g]; FO = MGFO[][g]
        ke = M.xke[g]; em = M.emit[g]; erb = M.erb[g]; frb = M.frb[g]; av = M.avchi[g]
        qg = g == MGQG[]
        # ---- energy
        for k in ngrs:ngre
            aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
            w = (en[k] * dvoln[k] - EO[k] * c.dvolo[k]) / dt + rmu[k+1] * fg[k+1] - rmu[k] * fg[k] +
                fedd[k] * eg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * eg[k] * aur * dvol[k] +
                excap(cc, ke[k], d[k], dt) * (cc * ke[k] * eg[k] - em[k]) * d[k] * dvol[k]
            if qg
                isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
                isempty(GASQ[]) || (w += GASQ[][k] * dvol[k])
            end
            if k == ngrs
                w += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * eg[k]
            elseif k == ngre
                w += dmdt[k] * erb[k]
                w += P.leobc != 3 ? urel[k+1] * rmu[k+1] * eg[k] : -dmdt[k+1] * erb[k+1]
            else
                w += -dmdt[k+1] * erb[k+1] + dmdt[k] * erb[k]
            end
            R[iE, kb(k)] = w
        end
        R[iE, 1] = X(iE, k0)
        R[iE, kb(ngrs-1)] = en[ngrs] - en[ngrs-1]
        R[iE, kb(ngre+1)] = -en[ngre] + en[ngre+1]
        R[iE, kb(ngre+2)] = -en[ngre+1] + en[ngre+2]
        # ---- momentum
        qF, qFo = M.qF[g], M.qFo[g]; chi = M.chif[g]
        for k in ngrs+1:ngre
            inert = (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                     qFo[k] * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                    ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2)
            drag = 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * av[k] * fg[k] / cc
            w = if !sph
                inert + rmu[k] * (fedd[k] * eg[k] - fedd[k-1] * eg[k-1]) + 0.5 * rmu[k] * (fg[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * eg[k-1] * dvol[k-1] + (3fedd[k] - 1) * eg[k] * dvol[k]) + drag
            else
                qs = vefq(fedd, r, xmu, k)
                inert + rmu[k] * (qs * fedd[k] * eg[k] - fedd[k-1] * eg[k-1] / qs) + 0.5 * rmu[k] * (fg[k] / cc^2) * (u[k+1] - u[k-1]) + drag
            end
            wt = tie ? TD.TIEW[k] : 0.0
            if wt > 0
                rcm = (0.5 * (r[k-1] + r[k]))^c.mu; rcp = (0.5 * (r[k] + r[k+1]))^c.mu
                hd = SLABTIE[] ? st_hdonor(max(TD.HEDD[k-1], 0.0), fedd[k-1], chi[k-1], chi[k], 0.5 * d[k-1] * dvol[k-1] / rmu[k], 0.5 * d[k] * dvol[k] / rmu[k]) :
                                 max(TD.HEDD[k-1], 0.0)
                B = rmu[k] * fg[k] / cc - (rcm * hd * eg[k-1] + rcp * min(TD.HEDD[k], 0.0) * eg[k])
                w = (1 - wt) * w + wt * B
            end
            R[iFg, kb(k)] = w
        end
        R[iFg, kb(ngrs)] = X(iFg, ngrs)
        bcmb = MGBCMB_[g]; ximg = P.ximr * bcmb
        R[iFg, kb(ngre+1)] = P.cpi * (2P.geddr + 1) * ximg + fn[ngre+1] -
                             cc * P.geddr * (sph ? mg_esurf(rn, en, M.chifn[g], dn, fn, fedd, ngre, xmu, cc, 4P.cpi * ximg / cc) : en[ngre])
        R[iFg, 1] = X(iFg, k0)
        k = ngrs - 1
        R[iFg, kb(k)] = P.lribc == 2 ? fn[k+2] + fn[k] : rn[k+1]^c.mu * fn[k+1] - rn[k]^c.mu * fn[k]
        R[iFg, kb(ngre+2)] = -fn[ngre+1] + fn[ngre+2]
    end
    nothing
end
"b_g(T_cmb) of the outer background (set with the edges)."
const MGBCMB_ = Float64[]

"""
Set up MGROUPS: edges [eV] (interior, ascending), the eV/k_B constant, T_cmb, the unknown offset i0 (9 with COUPLED). Allocates the arrays.
"""
function mg_setup!(edges, evk, tcmb, i0, mgr)
    issorted(edges) && all(>(0), edges) || error("MGROUPS: edges must be positive and ascending (eV): $edges")
    MG[] = true; MGN[] = length(edges) + 1; MGEDGE[] = collect(Float64, edges); MGEVK[] = evk; MGTCMB[] = tcmb; MGI0[] = i0
    empty!(MGBCMB_); for g in 1:MGN[]; push!(MGBCMB_, tcmb > 0 ? mg_bfrac(g, tcmb) : (g == 1 ? 1.0 : MGBFLOOR)); end
    mg_alloc!(mgr); empty!(MGWSC)
    nothing
end

"""
Initial group fields from the grey field (a restart from a grey checkpoint, or the IC): the CMB part partitioned at T_cmb, the excess at the
shape temperature Ts[k] (the driver passes T_r inside the photosphere and the star's T_eff outside), F_g = F_r b_g(Ts at the interface).
"""
function mg_init_from_grey!(Ts)
    a = gi(:ngrs); b = gi(:ngre); er = v(:ern); fr = v(:frn)
    arad = 4gf(:csigr) / gf(:cc); Ec = MGTCMB[] > 0 ? arad * MGTCMB[]^4 : 0.0
    for g in 1:MGN[]
        E = MGEN[][g]; F = MGFN[][g]; fill!(E, 0.0); fill!(F, 0.0)
        for k in a-2:b+2
            kk = clamp(k, a, b)
            Ex = max(er[k] - Ec, 0.0)
            E[k] = max(MGBCMB_[g] * min(er[k], Ec) + mg_bfrac(g, Ts[kk]) * Ex, MGBFLOOR * er[k])
        end
        for k in a+1:b+1; F[k] = fr[k] * mg_bfrac(g, Ts[min(k, b)]); end
        mg_ghosts!(E, F, a, b)
    end
    # exact sums (the rows E_r = Σ E_g, F_r = Σ F_g): renormalise
    for k in a-2:b+2
        s = sum(MGEN[][g][k] for g in 1:MGN[]); s > 0 && for g in 1:MGN[]; MGEN[][g][k] *= er[k] / s; end
    end
    for k in a+1:b+1
        s = sum(MGFN[][g][k] for g in 1:MGN[])
        d = fr[k] - s; MGFN[][1][k] += d          # (the remainder, round-off, to the IR group)
    end
    for g in 1:MGN[]; mg_ghosts!(MGEN[][g], MGFN[][g], a, b); end
    for g in 1:MGN[]; copyto!(MGEO[][g], MGEN[][g]); copyto!(MGFO[][g], MGFN[][g]); end
    mg_newstep!()
    nothing
end

"(diagnostic) group luminosities L_g = 4π r² F_g at interface k (values of the new level)."
mg_lum(k) = [4π * v(:rn)[k]^2 * MGFN[][g][k] for g in 1:MGN[]]
