# Outer boundary conditions.
#
# (A) No infall: Eddington grey photosphere at τ = 2/3 (same as MESA's default T–τ atmosphere).
#     T_b = T_eff,  P_tot,b = (2/3)(g/κ + F/c)   ⇔   P_gas,b = (2/3)(g/κ)(1 − Γ_κ),  Γ_κ = κF/(c g)
#
# (B) Spherical infall at free fall onto a radiative accretion shock at r = R.
#     Upstream: ρ_in(r) = Ṁ / (4π r² v_ff(r)), v_ff = sqrt(2GM/r)  (∝ r^-3/2).
#     The post-shock gas cools in a layer thin compared with R (checked separately) and settles
#     onto the envelope. Momentum balance across shock + settling layer:
#        P_tot,b = ρ_in v_ff² + (2/3) g/κ_b + (F_tot/c)(τ_in + 2/3)
#     Temperature (Eddington approximation through infall column τ_in + stellar τ = 2/3):
#        T_b⁴ = (3/4) T_L⁴ (τ_in + 4/3),   T_L⁴ = L_tot / (4πσR²),  L_tot = L + G M Ṁ / R
#     τ_in is computed from the LTE opacity of the infall with the Lucy (1971) temperature
#     T⁴(r) = T_L⁴ [W(r) + (3/4) τ(r)], W = ½(1 − sqrt(1 − (R/r)²)), integrated inward from ∞.
#     In the limit Ṁ → 0, (B) reduces exactly to (A).

vff(M, r) = sqrt(2G * M / r)
rho_infall(Mdot, M, r) = Mdot / (4π * r^2 * vff(M, r))

"LTE optical depth of the infall from ∞ to R and the photosphere radius within it (if any)."
function infall_tau(Mdot, M, R, Ltot, tab; rmax_fac = 1e4, n = 4000, κ_override = nothing)
    Mdot <= 0 && return (τ = 0.0, r_ph = R)
    TL4 = Ltot / (4π * σSB * R^2)
    rs = exp.(range(log(rmax_fac * R), log(R), length = n))
    τ = 0.0; r_ph = R; found = false
    for i in 1:n-1
        r1, r2 = rs[i], rs[i+1]; rm = sqrt(r1 * r2)
        ρ = rho_infall(Mdot, M, rm)
        W = 0.5 * (1 - sqrt(max(1 - (R / rm)^2, 0.0)))
        T = (TL4 * (W + 0.75τ))^0.25
        κ = κ_override === nothing ? kap(tab, ρ, T) : κ_override
        dτ = κ * ρ * (r1 - r2)
        if !found && τ + dτ >= 2 / 3
            r_ph = r1 - (2 / 3 - τ) / (κ * ρ); found = true
        end
        τ += dτ
    end
    return (τ = τ, r_ph = found ? r_ph : R)
end

"Solve for (P_b, T_b) at radius R given infall rate (Mdot = 0 gives the plain photosphere)."
function outer_bc(R, p::QSParams; Mdot = 0.0, κ_in = nothing)
    Lacc = Mdot > 0 ? G * p.M * Mdot / R : 0.0
    Ltot = p.L + Lacc
    inf = infall_tau(Mdot, p.M, R, Ltot, p.tab; κ_override = κ_in)
    TL4 = Ltot / (4π * σSB * R^2)
    Tb = (0.75 * TL4 * (inf.τ + 4 / 3))^0.25
    g = G * p.M / R^2
    F = Ltot / (4π * R^2)
    ram = Mdot > 0 ? rho_infall(Mdot, p.M, R) * vff(p.M, R)^2 : 0.0
    Prad_b = Prad(Tb)
    # f(P) = P − [ram + (2/3)g/κ(ρ(P,Tb)) + (F/c)(τ_in + 2/3)];  P must exceed Prad(Tb)
    fP(P) = (ρ = rho_PT(P, Tb); isnan(ρ) ? -Inf : P - (ram + (2 / 3) * g / kap(p.tab, ρ, Tb) + F / c * (inf.τ + 2 / 3)))
    lo = log(Prad_b * (1 + 1e-10)); hi = log(Prad_b + 1e8)
    fP(exp(hi)) < 0 && return nothing
    if fP(exp(lo)) > 0
        return nothing
    end
    for _ in 1:200
        mid = 0.5(lo + hi)
        fP(exp(mid)) > 0 ? (hi = mid) : (lo = mid)
        hi - lo < 1e-12 && break
    end
    Pb = exp(0.5(lo + hi))
    κb = kap(p.tab, rho_PT(Pb, Tb), Tb)
    return (Pb = Pb, Tb = Tb, τ_in = inf.τ, r_ph_infall = inf.r_ph, Ltot = Ltot, ram = ram,
            Teff_L = TL4^0.25, Γκ = κb * F / (c * g), Pgas_b = Pb - Prad_b)
end

"Mismatch for the shooting on R: +1 if envelope mass left at r_in exceeds MBH, −1 if exhausted."
function shoot_sign(R, p; Mdot = 0.0, κ_in = nothing)
    bc = outer_bc(R, p; Mdot = Mdot, κ_in = κ_in)
    bc === nothing && return (0, NaN, bc)
    st, mrin, _ = integrate_inward(R, bc.Pb, bc.Tb, p)
    st == :failed && return (0, NaN, bc)
    st == :mass_exhausted && return (-1, mrin, bc)
    return (mrin > p.MBH ? 1 : -1, mrin, bc)
end

"""
Find R such that m(ε_in R) = MBH. Scans log R on a grid to bracket a sign change, then bisects.
Returns the radius and boundary solution, or nothing.
"""
function solve_R(p::QSParams; Mdot = 0.0, κ_in = nothing, lRmin = 14.5, lRmax = 17.5, nscan = 31, verbose = false)
    lRs = range(lRmin, lRmax, length = nscan)
    sg = [shoot_sign(10^lR, p; Mdot = Mdot, κ_in = κ_in)[1] for lR in lRs]
    verbose && println("scan signs: ", sg)
    idx = findfirst(i -> sg[i] != 0 && sg[i+1] != 0 && sg[i] != sg[i+1], 1:nscan-1)
    idx === nothing && return nothing
    a, b = lRs[idx], lRs[idx+1]; sa = sg[idx]
    for _ in 1:60
        mid = 0.5(a + b)
        s = shoot_sign(10^mid, p; Mdot = Mdot, κ_in = κ_in)[1]
        s == 0 && (b = mid; continue)
        s == sa ? (a = mid) : (b = mid)
        b - a < 1e-7 && break
    end
    R = 10^(0.5(a + b))
    bc = outer_bc(R, p; Mdot = Mdot, κ_in = κ_in)
    _, _, prof = integrate_inward(R, bc.Pb, bc.Tb, p; save = true)
    return (R = R, bc = bc, prof = prof, signs = sg)
end
