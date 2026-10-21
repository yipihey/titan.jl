# S6b, Part 2: the optically thick radiative precursor of the accretion shock (steady state, free-fall flow, gas = radiation temperature),
# the SPS86 (eqs. 7) / OP03 treatment of the precursor as a diffusion region with advection, with the photosphere in the infalling envelope.
# Included by s6b_run.jl (option PRECTHICK=1; default off) after s6b_setup.jl.
#
# Steady spherical flow ρ v r² = −Ṁ/4π, v = −v_ff, total energy flux conservation (kinetic and gravitational energy cancel in free fall):
#     L_co(r) = L_∞ + Ṁ (h(r) − h_∞) + (16π/3) r² E |v|       (diffusive luminosity in the comoving frame)
# h = gas enthalpy (with the ionisation energy, StellarEOS), E = a T⁴ the radiation energy density (T_gas = T_rad where the precursor is thick). The last
# term is the radiation enthalpy flux 4/3 E v of the inflow (trapping: the radiation is carried inward faster than it diffuses out when τ ≳ c/v).
# At the surface (shock at R) the energy balance of a star that accretes cold gas of enthalpy h_s = h(R) is L_∞ = L* + Ṁ (G M/R − h_s): the energy used to
# heat and ionise the gas is not radiated (for M = 0.22 M☉, R = 27 R☉ the ionisation energy of the accreted gas is ≈ G M/R).
# Closure: diffusion, dE/dr = −3κρ L_co/(4πr²c), integrated inward from the Lucy photosphere (τ = 2/3 from infinity, where E = E_Lucy). Outside the photosphere
# and in the thin limit the radiation field is Lucy's law (s6b_setup.jl `precursor`); the thick-precursor physics enters as the difference of two diffusion integrals,
#     E_s = E_Lucy(R) × E_diff(R; advection, enthalpy) / E_diff(R; neither),
# which is 1 where advection and enthalpy do not matter (the S6a result is unchanged exactly there; no flux-limited closure needed in the thin part).

"Gas enthalpy table h(ρ, T) = e + P/ρ (gas only, per gram, relative to the neutral gas at T_ref), bilinear in (ln ρ, ln T)."
struct HTable; lρ::Vector{Float64}; lT::Vector{Float64}; h::Matrix{Float64}; end
const HTAB = Dict{Float64,HTable}()
const HT_REF = 2000.0
function htable(X)
    get!(HTAB, round(X, digits = 4)) do
        lρ = collect(range(log(1e-15), log(1e-1); length = 57)); lT = collect(range(log(2000.0), log(3e6); length = 90))
        h = zeros(length(lρ), length(lT))
        for (i, lr) in enumerate(lρ), (j, lt) in enumerate(lT)
            p, e = StellarEOS.gas_pe(exp(lr), exp(lt); X, Y = 1 - X, Z = 0.0, pairs = false); h[i, j] = e + p / exp(lr)
        end
        for i in 1:length(lρ); h[i, :] .-= h[i, 1]; end          # zero at T_ref (neutral gas)
        HTable(lρ, lT, h)
    end
end
function hgas(tb::HTable, ρ, T)
    x = (log(ρ) - tb.lρ[1]) / (tb.lρ[2] - tb.lρ[1]); y = (log(T) - tb.lT[1]) / (tb.lT[2] - tb.lT[1])
    i = clamp(floor(Int, x) + 1, 1, length(tb.lρ) - 1); j = clamp(floor(Int, y) + 1, 1, length(tb.lT) - 1)
    wx = clamp(x - (i - 1), 0.0, 1.0); wy = clamp(y - (j - 1), 0.0, 1.0)
    (1 - wx) * ((1 - wy) * tb.h[i, j] + wy * tb.h[i, j+1]) + wx * ((1 - wy) * tb.h[i+1, j] + wy * tb.h[i+1, j+1])
end

"""
    precursor_E(R, M, Mdot, Lin, hR; X, advect, rmax, dEmax, dτmax) -> (; Es, Ts, tau, rph, hgas_R, Linf, nsub)

Radiation energy density at the shock, E_s = E_thin(R) + E_d(R), by inward integration from rmax·R with adaptive substeps. E_thin(r) = a T_L⁴ W(r), W = ½(1 − √(1 − R²/r²)),
T_L⁴ = L_∞/(4πσR²): the free-streaming (geometric dilution) part of Lucy's law; E_d the diffusive part, dE_d/dr = −(3κρ/4πr²c)[L_∞ + Ṁ h(r)] − (4κρ|v|/c)(E_thin + E_d),
E_d(rmax) = 0 (Lucy's law is the case Ṁ h = 0, no advection and the plain optical depth instead of ∫(R/r)²κρ dr: the weighting here is the physical one). L_∞ = Lin − Ṁ h_R (`hR = nothing`:
L_∞ = Lin, no enthalpy term). κ at T = max((E/a)^¼, 2.5e3 K). `rph`: τ = 2/3 from infinity (Rosseland, same T), `hgas_R` = h_gas(ρ_s, T_s).
"""
function precursor_E(R, M, Mdot, Lin, hR; X = 0.75, advect = true, rmax = 50.0, dEmax = 0.2, dτmax = 0.3, tfloor = 2.5e3)
    c = Z0Physics.c; a = Z0Physics.arad; σ = Z0Physics.σSB; tb = htable(X)
    hterm = hR !== nothing
    Linf = hterm ? max(Lin - Mdot * hR, 1e-6 * Lin) : Lin
    TL4 = Linf / (4π * σ * R^2)
    ρf(ri) = Mdot / (4π * ri^2 * sqrt(2GRAV * M / ri)); vf(ri) = sqrt(2GRAV * M / ri)
    Wf(ri) = 0.5 * (1 - sqrt(max(1 - (R / ri)^2, 0.0))); Eth(ri) = a * TL4 * Wf(ri)
    rc = rmax * R; Ed = 0.0; τ = 0.0; rph = NaN; nsub = 0
    while rc > R * (1 + 1e-12)
        E0 = Eth(rc) + Ed; ρ0 = ρf(rc); T0 = max((E0 / a)^0.25, tfloor)
        κ0 = StellarOpacity.kappa_R(ρ0, T0, X); h0_ = hterm ? hgas(tb, ρ0, T0) : 0.0
        Lb0 = Linf + Mdot * h0_
        A0 = 3κ0 * ρ0 / (4π * rc^2 * c); B0 = advect ? 4κ0 * ρ0 * vf(rc) / c : 0.0
        slope = A0 * Lb0 + B0 * E0
        Δs = min(rc - R, dEmax * E0 / max(slope, 1e-300), dτmax / max(κ0 * ρ0, 1e-300), 0.05 * rc)
        rm = rc - 0.5Δs; ρm = ρf(rm); Eth_m = Eth(rm)
        Edm = Ed + 0.5Δs * slope; Tm = max(((Eth_m + Edm) / a)^0.25, tfloor)
        κm = StellarOpacity.kappa_R(ρm, Tm, X); hm = hterm ? hgas(tb, ρm, Tm) : 0.0
        A = 3κm * ρm / (4π * rm^2 * c); B = advect ? 4κm * ρm * vf(rm) / c : 0.0
        src = A * (Linf + Mdot * hm) + B * Eth_m
        Ednew = B * Δs < 1e-8 ? Ed + src * Δs : Ed * exp(B * Δs) + (src / B) * expm1(B * Δs)
        dτ = κm * ρm * Δs
        isnan(rph) && τ + dτ >= 2 / 3 && (rph = rc - (2 / 3 - τ) / dτ * Δs)
        τ += dτ; Ed = Ednew; rc -= Δs; nsub += 1
        nsub > 50000 && break
    end
    Es = Eth(R) + Ed; Ts = (Es / a)^0.25
    (; Es, Ts, tau = τ, rph = isnan(rph) ? R : rph, hgas_R = hterm ? hgas(tb, ρf(R), max(Ts, tfloor)) : 0.0, Linf, nsub)
end

"""
Self-consistent inflow enthalpy h_R = h_gas(ρ_s, T_s(h_R)) and the surface radiation temperature. g(h) = h_gas(T_s(h)) − h is strictly decreasing, but T_s(h) has a **jump**: for h below
h_c the ionisation front of the precursor lies outside the shock (hot, τ ~ 10²–10³, T_s ~ 3e4 K), above it inside (cold, thin precursor, T_s ~ 6e3 K), and h_c ≈ 0.6 G M/R at M = 0.22 M☉
(there the energy of the shock just ionises the infall: a first-order ionisation-front transition). The root of g is at the jump; the steady structure there is a front at the shock,
with any T_s in between. The jump moves smoothly with R, L*, Ṁ, but E_s(R) taken from the structure on either side is a cliff (S6a's stall). Here the root h_R (bisection on the sign of g, continuous)
is combined with the thermostat condition h_gas(ρ_s, T_s) = h_R: **T_s = T(h_R)** (inverse of the gas enthalpy at the pre-shock density), smooth in all parameters and equal to the
structure's T_s wherever g has a smooth root. Returns (the structure on the cold side of the root, h_R, T_s).
"""
function solve_thick(R, M, Mdot, Lin; X = 0.75, advect = true, rmax = 50.0, tol = 2e-3, nbis = 40)
    tb = htable(X); ρR = Mdot / (4π * R^2 * sqrt(2GRAV * M / R)); a = Z0Physics.arad
    f(h) = (d = precursor_E(R, M, Mdot, Lin, h; X, advect, rmax); (hgas(tb, ρR, max((d.Es / a)^0.25, 2.5e3)) - h, d))
    hmax = 0.995 * Lin / Mdot                      # L_inf > 0
    a_ = 0.0; fa, da = f(a_)
    hb = 0.0; dcold = da
    if fa > 0
        b_ = hmax; fb, db = f(b_); dcold = db; hb = b_
        if fb < 0
            for _ in 1:nbis
                h = 0.5 * (a_ + b_); fh, dh = f(h)
                fh > 0 ? (a_ = h) : (b_ = h; dcold = dh)
                hb = 0.5 * (a_ + b_)
                b_ - a_ < tol * max(h, 1e-2 * GRAV * M / R) && break
            end
        end
    end
    Ts = temperature_of_h(tb, ρR, hb)
    (dcold, hb, Ts)
end
"T with h_gas(ρ, T) = h (bisection; h_gas increases with T); clamped to [2.5e3 K, 3e6 K]."
function temperature_of_h(tb, ρ, h)
    lo = log(2.5e3); hi = log(3e6)
    hgas(tb, ρ, exp(lo)) >= h && return exp(lo)
    for _ in 1:50
        m_ = 0.5 * (lo + hi); hgas(tb, ρ, exp(m_)) < h ? (lo = m_) : (hi = m_)
        hi - lo < 1e-6 && break
    end
    exp(0.5 * (lo + hi))
end

"""
    precursor_thick(R, M, Mdot, Lstar; X, ...) -> (; Ts, tau, Tph, rph, Es, Es_lucy, Linf, hR, w)

Lin = L* + G M Ṁ/R. Two limits, blended in ln E_s with the weight w = smoothstep((ln f − ln f1)/(ln f2 − ln f1)), f = h_R/(G M/R) the fraction of the shock energy that the
accreted gas keeps as enthalpy (f1 = 0.01, f2 = 0.05):
  (a) f ≫ 1: ionisation-limited (M ≲ 3 M☉ at Ṁ ≳ 1e-3, or Ṁ ≳ 4e-3): `solve_thick`: T_s = T(h_R) with the root h_R of the energy balance (L_∞ = Lin − Ṁ h_R), see there;
  (b) f ≪ 1: the energy of the shock is not used up by the gas: S6a's Lucy law for Lin times the trapping factor E_adv/E_off of `precursor_E` (advection of the radiation with and without; both
      with the same weighting, so the factor is 1 for a thin precursor and ≈ e^{4τ v/c} for a thick one). The weighting 1/r² of the diffusive flux in `precursor_E` is more accurate than Lucy's
      plain optical depth, but it lowers E_s of the thick extended precursors (M ≳ 30 M☉ at Ṁ = 1.1e-3) by 4 (T_s by 0.7) and with it the ZAMS radius changes (a first version used it everywhere; R*(40 M☉) = 9.6 vs 4.6 R☉
      of S6a): kept as Lucy there so that the validated S6a results do not change.
T_eff,ph = (L_∞/4πr_ph²σ)^¼, r_ph: τ = 2/3 from infinity (`precursor_E`). `hterm = false`: only (b); `advect = false`: no trapping.
"""
function precursor_thick(R, M, Mdot, Lstar; X = 0.75, rmax = 50.0, n = 300, advect = true, hterm = true, h0 = nothing, f1 = 0.01, f2 = 0.05, kw...)
    Lin = Lstar + GRAV * M * Mdot / R
    pl = precursor(R, M, Mdot, Lin; X, rmax, n)
    (advect || hterm) || return merge(pl, (; Es_lucy = pl.Es, Linf = Lin, hR = 0.0, w = 0.0))
    # (b) Lucy × trapping
    Elb = pl.Es
    if advect
        eoff = precursor_E(R, M, Mdot, Lin, nothing; X, advect = false, rmax); eadv = precursor_E(R, M, Mdot, Lin, nothing; X, advect = true, rmax)
        Elb = pl.Es * eadv.Es / eoff.Es
    end
    w = 0.0; hR = 0.0; on = nothing; Ts = (Elb / Z0Physics.arad)^0.25; Es = Elb
    if hterm
        on, hR, Ts_t = solve_thick(R, M, Mdot, Lin; X, advect, rmax)
        f = hR / (GRAV * M / R)
        x = clamp((log(max(f, 1e-30)) - log(f1)) / (log(f2) - log(f1)), 0.0, 1.0); w = x * x * (3 - 2x)
        Es = exp((1 - w) * log(Elb) + w * log(Z0Physics.arad * Ts_t^4)); Ts = (Es / Z0Physics.arad)^0.25
    else
        on = precursor_E(R, M, Mdot, Lin, nothing; X, advect, rmax)
    end
    Linf = max(Lin - Mdot * hR, 1e-6 * Lin)
    Tph = (Linf / (4π * on.rph^2 * Z0Physics.σSB))^0.25
    (; Ts, tau = on.tau, Tph, rph = on.rph, Es, Es_lucy = pl.Es, Linf, hR, w)
end
