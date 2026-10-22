# S7a: pure-hydrogen H II region chemistry/energy coupled to the photon-conserving direct transfer of RadTransfer (inside-out sweep, time-averaged
# neutral fraction per zone as C²-Ray, Mellema et al. 2006). Used by the Iliev et al. (2006, 2009) tests in s7a_test.jl and the hydro driver s7a_hydro.jl.
# Assumes ChemNet, RadTransfer are loaded.
module HII
import ..ChemNet, ..RadTransfer
const CN = ChemNet; const RT = RadTransfer
const kB = CN.kB; const mH = CN.mH; const hP = CN.hP; const χ = CN.χH

"Case B recombination coefficient of H, Hui & Gnedin (1997): 2.753e-14 λ^1.5/(1 + (λ/2.74)^0.407)^2.242, λ = 315614/T (2.59e-13 at 1e4 K)."
aB_HG(T) = (λ = 315614.0 / T; 2.753e-14 * λ^1.5 / (1 + (λ / 2.74)^0.407)^2.242)
const AB_FIXED = Ref(0.0)             # > 0: constant α_B (Iliev Test 1: 2.59e-13)
aB(T) = AB_FIXED[] > 0 ? AB_FIXED[] : aB_HG(T)
"Collisional ionisation of H by electrons, Cen (1992): 5.85e-11 √T/(1 + √(T/1e5)) e^{−157809.1/T} [cm³ s⁻¹]."
kci(T) = 5.85e-11 * sqrt(T) / (1 + sqrt(T / 1e5)) * exp(-157809.1 / T)

# mean kinetic energy radiated per case-B recombination, from the Milne integrals of the H levels n = 2…15 (Johnson/Kramers σ of AtomicHHe):
# ⟨E_kin⟩(T) = Σ_n ∫ (hν − χ/n²) σ_n(ν) ν² e^{−hν/kT} dν / Σ_n ∫ σ_n ν² e^{−hν/kT} dν  (T-table, log T)
const EKTAB = Ref{Any}(nothing)
function ektab(levels)
    lT = collect(range(1.0, 8.0, length = 141)); out = zeros(length(lT))
    for (i, x) in enumerate(lT)
        T = 10^x; num = 0.0; den = 0.0
        for (Ex, g, Eth, σf) in levels
            ν0 = Eth / hP; νs = [ν0 * exp(y) for y in range(0.0, log(1 + 60kB * T / Eth + 30.0), length = 401)]
            for j in 1:400
                ν = 0.5 * (νs[j] + νs[j+1]); dν = νs[j+1] - νs[j]
                w = g * σf(ν) * ν^2 * exp(-hP * (ν - ν0) / (kB * T)) * dν
                num += w * hP * (ν - ν0); den += w          # kinetic energy of the captured electron = hν − E_th(n)
            end
        end
        out[i] = num / den
    end
    (lT, out)
end
"⟨E_kin⟩/kT of a case-B recombination at T (levels n ≥ 2 passed once to `init_ek`)."
function ekin(T)
    lT, e = EKTAB[]; x = clamp(log10(T), lT[1], lT[end]); i = clamp(searchsortedlast(lT, x), 1, length(lT) - 1)
    w = (x - lT[i]) / (lT[i+1] - lT[i]); ((1 - w) * e[i] + w * e[i+1])
end
init_ek(levels_n2up) = (EKTAB[] = ektab(levels_n2up); nothing)

"Cooling [erg cm⁻³ s⁻¹] of pure H at T, ionised fraction x: Cen (1992) collisional excitation and ionisation, case-B recombination ⟨E_kin⟩ α_B, bremsstrahlung (g_ff = 1.3), optional Compton."
function Λ(T, x, nH; compton = 0.0)
    ne = x * nH; nHI = (1 - x) * nH
    exc = 7.5e-19 / (1 + sqrt(T / 1e5)) * exp(-118348.0 / T) * ne * nHI
    ion = 1.27e-21 * sqrt(T) / (1 + sqrt(T / 1e5)) * exp(-157809.1 / T) * ne * nHI
    rec = aB(T) * ekin(T) * ne * ne
    ff = 1.42e-27 * 1.3 * sqrt(T) * ne * ne
    exc + ion + rec + ff
end

"Root of an increasing F on [lo, hi] (F(lo) < 0 < F(hi) given): Illinois regula falsi, returns the last iterate."
function illinois(F, lo, hi, flo, fhi)
    mid = 0.5 * (lo + hi); side = 0
    for _ in 1:100
        mid = (lo * fhi - hi * flo) / (fhi - flo); fm = F(mid)
        abs(fm) < 1e-15 && return mid
        if fm < 0
            lo = mid; flo = fm; side == -1 && (fhi *= 0.5); side = -1
        else
            hi = mid; fhi = fm; side == 1 && (flo *= 0.5); side = 1
        end
        hi - lo < 1e-14 * max(hi, 1e-30) && return mid
    end
    mid
end
mutable struct State
    x::Vector{Float64}; T::Vector{Float64}
end

"""
    xode(x0, Γ, c, a, dt) -> (x_new, x̄)

Exact solution and time average of dx/dt = (1 − x)(Γ + c x) − a x² = S0 + S1 x − S2 x² (constant coefficients) over dt: roots x₊ > 0 > x₋ of
S2 x² − S1 x − S0, x(t) = x₋ + (x₊ − x₋)/(1 − C e^{−Dt}), C = (x₀ − x₊)/(x₀ − x₋), D = √(S1² + 4 S2 S0), and the time average
x̄ = x₊ + ln[(1 − C e^{−D dt})/(1 − C)]/(S2 dt) (C²-Ray, Mellema et al. 2006).
"""
function xode(x0, Γ, c, a, dt)
    S0 = Γ; S1 = c - Γ; S2 = c + a; D = sqrt(S1^2 + 4S2 * S0)
    if S0 == 0 && S1 == 0                                                   # pure recombination: dx/dt = −S2 x²
        return (x0 / (1 + S2 * x0 * dt), log1p(S2 * x0 * dt) / (S2 * dt))
    end
    xp = S1 >= 0 ? (S1 + D) / (2S2) : 2S0 / (D - S1)                   # stable forms of the two roots
    xm = S1 >= 0 ? -2S0 / (S1 + D) : (S1 - D) / (2S2)
    C = (x0 - xp) / (x0 - xm); E = exp(-D * dt)
    xn = xm + (xp - xm) / (1 - C * E)
    xb = xp + log((1 - C * E) / (1 - C)) / (S2 * dt)
    (clamp(xn, 0.0, 1.0), clamp(xb, 0.0, 1.0))
end

"""
    sweep!(st, g, σ, src, r, nH, dt; energy, collisions, ...) -> NamedTuple

Advance every zone over dt, inside out, with the photon-conserving C²-Ray scheme (Mellema et al. 2006): σ[f] is the bin-averaged H I cross section on the
grid g, src.L the source power per bin. For each zone the time-averaged ionised fraction x̄ is the fixed point of x̄ = ⟨x⟩(Γ(x̄)), where
Γ = Σ_f (L_in,f/hν_f) (1 − e^{−τ_f})/(n_HI V), τ_f = n_H(1 − x̄) σ_f Δr, is the photo-rate per neutral atom that absorbs exactly the photons entering the zone
(so the ionisations over dt equal the photons absorbed over dt) and ⟨x⟩ the exact time average of the ODE with constant Γ, c = n_H k_ci, a = α_B n_H (`xode`).
The transmitted power is L_in e^{−τ(x̄)}. Energy: (3/2) n_H k [(1 + x_new) T_new − (1 + x_old) T_old] = dt (G − Λ(T_new, x̄)), G the photo-heating
of the absorbed photons (hν − χ). Returns the photon ledger (per second, time-averaged over the step).
"""
function sweep!(st::State, g, σ, src, r, nH, dt; energy = true, collisions = true, Tfloor = 1.0, xmin = 1e-12, nref = 3)
    n = length(st.x); nf = length(σ); Lin = copy(src.L)
    Nabs = 0.0; Gtot = 0.0
    xold = copy(st.x); Told = copy(st.T)
    for k in 1:n
        V = 4π / 3 * (r[k+1]^3 - r[k]^3); Δr = r[k+1] - r[k]; nHk = nH[k]
        x0 = max(xold[k], xmin); T0 = Told[k]
        Tn = T0; xn = x0; xbar = x0; Sphot = 0.0; G = 0.0; Lout = similar(Lin)
        for it in 1:nref
            Tr = 0.5 * (T0 + Tn)
            a = aB(Tr) * nHk; c = collisions ? kci(Tr) * nHk : 0.0
            # Γ per neutral atom (stable at x̄ → 1: (1 − e^{−τ})/τ)
            Γof(xb) = begin
                nHI = nHk * (1 - xb); s = 0.0
                for f in 1:nf
                    τ = nHI * σ[f] * Δr; s += Lin[f] / (hP * g.ν[f]) * σ[f] * Δr / V * (τ < 1e-12 ? 1.0 - 0.5τ : -expm1(-τ) / τ)
                end
                s
            end
            if all(f -> Lin[f] / (hP * g.ν[f]) < 1e-30 * nHk * V, 1:nf)
                xn, xbar = xode(x0, 0.0, c, a, dt)
            else
                F(xb) = xb - xode(x0, Γof(xb), c, a, dt)[2]
                lo = 0.0; hi = 1.0; flo = F(lo); fhi = F(hi)
                xbar = flo >= 0 ? lo : fhi <= 0 ? hi : illinois(F, lo, hi, flo, fhi)
                xn = xode(x0, Γof(xbar), c, a, dt)[1]
            end
            nHI = nHk * (1 - xbar)
            Sphot = 0.0; G = 0.0
            for f in 1:nf
                τ = nHI * σ[f] * Δr; ab = Lin[f] * (-expm1(-τ)); Sphot += ab / (hP * g.ν[f]); G += ab * (1 - χ / (hP * g.ν[f]))
            end
            G /= V
            if energy
                E0 = 1.5 * nHk * kB * (1 + x0) * T0
                resid(T) = 1.5 * nHk * kB * (1 + xn) * T - E0 - dt * (G - Λ(T, xbar, nHk))
                lo = Tfloor; hi = max(2T0, 1e5)
                while resid(hi) < 0 && hi < 1e9; hi *= 2; end
                resid(lo) > 0 && (hi = lo)
                if hi > lo
                    for _ in 1:80
                        mid = sqrt(lo * hi); resid(mid) < 0 ? (lo = mid) : (hi = mid)
                        hi / lo - 1 < 1e-13 && break
                    end
                end
                Tn = sqrt(lo * hi)
            end
        end
        nHI = nHk * (1 - xbar)
        for f in 1:nf
            τ = nHI * σ[f] * Δr; Lout[f] = Lin[f] * exp(-τ)
        end
        Nabs += Sphot; Gtot += G * V
        st.x[k] = max(xn, xmin); st.T[k] = Tn
        Lin = Lout
    end
    Nesc = sum(Lin[f] / (hP * g.ν[f]) for f in 1:nf)
    (Nsrc = sum(src.L[f] / (hP * g.ν[f]) for f in 1:nf), Nabs = Nabs, Nesc = Nesc, Eheat = Gtot)
end

"Position of the I-front: radius where x = 0.5 (linear interpolation of x between zone centres, outermost crossing from inside)."
function rfront(st::State, r)
    n = length(st.x); rc = [0.5 * (r[k] + r[k+1]) for k in 1:n]
    kf = 0
    for k in 1:n; st.x[k] >= 0.5 && (kf = k); end
    kf == 0 && return r[1]; kf == n && return r[end]
    w = (st.x[kf] - 0.5) / (st.x[kf] - st.x[kf+1]); rc[kf] + w * (rc[kf+1] - rc[kf])
end
"I-front radius from the position where the ionised-photon budget closes: the radius containing ∫ x dV (volume-integrated ionised fraction → equivalent sharp-front radius)."
function rfront_vol(st::State, r, nH)
    n = length(st.x); Ni = sum(st.x[k] * 4π / 3 * (r[k+1]^3 - r[k]^3) for k in 1:n); cbrt(3Ni / (4π))
end
end # module
