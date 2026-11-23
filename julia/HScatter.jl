# Rayleigh and Raman scattering of H(1s) between Lyman α and the Lyman limit (docs/POPIII_MULTIGROUP.md §16).
#
# Source: Kokubo 2024, MNRAS 529, 2131 (2024MNRAS.529.2131K): the Kramers–Heisenberg–Waller cross sections with the analytic radial
# integrals of the bound and continuum intermediate p states. With ε = λ_LyLimit/λ = hν/χ_H:
#   σ_1s/σ_T  = | Σ_{n≥2} g_n [ε/(a_n − ε) − ε/(a_n + ε)] + ∫₀^∞ dn′ g′(n′) [ε/(b − ε) − ε/(b + ε)] |²      (his eq. for σ_1s, second form)
#     a_n = 1 − 1/n², b = 1 + 1/n′², g_n = 2⁷3⁻¹ n⁵ (n²−1)⁻⁴ ((n−1)/(n+1))^{2n},
#     g′ = 2⁷3⁻¹ n′⁵ (n′²+1)⁻⁴ e^{−4n′ arccot n′}/(1 − e^{−2πn′})       (2 g_n = f_n, 2 g′ = df/dn′: the oscillator strengths)
#   σ_nBs/σ_T = (1 − (1−n_B⁻²)/ε) | Σ_n c_n [1/(a_n − ε) + 1/(n_B⁻² − 1/n² + ε)] + ∫ dn′ c′ [1/(b − ε) + 1/(n_B⁻² + 1/n′² + ε)] |²,
#   the same for σ_nBd (n_B = 3, 4) with their coefficients, each zero below its threshold ε = 1 − n_B⁻².
# Raman branches to n_B ≥ 5 (open only above Ly δ, ε > 0.96) are not in the source and are left out.
# The resonances are undamped in the source; here every bound term carries the radiative width of np → 1s (Γ_1n of the source's footnote,
# Bethe & Salpeter 1957) as 1/(a_n − ε − iγ_n/2): it only regularises the poles, which a frequency grid does not resolve anyway
# (|Δλ| ≲ 1e-4 Å; the line cores are the business of the line treatment, not of this continuum).
# Checks (hscatter_check): Thomas–Reiche–Kuhn sum Σ f_n + ∫ df = 1; ε → 0: σ_1s/σ_T → (9/8)² ε⁴ (static polarisability 9/2 a₀³).
module HScatter

const σT = 6.6524587e-25
const ALPHA = 7.2973525693e-3
const NMAX = 1000                                  # bound intermediate states summed explicitly (the source: 1000)
const NQ = 6000                                    # quadrature points of the continuum integral
const EPS0 = 0.80; const EPS1 = 0.9996; const NE = 8000   # table: uniform in ε (LW band and a margin)

acot(x) = atan(1 / x)
"Continuum quadrature: n′ = t/(1−t) on midpoints of t ∈ (0,1); returns (n′, weight dn′)."
function cquad()
    n = zeros(NQ); w = zeros(NQ)
    for i in 1:NQ; t = (i - 0.5) / NQ; n[i] = t / (1 - t); w[i] = 1 / (NQ * (1 - t)^2); end
    n, w
end
gb(n) = 2.0^7 / 3 * n^5 / (n^2 - 1)^4 * ((n - 1) / (n + 1))^(2n)
gc(n) = 2.0^7 / 3 * n^5 / (n^2 + 1)^4 * exp(-4n * acot(n)) / (-expm1(-2π * n))
"Γ(np → 1s)/ω_LyLimit (Kokubo 2024 footnote to his eq. for the Peebles formula)."
γ1(n) = 2.0^8 / 9 * ALPHA^3 * n * ((n - 1) / (n + 1))^(2n) / ((n - 1)^2 * (n + 1)^2)
# Raman coefficients (bound: c(n), continuum: c′(n′)) for final states 2s, 3s, 4s, 3d, 4d
const NB = (2, 3, 4, 3, 4)
cb(k, n) = k == 1 ? 2.0^9.5 / 3 * n^3 / (n^2 - 1) / (n^2 - 4)^2 * ((n - 1) / (n + 1))^n * ((n - 2) / (n + 2))^n :
           k == 2 ? 2.0^7 * sqrt(3) * n^3 / (n^2 - 1) / (n^2 - 9)^3 * (7n^2 - 27) * ((n - 1) / (n + 1))^n * ((n - 3) / (n + 3))^n :
           k == 3 ? 2.0^12 / 9 * n^3 / (n^2 - 1) / (n^2 - 16)^4 * (23n^4 - 288n^2 + 768) * ((n - 1) / (n + 1))^n * ((n - 4) / (n + 4))^n :
           k == 4 ? 2.0^9 * sqrt(3 / 5) * n^5 / (n^2 - 1) / (n^2 - 9)^3 * ((n - 1) / (n + 1))^n * ((n - 3) / (n + 3))^n :
                    2.0^13.5 / 9 / sqrt(5) * n^5 / (n^2 - 1) / (n^2 - 16)^4 * (7n^2 - 48) * ((n - 1) / (n + 1))^n * ((n - 4) / (n + 4))^n
function cc(k, n)
    m = NB[k]; e = exp(-2n * (acot(n) + acot(n / m))) / (-expm1(-2π * n))
    k == 1 ? 2.0^9.5 / 3 * n^3 / (n^2 + 1) / (n^2 + 4)^2 * e :
    k == 2 ? 2.0^7 * sqrt(3) * n^3 / (n^2 + 1) / (n^2 + 9)^3 * (7n^2 + 27) * e :
    k == 3 ? 2.0^12 / 9 * n^3 / (n^2 + 1) / (n^2 + 16)^4 * (23n^4 + 288n^2 + 768) * e :
    k == 4 ? 2.0^9 * sqrt(3 / 5) * n^5 / (n^2 + 1) / (n^2 + 9)^3 * e :
             2.0^13.5 / 9 / sqrt(5) * n^5 / (n^2 + 1) / (n^2 + 16)^4 * (7n^2 + 48) * e
end
const WK = Ref{Any}(nothing)
function work()
    WK[] === nothing || return WK[]
    nq, wq = cquad()
    ns = collect(2.0:NMAX)
    WK[] = (; nq, wq, bq = 1 .+ 1 ./ nq .^ 2, gq = gc.(nq) .* wq, ns, an = 1 .- 1 ./ ns .^ 2, gn = gb.(ns), γn = γ1.(ns),
              cbn = [[n == NB[k] ? 0.0 : cb(k, n) for n in ns] for k in 1:5], ccq = [cc.(k, nq) .* wq for k in 1:5])
end
"σ_1s/σ_T (Rayleigh) at ε = hν/χ_H < 1."
function ray(ε)
    W = work(); A = 0.0 + 0.0im
    @inbounds for i in eachindex(W.ns)
        A += W.gn[i] * (ε / (W.an[i] - ε - 0.5im * W.γn[i]) - ε / (W.an[i] + ε))
    end
    # bound states above NMAX: g_n → g_∞ n⁻³ with a_n → 1 (g_∞ = 2⁷3⁻¹ e⁻⁴)
    A += 2.0^7 / 3 * exp(-4.0) / (2 * (NMAX + 0.5)^2) * (ε / (1 - ε) - ε / (1 + ε))
    @inbounds for i in eachindex(W.nq)
        A += W.gq[i] * (ε / (W.bq[i] - ε) - ε / (W.bq[i] + ε))
    end
    abs2(A)
end
"σ/σ_T of the Raman branch k (1: 2s, 2: 3s, 3: 4s, 4: 3d, 5: 4d) at ε."
function raman(k, ε)
    m = NB[k]; th = 1 - 1 / m^2; ε > th || return 0.0
    W = work(); A = 0.0 + 0.0im; c = W.cbn[k]; q = W.ccq[k]
    @inbounds for i in eachindex(W.ns)
        c[i] == 0 && continue
        A += c[i] * (1 / (W.an[i] - ε - 0.5im * W.γn[i]) + 1 / (1 / m^2 - 1 / W.ns[i]^2 + ε))
    end
    @inbounds for i in eachindex(W.nq)
        A += q[i] * (1 / (W.bq[i] - ε) + 1 / (1 / m^2 + W.bq[i] - 1 + ε))
    end
    (1 - th / ε) * abs2(A)
end
"Total Raman σ/σ_T (1s → 2s, 3s, 3d, 4s, 4d)."
ramtot(ε) = raman(1, ε) + raman(2, ε) + raman(3, ε) + raman(4, ε) + raman(5, ε)

const TAB = Ref{Any}(nothing)
"Table on ε = EPS0 … EPS1 (NE points, uniform): σ_1s and σ_Raman [cm²] at the mid-points of NE cells."
function table()
    TAB[] === nothing || return TAB[]
    work(); dε = (EPS1 - EPS0) / NE; ε = [EPS0 + (i - 0.5) * dε for i in 1:NE]
    s1 = zeros(NE); sr = zeros(NE)
    Threads.@threads for i in 1:NE; s1[i] = σT * ray(ε[i]); sr[i] = σT * ramtot(ε[i]); end
    TAB[] = (; ε, dε, s1, sr)
end
"Checks of the module against exact limits; returns (TRK sum, σ_1s/(σ_T ε⁴) at ε = 1e-3)."
function hscatter_check()
    W = work()
    trk = 2 * sum(W.gn) + 2 * 2.0^7 / 3 * exp(-4.0) / (2 * (NMAX + 0.5)^2) + 2 * sum(W.gq)
    (trk, ray(1e-3) / 1e-12)
end

end # module
