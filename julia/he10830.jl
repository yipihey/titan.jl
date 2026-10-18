# He I λ10830 in the infall of a non-LTE solution: the population of the metastable 2³S level and the line profile.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env he10830.jl NAME [s3=1] [xsec=1] [out=TAG]
# Three-state helium (1¹S, 2³S, He⁺) after Oklopčić & Hirata (2018, ApJL 855, L11), in the radiation field of the full
# solution (J_ν at every non-LTE point from nlte_rates.jl, lib/spectra/nlterates_NAME.jls) at the solution's
# temperatures and electron and neutral-hydrogen densities:
#   He⁺:  n₊ n_e (α₁ + α₃) = n₁ Γ₁ + n₃ Γ₃
#   2³S:  n₊ n_e α₃ + n₁ n_e q₁₃ = n₃ (Γ₃ + A₃₁ + n_e (q₃₁ₐ + q₃₁ᵦ) + n_HI Q₃₁)
#   n₁ + n₃ + n₊ = n_He = (Y/4X) n_H
# Γ₁ = 4π∫J σ₁/hν (24.6 eV and up, X-rays included; σ₁ from Spectrum.σHe1), Γ₃ = 4π∫J σ₃/hν (above 4.77 eV, λ < 2600 Å).
# Rates (approximate, as stated): α₁ = 1.54e-13 (T/10⁴)^−0.486, α₃ = 2.10e-13 (T/10⁴)^−0.778 cm³ s⁻¹ (case B fits to
# Osterbrock & Ferland 2006); q₃₁ₐ = 2.6e-8, q₃₁ᵦ = 4.0e-9 (2³S → 2¹S, 2¹P) and q₁₃ = 4.5e-20 cm³ s⁻¹ as quoted by
# Oklopčić & Hirata at 10⁴ K, scaled with T as their Eq. 16 (T^−½ e^{−ΔE/kT}, constant collision strengths);
# Q₃₁ = 5e-10 cm³ s⁻¹ (neutral H); A₃₁ = 1.272e-4 s⁻¹. σ₃: 5.5e-18 cm² at threshold, ∝ ν^−2 (an approximation to
# Norcross 1971; xsec=f scales it, s3=f scales α₃, both sensitivity tests). Secondary ionizations of He by the X-rays
# are neglected.
# The line (2³S → 2³P, the three components treated as one at 10830.3 Å, f = 0.539 from A = 1.0216e7 s⁻¹): Sobolev
# depths along rays as in nlte_spectrum.jl, source function = the local continuum J at 10830 Å (pure scattering: 2³P
# decays only to 2³S). Writes lib/spectra/he10830_NAME_TAG.jls.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
s3f = parse(Float64, get(opt, :s3, "1")); xsf = parse(Float64, get(opt, :xsec, "1")); tag = get(opt, :out, "s3$(s3f)_x$(xsf)")
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"))
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
jr = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlterates_$name.jls"))
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(d.run)_relax_snaps.jls"))
s = S[argmin([abs(x.t - d.t) for x in S])]
atm = haskey(dn, :Tstar) ? merge(d.atm, (T = dn.Tstar,)) : d.atm
c, h, kB, eV = Z0.c, Z0.h, Z0.kB, Z0.eV; R0 = 8.4352e15
ν = jr.ν; wν = [i == 1 ? 0.5(ν[2] - ν[1]) : i == length(ν) ? 0.5(ν[end] - ν[end-1]) : 0.5(ν[i+1] - ν[i-1]) for i in eachindex(ν)]
X, Y = 0.76, 0.24
σ3(nu) = (E = h * nu / eV; E < 4.768 ? 0.0 : xsf * 5.5e-18 * (E / 4.768)^(-2))
σ1(nu) = (E = h * nu / eV; E < 24.587 ? 0.0 : SP.σHe1(E))
qscale(q, ΔE, T) = q * sqrt(1e4 / T) * exp(-ΔE * eV / kB * (1 / T - 1e-4))
res = map(eachindex(dn.nl)) do qd
    st = dn.nl[qd]; T = dn.T[qd]; J = jr.J[qd, :]
    nH = st.np + sum(st.n); nHI = sum(st.n); ne = st.ne; nHe = Y / (4X) * nH
    Γ1 = sum(4π * wν[k] * σ1(ν[k]) * J[k] / (h * ν[k]) for k in eachindex(ν))
    Γ3 = sum(4π * wν[k] * σ3(ν[k]) * J[k] / (h * ν[k]) for k in eachindex(ν))
    α1 = 1.54e-13 * (T / 1e4)^(-0.486); α3 = s3f * 2.10e-13 * (T / 1e4)^(-0.778)
    q31 = qscale(2.6e-8, 0.796, T) + qscale(4.0e-9, 1.398, T); q13 = qscale(4.5e-20, 19.82, T)
    out3 = Γ3 + 1.272e-4 + ne * q31 + nHI * 5e-10
    # linear system for (n1, n3, n+): rows = He⁺ balance, 2³S balance, number conservation
    A = [Γ1 Γ3 -ne*(α1+α3); ne*q13 -out3 ne*α3; 1.0 1.0 1.0]
    sol = A \ [0.0, 0.0, nHe]
    (r = dn.r[qd], T = T, n1 = sol[1], n3 = max(sol[2], 0.0), np = sol[3], nHe = nHe, Γ1 = Γ1, Γ3 = Γ3, ne = ne, out3 = out3)
end
ν10 = c / 10830.3e-8; λ10 = 10830.3e-8
f10 = 1.4992 * λ10^2 * 3 * 1.0216e7      # f_lu = (m_e c λ²/8π²e²)(g_u/g_l) A
J10 = [ (k = clamp(searchsortedlast(ν, ν10), 1, length(ν) - 1); jr.J[q, k] + (jr.J[q, k+1] - jr.J[q, k]) * (ν10 - ν[k]) / (ν[k+1] - ν[k])) for q in eachindex(dn.nl)]
K10 = [NH.πe2mc * f10 * λ10 * x.n3 for x in res]
@printf("%s: He I 2³S (f = %.3f for λ10830)\n", name, f10)
@printf("%8s %6s %10s %10s %10s %10s %10s %10s\n", "r/R0", "T", "x(He⁺)", "n(2³S)/nHe", "Γ₁ [s⁻¹]", "Γ₃ [s⁻¹]", "τ_Sob", "n(2³S)")
for xr in (1.0, 1.2, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0)
    q = argmin(abs.(dn.r ./ R0 .- xr)); x = res[q]; vff = dn.vabs[q]
    @printf("%8.2f %6.0f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", x.r / R0, x.T, x.np / x.nHe, x.n3 / x.nHe, x.Γ1, x.Γ3, K10[q] / (vff / x.r), x.n3)
end
N3 = sum(0.5 * (res[i].n3 + res[i+1].n3) * (res[i+1].r - res[i].r) for i in 1:length(res)-1)
@printf("radial column of 2³S: %.3e cm⁻²\n", N3)
# profile: ray integration (as nlte_spectrum.jl), K = πe²/mc f λ n₃, S = J_c(10830)
rsh = s.r[d.ksh+1]; inf = dn.sel; rI = atm.r[inf]; TI = atm.T[inf]; rmax = rI[end]; uabs = dn.vabs
function interp(x, xs, ys)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (log(x) - log(xs[i])) / (log(xs[i+1]) - log(xs[i])); ys[i] + t * (ys[i+1] - ys[i])
end
nlfull = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(inf); nlfull[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nlfull)
_, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
_, _, _, ci = SP.monochromatic(atm, op, R, ν10; moments = :I); ps, Ic = ci.p, ci.I
Fc = sum(0.5 * (Ic[j] * ps[j] + Ic[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1)
vgrid = collect(range(-2500e5, 2500e5, length = 251)); nz = 2000; gq = [(exp(5 * (i - 1) / (nz - 1)) - 1) / (exp(5) - 1) for i in 1:nz]
mHe = 4.0026 * 1.6605e-24
F = zeros(length(vgrid))
Threads.@threads for iv in eachindex(vgrid)
    v = vgrid[iv]; Iv = similar(Ic)
    for j in eachindex(ps)
        p = ps[j]; I = Ic[j]
        znear = sqrt(max(rmax^2 - p^2, 0.0)); zfar = p < rsh ? sqrt(rsh^2 - p^2) : -znear
        if znear > zfar && p < rmax
            zs = p < rsh ? zfar .+ (znear - zfar) .* gq : (zp = znear .* gq; vcat(-reverse(zp), zp[2:end]))
            for i in 1:length(zs)-1
                r1, r2 = sqrt(p^2 + zs[i]^2), sqrt(p^2 + zs[i+1]^2)
                (r1 < rI[1] || r2 < rI[1]) && continue
                w1, w2 = interp(r1, rI, uabs) * zs[i] / r1, interp(r2, rI, uabs) * zs[i+1] / r2
                D1, D2 = sqrt(2kB * interp(r1, rI, TI) / mHe + 1e12), sqrt(2kB * interp(r2, rI, TI) / mHe + 1e12)
                lo, hi = minmax(w1, w2); (v < lo - 5max(D1, D2) || v > hi + 5max(D1, D2)) && continue
                K1, K2 = interp(r1, rI, K10), interp(r2, rI, K10); S1, S2 = interp(r1, rI, J10), interp(r2, rI, J10)
                m = max(1, ceil(Int, (hi - lo) / (0.25 * min(D1, D2)))); dz = (zs[i+1] - zs[i]) / m
                for qq in 1:m
                    t = (qq - 0.5) / m; wt = w1 + t * (w2 - w1); Dt = D1 + t * (D2 - D1); x = (v - wt) / Dt
                    abs(x) > 6 && continue
                    τ = (K1 + t * (K2 - K1)) * exp(-x^2) / (sqrt(π) * Dt) * dz
                    I = I * exp(-τ) + (S1 + t * (S2 - S1)) * (-expm1(-τ))
                end
            end
        end
        Iv[j] = I
    end
    F[iv] = sum(0.5 * (Iv[j] * ps[j] + Iv[j+1] * ps[j+1]) * (ps[j+1] - ps[j]) for j in 1:length(ps)-1) / Fc
end
EW = sum(0.5 * ((1 - F[i]) + (1 - F[i+1])) * (vgrid[i+1] - vgrid[i]) for i in 1:length(vgrid)-1) / c * 10830.3
k = argmin(F); kx = argmax(F)
@printf("λ10830 profile: minimum %.4f at %+.0f km/s, maximum %.4f at %+.0f km/s, EW %+.3f Å (positive = absorption)\n",
        F[k], vgrid[k] / 1e5, F[kx], vgrid[kx] / 1e5, EW)
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "he10830_$(name)_$tag.jls"), (res = res, v = vgrid ./ 1e5, F = F, EW = EW, N3 = N3, f10 = f10))
println("saved lib/spectra/he10830_$(name)_$tag.jls")
