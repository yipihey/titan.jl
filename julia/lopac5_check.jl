# Size of the opacity derivative TITAN drops in lopac = 5 (bound-free term: the
# derivative of E1(b) is "ignored" in opac.f).  Analytic d ln κ / d ln T from
# Fortran opac vs a central FD of the same Fortran opac over T = 1e4..1e7 K.
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); using .TitanDriver; const TD = TitanDriver; const TF = TD.TF
TF.init(8; rundir = joinpath(TF.ROOT, "lib", "p3_runs", "lopac5"))
TF.setscalar!(:lopac, Int32(5))
a = Int(TF.getscalar(:ngrs)) - 1; b = Int(TF.getscalar(:ngre)) + 2
tn, dn = TF.var(:tn), TF.var(:dn)
Ts = exp10.(range(4, 7, length = b - a + 1))
for (i, k) in enumerate(a:b); tn[k] = Ts[i]; dn[k] = 1e-9; end
function op(); TD.fcall(:opac_f_); (copy(TF.var(:chifn)[a:b]), copy(TF.var(:xken)[a:b]), copy(TF.var(:xkpn)[a:b]),
     copy(TF.var(:dlcfdltn)[a:b]), copy(TF.var(:dlkedltn)[a:b]), copy(TF.var(:dlkpdltn)[a:b])); end
cf, ke, kp, dcf, dke, dkp = op()
h = 1e-6; T0 = copy(tn)
tn[a:b] .= T0[a:b] .* exp(h); cfp, kep, kpp = op()
tn[a:b] .= T0[a:b] .* exp(-h); cfm, kem, kpm = op()
fd(p, m, v) = (log.(p) .- log.(m)) ./ (2h)
println("lopac=5 at ρ = 1e-9 g/cm³: analytic vs FD d ln κ/d ln T (κ_E includes Compton term with E_r of the state)")
for (lbl, an, f) in (("χ_F (flux mean)", dcf, fd(cfp, cfm, 0)), ("κ_E", dke, fd(kep, kem, 0)), ("κ_P (Planck)", dkp, fd(kpp, kpm, 0)))
    err = abs.(an .- f)
    i = argmax(err)
    @printf("  %-16s max |analytic − FD| = %.3f at T = %.3g K (analytic %.3f, FD %.3f); median %.3f\n",
            lbl, err[i], Ts[i], an[i], f[i], sort(err)[end ÷ 2])
end
for T in (1e4, 3e4, 1e5, 3e5, 1e6, 1e7)
    i = argmin(abs.(Ts .- T))
    @printf("     T=%.1e  κ_P: analytic %+.3f  FD %+.3f\n", Ts[i], dkp[i], fd(kpp, kpm, 0)[i])
end
TF.finalize()
