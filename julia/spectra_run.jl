# Emergent spectra of TITAN snapshots (Spectrum.jl).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env spectra_run.jl RUN TIME[s|end] [ext=1] [rmax=100] [n=500] [lucy=0] [rmin=0.6] [out=NAME]
# RUN: lib/qs_runs/RUN_relax_snaps.jls; TIME: snapshot nearest to that time. ext=1 continues the infall as free fall
# beyond the grid to rmax·r_out. Writes lib/spectra/NAME.jls with ν, L_ν and the atmosphere, and prints the
# frequency-integrated L against TITAN's L_rad at the outer boundary.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics
run = ARGS[1]; tq = ARGS[2]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[3:end])
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(run)_relax_snaps.jls"))
s = tq == "end" ? S[end] : S[argmin([abs(x.t - parse(Float64, tq)) for x in S])]
M = 1e6 * Z0.Msun
Mdot = -4π * s.r[end]^2 * s.d[end] * s.u[end]
ext = get(opt, :ext, "1") == "1" ? (Mdot, M, parse(Float64, get(opt, :rmax, "100")), 1000.0) : nothing
atm = SP.atmosphere(s; ext, rmin = parse(Float64, get(opt, :rmin, "0.6")) * 8.4352e15, split = parse(Float64, get(opt, :split, "1.5")))
# ion=1: the infall (every zone outside the shock, and the extension) is taken as fully ionised instead of in LTE at
# its local T (the opposite limit to LTE for gas irradiated by the star and the shock; see shock_budget.jl)
T_titan = copy(atm.T)
ν = SP.grid(; n = parse(Int, get(opt, :n, "500")))
t0 = time()
# Unsöld–Lucy correction of the hydrostatic atmosphere below the shock (lucy=N iterations): TITAN's temperatures are in
# radiative equilibrium for its grey (Planck/Rosseland-mean) radiation field, not for the non-grey LTE opacities used
# here; without the correction the frequency-integrated flux comes out up to 24 % too high (UV window, 912–3646 Å).
nlucy = parse(Int, get(opt, :lucy, "0"))
get_lucytop(opt) = parse(Int, get(opt, :ltop, "2"))     # Lucy up to zone kc − ltop (2: the zones right under the hot layer are irradiated by it and the correction oscillates there)
rmin = parse(Float64, get(opt, :rmin, "0.6")) * 8.4352e15
keep = findall(k -> 0.5 * (s.r[k] + s.r[k+1]) >= rmin, eachindex(s.d))
du = [s.u[k+1] - s.u[k] for k in 1:length(s.u)-1]; ksh = argmin(du)          # shock zone (largest velocity drop)
# no shock: no zone drops the velocity by 100 km/s, and the 21 zones around the largest drop by less than 500 km/s together
# (the 300-zone wide-domain grid, W30, spreads the jump over ~10 very fine zones, each dropping by less than 100 km/s)
noshock = du[ksh] > -1e7 && sum(du[max(1, ksh - 10):min(end, ksh + 10)]) > -5e7
noshock && (ksh = length(s.d))
# bottom of the hot post-shock cooling layer: the last zone below the shock with T < 1e4 K. The cooling layer is
# heated by the shock, not in radiative equilibrium, so it is excluded from the Unsöld–Lucy correction.
kc = noshock ? length(s.d) : something(findlast(k -> k < ksh && s.T[k] < 1e4, 1:length(s.d)), ksh)
Lz = [0.5 * (s.Lrad[k] + s.Lconv[k] + s.Lrad[k+1] + s.Lconv[k+1]) for k in eachindex(s.d)]   # TITAN L at zone centres
Ltarget = [z > 0 ? Lz[z] : s.Lrad[end] for z in atm.zone]
# hotion=1: the post-shock cooling layer (and any shock-heated zone above it) at T > 2e4 K is treated as fully
# ionised plasma (free–free and Thomson only) instead of in LTE (a test of how much its bound–free opacity matters).
# bare=1: drop the infall (every point outside the shock zone) to isolate what the infall itself adds
if get(opt, :bare, "0") == "1"
    kp = findall(z -> 0 < z <= ksh, atm.zone)
    atm = (r = atm.r[kp], ρ = atm.ρ[kp], T = atm.T[kp], zone = atm.zone[kp], jump = atm.jump[kp[1:end-1]])
    T_titan = T_titan[kp]
    Ltarget = Ltarget[kp]
end
# tfrom=NAME: take the (Unsöld–Lucy corrected) temperatures of an earlier spectrum lib/spectra/NAME.jls at the same
# points (use with lucy=0, e.g. bare=1 tfrom=A220_late_ext: the corrected star of that run with the infall removed)
if haskey(opt, :tfrom)
    dT = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$(opt[:tfrom]).jls")).atm
    Tn = map(eachindex(atm.r)) do j
        i = findfirst(x -> abs(x / atm.r[j] - 1) < 1e-12, dT.r); i === nothing ? atm.T[j] : dT.T[i]
    end
    println("temperatures from $(opt[:tfrom]): ", count(j -> Tn[j] != atm.T[j], eachindex(Tn)), " points changed")
    atm = merge(atm, (T = Tn,))
end
hot = [z > kc && T_titan[j] > 2e4 && get(opt, :hotion, "0") == "1" for (j, z) in enumerate(atm.zone)]
infall_ion = get(opt, :ion, "0") == "1"
atm = merge(atm, (ion = [hot[j] || (infall_ion && (z == 0 || z > ksh)) for (j, z) in enumerate(atm.zone)],))
println("zones treated as fully ionised: ", count(atm.ion), " (hot layer ", count(hot), ")")
# corrected zones: from 0.95 R up to zone kc − ltop, just below the bottom of the cooling layer (all zones above 0.95 R if no shock)
# ldr=x (in R, default 0): also leave out the zones within x·R below the bottom of the cooling layer. On the 300-zone
# wide-domain grid (W30) the bottom of the layer is resolved by zones of 1e-8 R, so two zones are a far thinner slice
# than on A220 (5e-6 R) and the correction there hits its ±30 % clamp; ldr = 5e-6 leaves out the same slice as ltop = 2 on A220.
ldr = parse(Float64, get(opt, :ldr, "0")) * 8.4352e15
sel = findall(j -> 0 < atm.zone[j] <= kc - get_lucytop(opt) && atm.r[j] >= 0.95 * 8.4352e15 && atm.r[j] <= s.r[kc+1] - ldr, eachindex(atm.r))
# re=1 (experimental, does not converge for the A220 snapshot: the 1000 K extension is far from equilibrium with the
# stellar light and the 30 % per-pass clamp binds for all 20 passes): the infall is brought to local radiative equilibrium with the
# non-grey field instead of keeping TITAN's (grey) temperatures and the assumed T ∝ r^-1/2 of the extension
thin = get(opt, :re, "0") == "1" ? findall(z -> z == 0 || z > ksh, atm.zone) : Int[]
println("shock zone $ksh (r/R = $(round(s.r[ksh+1]/8.4352e15, digits = 5))), cooling layer down to zone $kc; Lucy on ", length(sel), " zones")
Tlog = Float64[]
# the Lucy passes use a coarser frequency grid (nl=…, default 200); the final spectrum (and the flux check) uses n
νl = SP.grid(; n = parse(Int, get(opt, :nl, "200")))
for it in 1:nlucy
    _, _, mom = SP.spectrum(atm, νl; moments = true)
    Tn, err = SP.lucy(atm, mom, sel, Ltarget; damp = parse(Float64, get(opt, :damp, "1")), thin)
    rerr = isempty(thin) ? 0.0 : maximum(abs(mom.κJ[k] * mom.J[k] / (mom.κB[k] * mom.B[k]) - 1) for k in thin)
    jmax = argmax(abs.(Tn ./ atm.T .- 1))
    @printf("Lucy %2d: max |ΔH/H| in the corrected zones %.4f, max |ΔT/T| %.4f (zone %d, r/R %.5f, T %.0f K); infall max |κJ J/κB B − 1| %.4f\n",
            it, err, maximum(abs.(Tn ./ atm.T .- 1)), atm.zone[jmax], atm.r[jmax] / 8.4352e15, atm.T[jmax], rerr)
    global atm = merge(atm, (T = Tn,))
end
L, its = SP.spectrum(atm, ν; verbose = true)
Ltot = sum(0.5 * (L[i] + L[i+1]) * (ν[i+1] - ν[i]) for i in 1:length(ν)-1)
@printf("%s t=%.4e: %d points (%d in extension), %d frequencies in %.0f s; ∫L_ν dν = %.4e erg/s, TITAN L_rad(out) = %.4e (ratio %.4f), Ṁ = %.3f Msun/yr\n",
        run, s.t, length(atm.r), count(==(0), atm.zone), length(ν), time() - t0, Ltot, s.Lrad[end], Ltot / s.Lrad[end], Mdot * Z0.yr / Z0.Msun)
mkpath(joinpath(@__DIR__, "..", "lib", "spectra"))
name = get(opt, :out, "$(run)_$(round(Int, s.t))_ext$(ext === nothing ? 0 : 1)")
serialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls"), (run = run, t = s.t, ν = ν, L = L, its = its, atm = atm, T_titan = T_titan, nlucy = nlucy, Mdot = Mdot, Ltitan = s.Lrad[end], ksh = ksh, kc = kc, keep = keep, ion = get(opt, :ion, "0") == "1"))
println("saved lib/spectra/$name.jls")
