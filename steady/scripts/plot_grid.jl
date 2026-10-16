using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie, CSV, DataFrames, Serialization
const OUT = joinpath(@__DIR__, "..", "output")
df = CSV.read(joinpath(OUT, "infall_grid.csv"), DataFrame)
profs = deserialize(joinpath(OUT, "infall_profiles.jls"))

blue, orange = colorant"#2a78d6", colorant"#eb6834"
ink2 = colorant"#52514e"
set_theme!(Theme(fontsize = 14, Axis = (xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
               xminorticksvisible = false, topspinevisible = false, rightspinevisible = false)))

sub(f) = sort(df[(df.fBH .== f) .& (df.Mdot_Msun_yr .> 0) .& (df.Mdot_Msun_yr .<= 1.0), :], :Mdot_Msun_yr)   # extended branch only (first-root grid; Ṁ>1 rows are the compact branch)
s1, s2 = sub(0.1), sub(0.01)
base(f) = df[(df.fBH .== f) .& (df.Mdot_Msun_yr .== 0), :]
b1, b2 = base(0.1), base(0.01)

fig = Figure(size = (1100, 820))
lab(f) = "M_BH/M = $(f)"
# (a) Teff of the star (no-infall reference as dashed)
ax = Axis(fig[1, 1], xscale = log10, xlabel = "Ṁ (M⊙/yr)  [= Ṁ/f for a stream covering fraction f]",
          ylabel = "T at stellar boundary (K)", title = "a  Boundary temperature (extended branch)")
for (s, b, col, f) in ((s1, b1, blue, 0.1), (s2, b2, orange, 0.01))
    scatterlines!(ax, s.Mdot_Msun_yr, s.Tb, color = col, linewidth = 2, markersize = 9, label = lab(f))
    hlines!(ax, b.Tb, color = col, linestyle = :dash, linewidth = 1)
end
axislegend(ax, position = :lt, framevisible = false)

# (b) ram pressure vs the photospheric gas pressure of the no-infall star
ax = Axis(fig[1, 2], xscale = log10, yscale = log10, xlabel = "Ṁ (M⊙/yr)",
          ylabel = "pressure (dyn cm⁻²)", title = "b  Ram pressure vs. no-infall photospheric P_gas")
lines!(ax, s1.Mdot_Msun_yr, s1.ram, color = ink2, linewidth = 2, label = "ρ_in v_ff² (both)")
hlines!(ax, b1.Pgas_b, color = blue, linestyle = :dash, linewidth = 1.5, label = "P_gas,phot  M_BH/M=0.1")
hlines!(ax, b2.Pgas_b, color = orange, linestyle = :dash, linewidth = 1.5, label = "P_gas,phot  M_BH/M=0.01")
axislegend(ax, position = :lt, framevisible = false)

# (c) optical depth of infall: LTE (self-consistent if not photoionized) vs fully ionized
ax = Axis(fig[2, 1], xscale = log10, yscale = log10, xlabel = "Ṁ (M⊙/yr)",
          ylabel = "τ of infall column", title = "c  Infall optical depth: LTE vs. fully ionized")
scatterlines!(ax, s1.Mdot_Msun_yr, s1.tau_in_LTE, color = blue, linewidth = 2, markersize = 9, label = "LTE (κ_R of neutral gas)")
lines!(ax, s1.Mdot_Msun_yr, s1.tau_in_es, color = blue, linewidth = 2, linestyle = :dot, label = "ionized (κ_es)")
hlines!(ax, [2/3], color = ink2, linewidth = 1)
text!(ax, 0.012, 0.75, text = "τ = 2/3", color = ink2, fontsize = 12)
axislegend(ax, position = :lt, framevisible = false)

# (d) photon budget: ionizations by shock radiation / recombinations in infall (≤10 R)
ax = Axis(fig[2, 2], xscale = log10, yscale = log10, xlabel = "Ṁ (M⊙/yr)",
          ylabel = "ionizations / recombinations", title = "d  Can the shock keep the infall ionized?")
band!(ax, s1.Mdot_Msun_yr, s1.ion_ratio_1keV, s1.ion_ratio_35eV, color = (blue, 0.25))
lines!(ax, s1.Mdot_Msun_yr, s1.ion_ratio_35eV, color = blue, linewidth = 2, label = "35 eV per ionization")
lines!(ax, s1.Mdot_Msun_yr, s1.ion_ratio_1keV, color = blue, linewidth = 2, linestyle = :dash, label = "1 keV per ionization")
hlines!(ax, [1.0], color = ink2, linewidth = 1)
axislegend(ax, position = :rt, framevisible = false)
save(joinpath(OUT, "fig_infall_boundary.png"), fig, px_per_unit = 2)

# ---- profiles: local Eddington factor and convective Mach number vs T, fBH = 0.1 and 0.01
fig2 = Figure(size = (1100, 430))
for (k, (f, col)) in enumerate(((0.1, blue), (0.01, orange)))
    ax = Axis(fig2[1, k], xscale = log10, yscale = log10, xlabel = "T (K)",
              ylabel = "Γ = κL/(4πcGm)", title = "$(k == 1 ? "a" : "b")  Envelope, M_BH/M = $(f), Ṁ = 0",
              xreversed = false)
    pr = profs[(f, 0.0)].prof
    lines!(ax, pr[:T], pr[:Γ], color = col, linewidth = 2, label = "Γ (radiation alone)")
    lines!(ax, pr[:T], max.(pr[:v] ./ pr[:cs], 1e-4), color = ink2, linewidth = 2, linestyle = :dash, label = "MLT Mach number")
    hlines!(ax, [1.0], color = (:black, 0.4), linewidth = 1)
    xlims!(ax, 5e3, 1e6); ylims!(ax, 1e-3, 300)
    axislegend(ax, position = :rt, framevisible = false)
end
save(joinpath(OUT, "fig_envelope_profiles.png"), fig2, px_per_unit = 2)

# ---- solution branches R(Ṁ) and temperatures (stellar boundary vs. observed infall photosphere)
br = CSV.read(joinpath(OUT, "branches.csv"), DataFrame)
br.Mplot = max.(br.Mdot, 0.08)              # put Ṁ = 0 at the left edge of the log axis
fig3 = Figure(size = (1100, 430))
axR = Axis(fig3[1, 1], xscale = log10, yscale = log10, xlabel = "Ṁ (M⊙/yr)   (Ṁ = 0 plotted at 0.08)",
           ylabel = "R of stellar boundary (cm)", title = "a  Steady solutions: two branches that merge and vanish")
axT = Axis(fig3[1, 2], xscale = log10, yscale = log10, xlabel = "Ṁ (M⊙/yr)   (Ṁ = 0 plotted at 0.08)",
           ylabel = "T (K)", title = "b  Boundary T (filled) vs. observed photosphere T (open)")
for (f, col) in ((0.1, blue), (0.01, orange))
    d = br[br.fBH .== f, :]
    # split branches: compact = R < 0.5 × R(Ṁ=0)
    R0 = d.R[d.Mdot .== 0][1]
    for (mask, ls, nm) in ((d.R .> 0.5R0, :solid, "extended"), (d.R .<= 0.5R0, :dash, "compact"))
        dd = sort(d[mask, :], :Mdot)
        nrow(dd) == 0 && continue
        scatterlines!(axR, dd.Mplot, dd.R, color = col, linestyle = ls, linewidth = 2, markersize = 8,
                      label = "M_BH/M=$(f), $(nm)")
        scatterlines!(axT, dd.Mplot, dd.Tb, color = col, linestyle = ls, linewidth = 2, markersize = 8)
        scatter!(axT, dd.Mplot, dd.T_obs, color = :white, strokecolor = col, strokewidth = 1.5, markersize = 9)
    end
end
axislegend(axR, position = :lb, framevisible = false, labelsize = 11)
save(joinpath(OUT, "fig_branches.png"), fig3, px_per_unit = 2)
println("saved")
