# Steady (t_KH ≪ M/Ṁ) quasi-star envelopes with a spherical free-fall accretion shock as outer boundary.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Serialization, Printf
tab = load_table()
const kB = 1.380649e-16; const mH = 1.6735575e-24; const eV = 1.602176634e-12
M = 1e6Msun; L = LEdd_es(M)
Mdots = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]   # Msun/yr (= Ṁ/f for streams)
rows = String[]; profs = Dict{Tuple{Float64,Float64},Any}()
push!(rows, "fBH,Mdot_Msun_yr,R_cm,logg,Teff_L,Tb,Pb,Pgas_b,ram,tau_in_LTE,tau_in_es,Lacc_over_L,Gamma_es_total,Gamma_phot,maxMach,T_at_maxMach,maxGamma,Tshock,lcool_over_R,ion_ratio_35eV,ion_ratio_1keV")
for fBH in (0.1, 0.01)
    p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
    for md in Mdots
        Mdot = md * Msun / yr
        s = solve_R(p; Mdot = Mdot, lRmin = 15.0, lRmax = 16.8, nscan = 37)
        if s === nothing
            println("fBH=$fBH Mdot=$md: NO SOLUTION"); push!(rows, "$fBH,$md" * ",NaN"^19); continue
        end
        R = s.R; bc = s.bc; pr = s.prof
        mach = pr[:v] ./ pr[:cs]; im = argmax(mach)
        # diagnostics of the infall
        v = vff(M, R); ρin = Mdot > 0 ? rho_infall(Mdot, M, R) : 0.0
        τes = 0.352 * ρin * 2R                                   # ∫ κ_es ρ dr, ρ ∝ r^-3/2
        Lacc = G * M * Mdot / R
        Γes = (L + Lacc) / LEdd_es(M)                            # radiation force / gravity on ionized gas
        Tsh = 3 / 16 * 0.59 * mH * v^2 / kB
        # post-shock cooling (bremsstrahlung, Z=0): t = 3 n k T / (n_e n_i Λ), compression ×4
        nps = 4ρin / (0.59mH); Λ = 2.4e-27 * sqrt(Tsh)
        lcool = Mdot > 0 ? (v / 4) * 3kB * Tsh / (nps * Λ) / R : 0.0
        # photon budget to keep infall ionized out to 10 R: recombinations vs ionizations from half of L_acc
        nR = 0.76ρin / mH; αB = 2.6e-13
        Nrec = αB * nR^2 * 4π * R^3 * log(10.0)
        ion35 = Mdot > 0 ? 0.5Lacc / (35eV) / Nrec : NaN
        ion1k = Mdot > 0 ? 0.5Lacc / (1000eV) / Nrec : NaN
        push!(rows, @sprintf("%g,%g,%.5e,%.4f,%.1f,%.1f,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.5f,%.4f,%.4f,%.4e,%.3f,%.3e,%.3e,%.3e,%.3e",
              fBH, md, R, log10(G * M / R^2), bc.Teff_L, bc.Tb, bc.Pb, bc.Pgas_b, bc.ram, bc.τ_in, τes, Lacc / L, Γes,
              bc.Γκ, mach[im], pr[:T][im], maximum(pr[:Γ]), Tsh, lcool, ion35, ion1k))
        profs[(fBH, md)] = (R = R, bc = bc, prof = pr)
        @printf("fBH=%.2f Ṁ=%6.2f: R=%.3e Teff=%6.0f Tb=%6.0f Pgas_b=%.3g ram=%.3g τ_LTE=%.2e τ_es=%.2e Γ_es=%.4f Mach=%.2f maxΓ=%.1f lcool/R=%.1e ion(35eV)=%.2e ion(1keV)=%.2e\n",
                fBH, md, R, bc.Teff_L, bc.Tb, bc.Pgas_b, bc.ram, bc.τ_in, τes, Γes, mach[im], maximum(pr[:Γ]), lcool, ion35, ion1k)
        flush(stdout)
    end
end
open(joinpath(@__DIR__, "..", "output", "infall_grid.csv"), "w") do io; foreach(r -> println(io, r), rows); end
serialize(joinpath(@__DIR__, "..", "output", "infall_profiles.jls"), profs)
