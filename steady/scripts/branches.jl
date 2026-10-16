# Map ALL steady solutions R(Ṁ): fine scan of the shooting sign over log R, bisect every sign change.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Printf, Serialization
tab = load_table(); M = 1e6Msun; L = LEdd_es(M)
Mdots = [0.0, 0.3, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 30.0, 100.0]
lRs = collect(14.0:0.025:16.6)
out = open(joinpath(@__DIR__, "..", "output", "branches.csv"), "w")
println(out, "fBH,Mdot,R,Teff_L,Tb,Pgas_b,ram,tau_LTE,r_ph_infall,T_obs,Gamma_kappa_b")
allscans = Dict()
for fBH in (0.1, 0.01)
    p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
    for md in Mdots
        Mdot = md * Msun / yr
        sg = [shoot_sign(10^lR, p; Mdot = Mdot)[1] for lR in lRs]
        allscans[(fBH, md)] = sg
        roots = Float64[]
        for i in 1:length(lRs)-1
            (sg[i] == 0 || sg[i+1] == 0 || sg[i] == sg[i+1]) && continue
            a, b, sa = lRs[i], lRs[i+1], sg[i]
            for _ in 1:40
                mid = 0.5(a + b); s = shoot_sign(10^mid, p; Mdot = Mdot)[1]
                s == 0 && (b = mid; continue)
                s == sa ? (a = mid) : (b = mid)
            end
            push!(roots, 0.5(a + b))
        end
        nfail = count(==(0), sg)
        @printf("fBH=%.2f Ṁ=%6.2f: %d root(s) %s   [%d/%d scan points had no BC/integration]\n", fBH, md, length(roots),
                join([@sprintf("%.3f", r) for r in roots], ", "), nfail, length(lRs))
        for lR in roots
            R = 10^lR; bc = outer_bc(R, p; Mdot = Mdot)
            Ltot = bc.Ltot; rph = bc.r_ph_infall
            Tobs = (Ltot / (4π * σSB * rph^2))^0.25
            @printf(out, "%g,%g,%.5e,%.1f,%.1f,%.4e,%.4e,%.4e,%.4e,%.1f,%.4f\n", fBH, md, R, bc.Teff_L, bc.Tb, bc.Pgas_b, bc.ram, bc.τ_in, rph, Tobs, bc.Γκ)
            @printf("      R=%.3e Tb=%6.0f Pgas_b=%.3g ram=%.3g τ_LTE=%.2f r_ph/R=%.2f T_obs=%6.0f Γκ_b=%.3f\n", R, bc.Tb, bc.Pgas_b, bc.ram, bc.τ_in, rph / R, Tobs, bc.Γκ)
        end
        flush(stdout); flush(out)
    end
end
close(out)
serialize(joinpath(@__DIR__, "..", "output", "branch_scans.jls"), (lRs, allscans))
