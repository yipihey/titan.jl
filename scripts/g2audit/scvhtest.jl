include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
# SCVH x_e on and between isotherms vs the network's Saha, at fixed log ρ
for lρ in (-8.0, -6.42, -4.0, -2.0, -1.0, 0.0)
    println("log rho = $lρ")
    for lT in 3.86:0.02:4.42
        ρ = 10^lρ; T = 10^lT
        s = try SCVH.state(ρ, T) catch; nothing end
        s === nothing && continue
        X2, X1 = s.xH2, s.xH; xp = max(1 - X2 - X1, 0.0) / (1 + X2)
        n = net_lte(ρ, T); ys, z1, z2, xe = EQ.saha_he(ρ, T)
        node = abs((lT - 2.10) / 0.08 - round((lT - 2.10) / 0.08)) < 1e-6
        u1, p1 = EQ.eos_he(ρ, T, ys, z1, z2)
        @printf("  lT %.2f %s xHp SCVH %.3e Saha(tab) %.3e net %.3e | p SCVH/Saha %.4f  u+UOFF SCVH/Saha %.4f\n", lT, node ? "*" : " ", xp, ys[1], n.nHp / n.NH, s.p / p1, s.u / (u1 + EQ.UOFF))
    end
end
