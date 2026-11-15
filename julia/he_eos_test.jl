using Printf, ForwardDiff
dir = ARGS[1]; mode = ARGS[2]
for f in ("PrimordialChem", "SCVH", "PopIIIEqEOS"); include(joinpath(dir, f * ".jl")); end
include("/home/tabel/Projects/titan.jl/julia/Z0Physics.jl")
const EQ = PopIIIEqEOS; const PC = PrimordialChem
function main(mode)
    pts = [(-12.0, 3.4), (-10.9, 4.6), (-8.0, 4.0), (-5.0, 3.6), (-14.0, 5.2), (-3.0, 5.5), (-11.0, 3.9)]
    if mode == "bits"
        for (lr, lT) in pts; println(join(string.(reinterpret(UInt64, collect(Float64, EQ.node(lr, lT)))), " ")); end
    else
        Z0 = Main.Z0Physics
        for (ρ, T) in ((1.2e-11, 4.0e4), (4.5e-11, 3.5e4), (2.3e-10, 1.3e5), (1e-8, 2e4), (1e-12, 9e4), (5e-9, 4.5e5))
            y, z1, z2, xe = EQ.saha_he(ρ, T); h0 = PC.yHe - z1 - z2
            s = Z0.saha(ρ, T)    # (ne, nH0, nHp, nHe0, nHe1, nHe2), Y = 0.24 composition
            ne = xe * ρ / PC.μH
            # compare the ionisation ratios at each solver's own n_e via the Saha products (n1 ne / n0)
            r1 = z1 / h0 * ne; r2 = z2 / z1 * ne; q1 = s[5] / s[4] * s[1]; q2 = s[6] / s[5] * s[1]
            u0, p0 = PC.eos_rhfc02(ρ, T, PC.y_saha(ρ, T)); u1, p1 = EQ.eos_he(ρ, T, y, z1, z2)
            @printf("rho %.1e T %.2e: He0 %.3e He+ %.3e He++ %.3e (of yHe) xe %.4f | (n1 ne/n0) ours/Z0 %.4f  (n2 ne/n1) %.4f | p_he/p_noHe %.4f  du %.3e erg/g (u_noHe %.3e)\n",
                    ρ, T, h0 / PC.yHe, z1 / PC.yHe, z2 / PC.yHe, xe, r1 / q1, r2 / q2, p1 / p0, u1 - u0, u0)
        end
    end
end
main(mode)
