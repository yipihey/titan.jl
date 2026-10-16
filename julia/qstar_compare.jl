# Compare a TITAN quasi-star snapshot with the steady QuasiStarInfall envelope it was built from.
#   julia qstar_compare.jl SNAPS.jls PROFILE.csv [index]
# Prints, at the steady-model points inside the TITAN domain: T, ρ, and the convective fraction
# L_conv/L (steady model: 1 − ∇/∇_rad; TITAN: 4πr²F_conv / L from the residual's MLT flux).
using Serialization, Printf
module QuasiStarSetup   # the two helpers of QuasiStarSetup.jl, without its TITAN dependencies
src = read(joinpath(@__DIR__, "QuasiStarSetup.jl"), String)
for name in ("function read_profile", "function linterp")
    i = findfirst(name, src)[1]; j = findnext("\nend\n", src, i)[end]
    eval(Meta.parse(src[i:j]))
end
end
snaps = deserialize(ARGS[1]); s = snaps[length(ARGS) >= 3 ? parse(Int, ARGS[3]) : lastindex(snaps)]
P = QuasiStarSetup.read_profile(ARGS[2]); C = P.cols; L = P.hdr[:L]; R = P.hdr[:R]
li(x, y, q) = QuasiStarSetup.linterp(x, y, q)
rc = [((s.r[i]^3 + s.r[i+1]^3) / 2)^(1 / 3) for i in eachindex(s.d)]
@printf("snapshot t = %.4e s\n%8s %10s %10s %7s %10s %10s %7s %9s %9s\n", s.t, "r/R", "T_steady", "T_titan", "dT/T",
        "ρ_steady", "ρ_titan", "dρ/ρ", "Lc/L st", "Lc/L ti")
for i in length(C[:r]):-1:1
    r = C[:r][i]; (r < s.r[1] || r > s.r[end]) && continue
    Tt = exp(li(rc, log.(s.T), r)); ρt = exp(li(rc, log.(s.d), r))
    lcs = max(0.0, 1 - C[:nabla][i] / C[:nabla_rad][i]); lct = li(s.r, s.Lconv, r) / L
    @printf("%8.5f %10.1f %10.1f %7.4f %10.3e %10.3e %7.4f %9.4f %9.4f\n", r / R, C[:T][i], Tt, Tt / C[:T][i] - 1,
            C[:rho][i], ρt, ρt / C[:rho][i] - 1, lcs, lct)
end
