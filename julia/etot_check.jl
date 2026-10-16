# Item 2: tier-(a)-style check of the total-energy residual: coloured AD Jacobian (6 colours)
# vs dense ForwardDiff, at a mid-Newton state; out-of-band dependence of the new row.
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "TitanResidFast.jl"))
using .TitanDriver; const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast
name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TD.bitwise!(true); TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p4_runs", "ec_$name"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update)
c = TX.context(; etotal = true)
TX.assemble!(c; ncol = 6)
E = [copy(TF.var(s)) for s in (:em2, :em1, :e00, :ep1, :ep2)]
neqn, nb, k0 = c.neqn, c.nb, c.k0; N = neqn * nb
f(x) = (R = zeros(eltype(x), neqn, nb); TX.residual!(R, TX.getws(eltype(x), c.mgr), c, x); vec(R))
J = ForwardDiff.jacobian(f, zeros(N))
it = c.ix[5]; inb = 0.0; outb = 0.0; worst = 0.0
for kb in 1:nb, i in 1:neqn
    row = (kb - 1) * neqn + i
    rowmax = maximum(abs, J[row, :])
    rowmax == 0 && continue
    for cb in 1:nb, j in 1:neqn
        v = J[row, (cb-1)*neqn+j]
        if abs(cb - kb) <= 2
            d = abs(v - E[cb-kb+3][i, j, k0+kb-1]) / rowmax
            global worst = max(worst, d)
        elseif i == it
            global outb = max(outb, abs(v) / rowmax)
        end
    end
end
@printf("%s step %d: etot rows coloured vs dense AD max |Δ|/rowmax = %.2e; energy-row dependence outside ±2 blocks = %.2e\n",
        name, nstep, worst, outb)
TF.finalize()
