# Phase 7a check: (1) Builtin material leaves the benchmark unchanged; (2) with the Z0 material the
# coloured-AD Jacobian equals the dense ForwardDiff Jacobian and central finite differences of the residual.
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl"))
include(joinpath(@__DIR__, "TitanMaterials.jl"))
using .TitanDriver; const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials
TF.load_state!(joinpath(TF.ROOT, "reference_runs", "init", "SUPC.tstate.gz"))
# move the state into the Z0 table's range (physics irrelevant here: this tests derivatives only)
a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
tn = TF.var(:tn); tn[a-1:b+2] .= 6000.0 .+ 3000.0 .* range(0, 1, length = b - a + 4)
ern = TF.var(:ern); ern[a-1:b+2] .= 7.56e-15 .* tn[a-1:b+2] .^ 4
TM.use!(TM.Z0Gas())
TD.EOS[](); TD.OPAC[]()
c = TX.context()
TX.assemble!(c; ncol = 6)
E = [copy(TF.var(s)) for s in (:em2, :em1, :e00, :ep1, :ep2)]
neqn, nb, k0 = c.neqn, c.nb, c.k0; N = neqn * nb
f(x) = (R = zeros(eltype(x), neqn, nb); TX.residual!(R, TX.getws(eltype(x), c.mgr), c, x); vec(R))
J = ForwardDiff.jacobian(f, zeros(N))
worst = 0.0
for kb in 1:nb, i in 1:neqn
    row = (kb - 1) * neqn + i; rowmax = maximum(abs, J[row, :]); rowmax == 0 && continue
    for cb in max(1, kb-2):min(nb, kb+2), j in 1:neqn
        global worst = max(worst, abs(J[row, (cb-1)*neqn+j] - E[cb-kb+3][i, j, k0+kb-1]) / rowmax)
    end
end
@printf("Z0 material: coloured AD vs dense AD Jacobian, max |Δ|/rowmax = %.2e  (N = %d unknowns)\n", worst, N)
# central finite differences on a sample of columns
h = 1e-6; fdw = 0.0
for col in 1:7:N
    e = zeros(N); e[col] = h
    d = (f(e) .- f(-e)) ./ (2h)
    for row in 1:N
        rowmax = maximum(abs, J[row, :]); rowmax == 0 && continue
        global fdw = max(fdw, abs(d[row] - J[row, col]) / rowmax)
    end
end
@printf("Z0 material: dense AD vs central FD (h = 1e-6, every 7th column), max |Δ|/rowmax = %.2e\n", fdw)
t = @elapsed for i in 1:20; TX.assemble!(c; ncol = 6); end
@printf("assemble with Z0 material: %.2f ms/call\n", 1e3 * t / 20)
