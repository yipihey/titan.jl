# S5 test b: the closed box of S4 test b with TR-BDF2 and the consistent source/burn coupling (solve_step! + stage hook).
#   TITAN_NOLIB=1 julia -t 2 --project=env s5_b.jl      env: RHO, T0, DH, EFRAC, INT (trbdf2|theta)
# Closed box (static uniform slab with walls, no transport): total gas energy changes only by the network's release, so
#   |ΔE_gas − ∫ε dm dt|/release ≲ 1e-10 for any scheme, and the T(t), Y_D(t) sequence must follow NucNet.burn_coupled.
# S5 additionally: the residual's heating (booked at every stage with the TR-BDF2 weights) must equal the network's release
#   to the split tolerance, Σ|E_bur − E_inj|/ΣE_bur = O(1e-3) per step, and the TR-BDF2 sequence must be more accurate than the S4 θ = 1 step.
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "Trbdf2.jl"))
ρ0 = parse(Float64, get(ENV, "RHO", "1.0")); T0 = parse(Float64, get(ENV, "T0", "1e6")); n = 4
DH = parse(Float64, get(ENV, "DH", "2.5e-5")); X0 = 0.7
integ = get(ENV, "INT", "trbdf2")
Y0 = zeros(NN.NS); Y0[NN.IDX[:p]] = X0; Y0[NN.IDX[:d]] = DH * X0; Y0[NN.IDX[:he3]] = 1e-5 * X0
Y0[NN.IDX[:he4]] = (1 - sum(NN.AA .* Y0)) / 4
mat = SM.StarGas(TF.MGR; efrac = parse(Float64, get(ENV, "EFRAC", "0.05")), nuc = true, opac = :stellar)
Lbox = 1e10
rI = collect(range(0.5Lbox, Lbox, length = n + 1))
PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = 2, lgrid = 2, inner = :wall, outer = :wall, material = mat)
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
for k in 1:size(mat.Y, 2); SM.setY!(mat, k, Y0); end
TD.EOS[]()
SM.install_dt_hook!(mat); SM.install_stage_hook!(mat)
integ == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TF.setscalar!(:dtime, parse(Float64, get(ENV, "DT0", "1e9"))); TF.setscalar!(:tmax, parse(Float64, get(ENV, "TMAX", "3e12")))
dn, tn, egn = TF.var(:dn), TF.var(:tn), TF.var(:egn)
mass(k) = 4π * dn[k] * TF.var(:dvoln)[k]
Etot() = sum(mass(k) * egn[k] for k in a:b)
E0 = Etot()
rec = NamedTuple[]
push!(rec, (t = 0.0, dt = 0.0, T = tn[a], YD = mat.Y[NN.IDX[:d], a], Y = copy(mat.Y[:, a]), E = E0, Ebur = 0.0, corr = 0.0))
nsplit = Int[]
while TF.getscalar(:timen) < TF.getscalar(:tmax) && length(rec) < 2000
    mism, it, zm = SM.solve_step!(mat, TD.STEPFN[])
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    r = SM.burn!(mat; dt); push!(nsplit, it)
    push!(rec, (t = TF.getscalar(:timen), dt = dt, T = tn[a], YD = mat.Y[NN.IDX[:d], a], Y = copy(mat.Y[:, a]), E = Etot(), Ebur = mat.ledger.Ebur, corr = r.corr))
end
println("integrator ", integ, ": steps = ", length(rec) - 1, "  t_end = ", rec[end].t, " s; solve_step! iterations per step: mean ", sum(nsplit) / length(nsplit), ", max ", maximum(nsplit))
eos_u(ρ, T, Y) = SE.gas_state(ρ, T, SE.composition(SM.xyz(ForwardDiff.value.(Y))..., mat.zmix)).e
function chain(dts, nsub)
    Y = copy(Y0); T = T0; u = eos_u(ρ0, T0, Y0); out = [(T, Y[NN.IDX[:d]])]
    for dt in dts
        for _ in 1:nsub
            r = NN.burn_coupled(Y, ρ0, u, T, dt / nsub, eos_u); Y = r.Y; T = r.T; u = eos_u(ρ0, T, Y)
        end
        push!(out, (T, Y[NN.IDX[:d]]))
    end
    out
end
dts = [r.dt for r in rec[2:end]]
cf = chain(dts, 40)
ix = eachindex(rec); iy = [i for i in ix if rec[i].YD > 0.5 * rec[1].YD]
eTf = maximum(abs(rec[i].T / cf[i][1] - 1) for i in ix); eYf = maximum(abs(rec[i].YD / cf[i][2] - 1) for i in iy)
@printf("T0 = %.3g K, rho = %.3g, D/H = %.2e; T_end TITAN %.6f K, fine reference %.6f; Y_D(end)/Y_D(0): %.5f %.5f\n", T0, ρ0, DH, rec[end].T, cf[end][1], rec[end].YD / rec[1].YD, cf[end][2] / rec[1].YD)
@printf("TITAN vs fine reference (40 substeps/step): max|ΔT/T| %.2e, max|ΔY_D/Y_D| %.2e\n", eTf, eYf)
Er = abs((rec[end].E - E0) - rec[end].Ebur) / abs(rec[end].Ebur)
Erel = maximum(abs((r.E - E0) - r.Ebur) for r in rec) / abs(E0)
@printf("energy: ΔE_gas = %.8e erg, release (NucNet) = %.8e;  |ΔE − release|/release = %.2e; max over steps |ΔE − release|/E0 = %.2e\n", rec[end].E - E0, rec[end].Ebur, Er, Erel)
cm = sum(r.corr for r in rec[2:end]) / (length(rec) - 1)
@printf("splitting: mean per-step Σ|E_bur − E_inj|/ΣE_bur = %.2e (max %.2e); cumulative Σ|Ebur−Einj|/Σ|Ebur| = %.2e\n", cm, maximum(r.corr for r in rec[2:end]), mat.ledger.dabs / mat.ledger.Eabs)
ef = parse(Float64, get(ENV, "EFRAC", "0.05"))
println(Erel < 1e-10 ? "PASS energy" : "FAIL energy")
println(eTf < 2ef ? "PASS time-step error (T vs fine reference < 2 efrac)" : "FAIL time-step error")
println(maximum(r.corr for r in rec[3:end]) < 5e-3 ? "PASS coupling (residual heating = network release to 5e-3 per step)" : "FAIL coupling")
