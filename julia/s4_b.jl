# S4 test b: closed box (static uniform slab with walls, no transport), TITAN step sequence vs NucNet.burn_coupled.
include(joinpath(@__DIR__, "s4_common.jl"))
ρ0 = parse(Float64, get(ENV, "RHO", "1.0")); T0 = parse(Float64, get(ENV, "T0", "1e6")); n = 4
DH = parse(Float64, get(ENV, "DH", "2.5e-5")); X0 = 0.7
Y0 = zeros(NN.NS); Y0[NN.IDX[:p]] = X0; Y0[NN.IDX[:d]] = DH * X0; Y0[NN.IDX[:he3]] = 1e-5 * X0
Y0[NN.IDX[:he4]] = (1 - sum(NN.AA .* Y0)) / 4
mat = SM.StarGas(TF.MGR; efrac = parse(Float64, get(ENV, "EFRAC", "0.05")), nuc = get(ENV, "NUC", "1") == "1", opac = :stellar)
Lbox = 1e10
rI = collect(range(get(ENV, "GEOM", "2") == "0" ? 0.0 : 0.5Lbox, Lbox, length = n + 1))
PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = parse(Int, get(ENV, "GEOM", "2")), lgrid = 2, inner = :wall, outer = :wall,
          material = mat)
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
for k in 1:size(mat.Y, 2); SM.setY!(mat, k, Y0); end
TD.EOS[]()
SM.install_dt_hook!(mat)
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TF.setscalar!(:dtime, parse(Float64, get(ENV, "DT0", "1e9"))); TF.setscalar!(:tmax, parse(Float64, get(ENV, "TMAX", "3e12")))
dn, tn, egn = TF.var(:dn), TF.var(:tn), TF.var(:egn)
mass(k) = 4π * dn[k] * TF.var(:dvoln)[k]
Etot() = sum(mass(k) * egn[k] for k in a:b)
E0 = Etot()
rec = NamedTuple[]
push!(rec, (t = 0.0, dt = 0.0, T = tn[a], YD = mat.Y[NN.IDX[:d], a], Y = copy(mat.Y[:, a]), E = E0, Ebur = 0.0))
TD.run_julia!(maxsteps = 100000, output = function (k)
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    SM.burn!(mat; dt)
    push!(rec, (t = TF.getscalar(:timen), dt = dt, T = tn[a], YD = mat.Y[NN.IDX[:d], a], Y = copy(mat.Y[:, a]), E = Etot(), Ebur = mat.ledger.Ebur))
end)
println("steps = ", length(rec) - 1, "  t_end = ", rec[end].t, " s")
# reference: NucNet.burn_coupled chain with the same dt sequence, and a fine reference (40 substeps per step)
eos_u(ρ, T, Y) = SE.gas_state(ρ, T, SE.composition(SM.xyz(ForwardDiff.value.(Y))..., mat.zmix)).e   # composition dependence of u not differentiated (quasi-Newton)
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
c1 = chain(dts, 1); cf = chain(dts, 40)
ix = eachindex(rec); iy = [i for i in ix if rec[i].YD > 0.5 * rec[1].YD]    # Y_D compared while D is above half its initial value (afterwards the relative error is a timing shift)
eT1 = maximum(abs(rec[i].T / c1[i][1] - 1) for i in ix); eY1 = maximum(abs(rec[i].YD / c1[i][2] - 1) for i in iy)
eTf = maximum(abs(rec[i].T / cf[i][1] - 1) for i in ix); eYf = maximum(abs(rec[i].YD / cf[i][2] - 1) for i in iy)
eTc = maximum(abs(c1[i][1] / cf[i][1] - 1) for i in ix); eYc = maximum(abs(c1[i][2] / cf[i][2] - 1) for i in iy)
@printf("T0 = %.3g K, rho = %.3g, D/H = %.2e; T_end TITAN %.6f K, burn_coupled chain %.6f, fine ref %.6f; Y_D(end)/Y_D(0): %.5f %.5f %.5f\n",
        T0, ρ0, DH, rec[end].T, c1[end][1], cf[end][1], rec[end].YD / rec[1].YD, c1[end][2] / rec[1].YD, cf[end][2] / rec[1].YD)
@printf("TITAN vs same-dt burn_coupled chain: max|ΔT/T| %.2e, max|ΔY_D/Y_D| %.2e\n", eT1, eY1)
@printf("TITAN vs fine reference (40 substeps): max|ΔT/T| %.2e, max|ΔY_D/Y_D| %.2e;  (chain vs fine: %.2e, %.2e)\n", eTf, eYf, eTc, eYc)
Er = abs((rec[end].E - E0) - rec[end].Ebur) / abs(rec[end].Ebur)
@printf("energy: ΔE_gas = %.8e erg, ∫ε dm dt (NucNet) = %.8e;  |ΔE − release|/release = %.2e;  E_gas/E0 − 1 = %.4f\n",
        rec[end].E - E0, rec[end].Ebur, Er, rec[end].E / E0 - 1)
Erel = maximum(abs((r.E - E0) - r.Ebur) for r in rec) / abs(E0)
@printf("max over steps |ΔE − release|/E0 = %.2e;  ledger: Einj %.6e Ebur %.6e splitting correction Σ|Ebur−Einj|/Σ|Ebur| %.2e, resid(max) %.2e, %d steps\n",
        Erel, mat.ledger.Einj, mat.ledger.Ebur, mat.ledger.dabs / mat.ledger.Eabs, mat.ledger.resid, mat.ledger.nsteps)
mtot = 4π * sum(dn[k] * TF.var(:dvoln)[k] for k in a:b)
dB = sum((rec[end].Y[i] - Y0[i]) * SM.BIND[i] for i in 1:NN.NS) * NN.NA * NN.MeV
@printf("mass-defect check: N_A Σ ΔY_i B_i = %.6e erg/g vs net release %.6e erg/g (+ neutrinos): ratio %.5f\n", dB, rec[end].Ebur / mtot, dB / (rec[end].Ebur / mtot))
ef = parse(Float64, get(ENV, "EFRAC", "0.05"))
println(Erel < 1e-10 ? "PASS energy" : "FAIL energy")
println(max(eT1, eY1) < 2ef ? "PASS sequence (TITAN vs burn_coupled, same dt: < 2 efrac)" : "FAIL sequence")
println(eTf < 2ef ? "PASS time-step error (T vs fine reference < 2 efrac)" : "FAIL time-step error")
if get(ENV, "TRACE", "0") == "1"
    for i in ix; @printf("%3d t %.4e dt %.3e T %.6e (chain %.6e fine %.6e)  YD/YD0 %.4e (%.4e %.4e)\n", i - 1, rec[i].t, rec[i].dt, rec[i].T, c1[i][1], cf[i][1], rec[i].YD / rec[1].YD, c1[i][2] / rec[1].YD, cf[i][2] / rec[1].YD); end
end
