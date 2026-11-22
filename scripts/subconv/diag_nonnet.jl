# unstable subsonic interfaces NOT in network zones (branch != 1) at a step: julia diag_nonnet.jl <snaps> <step>
using Serialization, Printf
v = Any[]; open(ARGS[1]) do io; while !eof(io); push!(v, deserialize(io)); end; end
s = v[findlast(x -> x.step == parse(Int, ARGS[2]), v)]; cd = s.cd; js = findfirst(<(-2e5), s.u)
@printf("step %s shock interface %d\n", ARGS[2], js)
for j in 2:js-1
    (cd.mach_i[j] < 0.5 && (cd.N2_L[j] < 0 || cd.N2_eq[j] < 0) && cd.branch[j] != 1) || continue
    @printf("j=%d r/R*=%.4f T=%.0f br=%d N2_L=%.2e N2_eq=%.2e nab=%.3f nad_used=%.3f w_th=%g Fconv/F=%.2e Fpot_th/F=%.2e\n", j, s.r[j]/s.r[js], 0.5*(s.T[j-1]+s.T[j]), cd.branch[j], cd.N2_L[j], cd.N2_eq[j], cd.nab[j], cd.nad_used[j], cd.w_th[j], cd.Fconv_over_Ftot[j], cd.Fpot_over_Ftot[j])
end
