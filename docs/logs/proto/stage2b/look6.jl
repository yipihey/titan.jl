using Serialization, Printf
d=ARGS[1]; S=Any[]
for f in sort(filter(x->occursin(r"^snaps_.*jls$",x),readdir(d)))
 open(joinpath(d,f)) do io; while !eof(io); try push!(S,deserialize(io)); catch; break; end; end; end
end
for k in (9,14,17,19)
 println("k=$k"); for s in S; c=s.cd
  @printf(" st %d dt=%.3g w=%.4g sup=%.4g a=%.3g rho=%.4g T=%.5g u=%.3g nabraw=%.3g\n",s.step,s.dt,c.tdc_w[k],c.sup[k],c.tdc_a[k],s.ρ[k],s.T[k],s.u[k],c.nab_raw[k])
 end; end
