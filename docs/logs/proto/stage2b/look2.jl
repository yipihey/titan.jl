using Serialization, Printf
d=ARGS[1]; S=Any[]
for f in sort(filter(x->occursin(r"^snaps_.*jls$",x),readdir(d)))
 open(joinpath(d,f)) do io; while !eof(io); push!(S,deserialize(io)); end; end
end
for si in (1,length(S)); s=S[si]; c=s.cd; println("snap ",si," step ",s.step," dt ",s.dt)
 for k in 2:8:42
  @printf("%3d r=%.3g T=%.3g rho=%.3g frad=%.3g Ftot=%.3g Fconv=%.3g nabrad=%.4g nabraw=%.4g sup=%.3g nadfr=%.4f dfr=%.3f nadeq=%.4f cpfr=%.3g g=%.3g\n",k,s.r[k],s.T[k],s.ρ[k],c.frad[k],c.Ftot[k],c.fconv[k],c.nab_rad[k],c.nab_raw[k],c.sup[k],c.nad_fr[k],c.del_fr[k],c.nad_eq[k],c.cp_fr[k],c.g[k])
 end; end
