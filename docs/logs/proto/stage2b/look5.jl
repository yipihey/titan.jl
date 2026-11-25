using Serialization, Printf
d=ARGS[1]; S=Any[]
for f in sort(filter(x->occursin(r"^snaps_.*jls$",x),readdir(d)))
 open(joinpath(d,f)) do io; while !eof(io); try push!(S,deserialize(io)); catch; break; end; end; end
end
for si in (1,length(S)); s=S[si]; c=s.cd; println("snap ",si," step ",s.step," dt ",s.dt)
 for k in 2:3:40
  @printf("%3d r=%.3g T=%.3g rho=%.3g w=%.3g sup=%.3g a=%.3g aζ~%.3g U=%.3g lmix=%.3g HP=%.3g del=%.3g G1=%.4f nabraw=%.3g nabrad=%.3g PtP=%.2g Fc/Ft=%.2g\n",k,s.r[k],s.T[k],s.ρ[k],c.tdc_w[k],c.sup[k],c.tdc_a[k],c.tdc_a[k]*(sup=c.sup[k]; sup>0 ? sup/(sqrt(c.tdc_U[k]^2+sup)+c.tdc_U[k]) : 0.0),c.tdc_U[k],c.lmix[k],c.HP[k],c.delb[k],c.gam1b[k],c.nab_raw[k],c.nab_rad[k],c.tdc_PtP[min(k,end)],c.Fconv_over_Ftot[k])
 end; end
