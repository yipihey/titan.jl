using Serialization, Printf
d=ARGS[1]; S=Any[]
for f in sort(filter(x->occursin(r"^snaps_.*jls$",x),readdir(d)))
 open(joinpath(d,f)) do io; while !eof(io); try push!(S,deserialize(io)); catch; break; end; end; end
end
s=S[parse(Int,ARGS[2])]; c=s.cd; ar=7.565723e-15
println(propertynames(s))
for k in 2:7
 P(j)=s.pg[j]+s.fedd[j]*s.er[j]
 dlr=log(s.ρ[k]/s.ρ[k-1]); dlT=log(s.T[k]/s.T[k-1]); dlP=log(P(k)/P(k-1))
 Pi=0.5*(P(k)+P(k-1)); ρi=0.5*(s.ρ[k]+s.ρ[k-1]); lph=-ρi*c.g[k]*0.5*(s.r[k+1]-s.r[k-1])/Pi
 @printf("k=%d dlnrho=%.4g dlnT=%.4g dlnP=%.4g lph~%.4g ratio=%.4g | Er/aT4-1: %.3g %.3g | fedd %.5f | sup=%.3g gam1b=%.5f delb=%.3f\n",k,dlr,dlT,dlP,lph,dlP/lph,s.er[k-1]/(ar*s.T[k-1]^4)-1,s.er[k]/(ar*s.T[k]^4)-1,s.fedd[k],c.sup[k],c.gam1b[k],c.delb[k])
 A=dlr-dlP/c.gam1b[k]; @printf("   A=%.4g  A/(del|lph|)=%.4g ; (nab_raw-nad)=%.4g  dlT/dlP=%.4g\n",A,A/(c.delb[k]*abs(lph)), (c.nab_raw[k]-c.nad_fr[k]), dlT/dlP)
end
