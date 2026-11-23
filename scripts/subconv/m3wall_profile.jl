using Serialization, Printf
f = ARGS[1]; zlo = parse(Int, ARGS[2]); zhi = parse(Int, ARGS[3]); which = parse(Int, ARGS[4])
open(f) do io
  i = 0
  while !eof(io)
    s = deserialize(io); i += 1
    i == which || continue
    println("step ", s.step, " dt ", s.dt)
    a=7.5657e-15
    @printf("%4s %10s %10s %10s %10s %10s %10s %10s %10s\n","k","rho","T","Er/aT4","Fr","chi","u","pg","tau")
    for k in zlo:zhi
      @printf("%4d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",k,s.ρ[k],s.T[k],s.er[k]/(a*s.T[k]^4),s.fr[k],s.chi[k],s.u[k],s.pg[k], get(s.cd,:tau,fill(NaN,200))[k])
    end
    break
  end
end
