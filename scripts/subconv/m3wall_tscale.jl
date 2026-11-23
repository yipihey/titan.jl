using Serialization, Printf
f = ARGS[1]; which = parse(Int, ARGS[2])
open(f) do io
  i = 0
  while !eof(io)
    s = deserialize(io); i += 1
    i == which || continue
    c = s.cd; g(sym) = c[sym]
    @printf("step %d dt %.3g\n", s.step, s.dt)
    @printf("%3s %9s %9s %9s %9s | %9s %9s %9s %9s %9s %9s\n","k","w","ell","a","sup","l/w","l/a","t_grow","dr/|u|","dr2/(at l w)","l/(2Ua)")
    for k in 31:41
      w = g(:tdc_w)[k]; l = g(:lmix)[k]; a = g(:tdc_a)[k]; sup = g(:sup)[k]; U = g(:tdc_U)[k]
      dr = s.r[k+1]-s.r[k]; u = abs(s.u[k])
      tg = sup>0 ? l/(a*sqrt(sup)) : Inf
      @printf("%3d %9.2e %9.2e %9.2e %9.2e | %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",k,w,l,a,sup, w>0 ? l/w : Inf, l/a, tg, dr/u, w>0 ? dr^2/(0.01*l*w) : Inf, U>0 ? l/(2U*a) : Inf)
    end
    break
  end
end
