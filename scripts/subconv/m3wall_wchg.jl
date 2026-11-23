using Serialization, Printf
global prev = nothing
open(ARGS[1]) do io
  while !eof(io)
    s = deserialize(io)
    w = s.cd.tdc_w
    if prev !== nothing
      @printf("%4d dt %.2e  w35 %.3e w36 %.3e w37 %.3e | dw/w over the %d steps: %.3f %.3f %.3f | Hp37 %.3e\n", s.step, s.dt, w[35], w[36], w[37], s.step-prev.step, (w[35]-prev.w[35])/max(w[35],1), (w[36]-prev.w[36])/max(w[36],1), (w[37]-prev.w[37])/max(w[37],1), s.Y[8,37])
    end
    global prev = (w = copy(w), step = s.step)
  end
end
