using WasmTarget
const CNT = Ref(0)
function bump(k::Int32)::Int32
    CNT[] = 0                      # reset
    for i in 1:k; CNT[] += 1; end
    Int32(CNT[])
end
const CNT2 = Ref(5)                # nonzero initial value
function bump2(k::Int32)::Int32
    CNT2[] = 0
    for i in 1:k; CNT2[] += 1; end
    Int32(CNT2[])
end
b = compile_multi([(bump, (Int32,)), (bump2, (Int32,))]); write("out/bump.wasm", b)
println("native: ", bump(Int32(3)), " ", bump(Int32(3)), " | ", bump2(Int32(3)), " ", bump2(Int32(3)))
