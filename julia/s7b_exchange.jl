# S7b: exchange records between the envelope and the protostar process (key = value text files written atomically). No dependencies.
module S7bX
"Write a record (key = value lines) atomically (temporary file + rename)."
function write_rec(path::AbstractString, d::AbstractDict)
    tmp = path * ".tmp$(getpid())"
    open(tmp, "w") do io
        for k in sort(collect(keys(d))); println(io, k, " = ", d[k]); end
    end
    mv(tmp, path; force = true); nothing
end
function read_rec(path::AbstractString)
    d = Dict{String,Float64}()
    for l in eachline(path)
        occursin('=', l) || continue
        k, v = split(l, "="; limit = 2); d[strip(k)] = parse(Float64, strip(v))
    end
    d
end
"Block until `path` exists (polling), then read it; error after `timeout` seconds."
function wait_rec(path::AbstractString; poll = 0.25, timeout = 7200.0)
    t0 = time()
    while !isfile(path)
        time() - t0 > timeout && error("S7b exchange: timeout waiting for $path")
        sleep(poll)
    end
    for _ in 1:40                         # (the rename is atomic; a read error would mean a foreign writer)
        try return read_rec(path) catch; sleep(0.05); end
    end
    read_rec(path)
end

end # module
