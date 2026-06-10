using Pkg

Pkg.activate(; temp = true)
Pkg.add(url = "https://github.com/mroughan/IncCSV.jl")
Pkg.develop(path = dirname(@__DIR__))
Pkg.add("Aqua")

using Aqua
using S5

Aqua.test_all(S5; stale_deps = (ignore = [:StableRNGs],))

