using Pkg

Pkg.activate(; temp = true)
Pkg.add(url = "https://github.com/mroughan/IncCSV.jl")
Pkg.develop(path = dirname(@__DIR__))
Pkg.add("JET")

using JET
using S5

JET.test_package(S5; target_modules = (S5,))

