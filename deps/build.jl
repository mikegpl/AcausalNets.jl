#=
build:
- Julia version: 1.0
- Author: marcin
- Date: 2018-07-30
=#
using Pkg

if !("QuantumInformation" in keys(Pkg.installed()))
    Pkg.clone("https://github.com/ZKSI/QuantumInformation.jl")
end