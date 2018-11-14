#=
build:
- Julia version: 1.0
- Author: marcin
- Date: 2018-07-30
=#
using Pkg

if !("QI" in keys(Pkg.installed()))
    Pkg.clone("https://github.com/ZKSI/QI.jl")
end