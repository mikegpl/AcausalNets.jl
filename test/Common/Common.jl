#=
Common:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-16
=#
using AcausalNets.Common

@testset "Common" begin
    for file in ["sets.jl", "variables.jl"]
        include(file)
    end
end