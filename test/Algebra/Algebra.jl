#=
Algebra:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-15
=#
using AcausalNets.Algebra

@testset "Algebra" begin
    for file in ["common", "distribution_operations"]
        include("$file.jl")
    end
end