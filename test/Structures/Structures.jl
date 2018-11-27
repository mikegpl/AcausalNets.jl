#=
Structures:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-16
=#
using AcausalNets.Structures

@testset "Structures" begin
    for file in ["bayes_net.jl"]
        include(file)
    end
end
