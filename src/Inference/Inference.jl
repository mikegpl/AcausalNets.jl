#=
Inference:
- Julia version: 
- Author: marcin
- Date: 2018-08-18
=#

module Inference
    include("clusterization.jl")
    include("join_tree.jl")
    include("evidence_propagation.jl")
    include("evidence.jl")
    include("inferrer.jl")
    include("naive_inferrer.jl")


    export
        Evidence,
        Inferrer,
        infer
end # module