#=
Inference:
- Julia version: 0.7
- Author: marcin
- Date: 2018-08-18
=#

module Inference
    include("clusterization.jl")
    include("join_tree.jl")
#     include("evidence_propagation.jl")
    include("evidence.jl")
    include("inference_join_tree.jl")
    include("inference_naive.jl")
    include("inference_belief.jl")
    include("inference_api.jl")

    export
        Evidence,
        infer,
        infer_join_tree,
        infer_naive,
        infer_belief
end # module