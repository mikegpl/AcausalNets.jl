#=
Inference:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-16
=#

module Inference
    for file in ["clusterization.jl", "evidence.jl", "inference_api.jl", "infer_naive.jl"] # todo belief and join_tree
        include(file)
    end
end