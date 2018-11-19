#=
inference_api:
- Julia version: 
- Author: marcin
- Date: 2018-11-14
=#
import AcausalNets.Inference:
    infer_join_tree, # TODO deprecate
    infer_naive, # The only flawless (albeit slow) version
    infer_belief # TODO make it work

import AcausalNets.Structures:
    DiscreteBayesNet

function infer(
        dbn::DiscreteBayesNet{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[],
        inference_strategy::Function = infer_naive # not parametrizable yet, should be (DiscreteBayesNet{S}, Vector{Variable}, Vector{E})::S
        )::S where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
        result, _ = inference_strategy(dbn, vars_to_infer, observations)
        return result
end