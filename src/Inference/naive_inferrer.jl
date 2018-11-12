#=
naive_inference:
- Julia version: 0.7
- Author: marcin
- Date: 2018-10-13
=#
import AcausalNets.Systems:
    DiscreteSystem,
    merge_systems,
    reduce_distribution,
    permute_system,
    sub_system

import AcausalNets.Inference:
    apply_evidence

function infer_naive(
        inferrer::Inferrer{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[],
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
    result, _ = infer_naive_debug(inferrer, vars_to_infer, observations)
    return result
end

function infer_naive_debug(
        inferrer::Inferrer{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[]
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
    """
    Naive version of inference, which merges all distributions into one and traces out the result.
    This is simple, but effectively creates a system of all possible variable combinations.
    # TODO this can be simplified, since there's no join tree involved
    """
    length(vars_to_infer) > 0 || error("At least one variable to infer must be specified!")
    dbn = inferrer.bayes_net
    cluster = merge_systems(systems(dbn))
    observations = merge_systems(observations)
    evidence_cluster = apply_evidence(cluster, observations)
    inferred_cluster = S(
        variables(evidence_cluster),
        distribution(evidence_cluster) / tr(distribution(evidence_cluster))
        )
    inference_result = sub_system(inferred_cluster, vars_to_infer)
    intermediate_elements = (
        dbn,
        cluster,
        observations,
        evidence_cluster,
        inferred_cluster
    )

#     inference_result = permute_system(inferred_system, new_variable_indexing)
    inference_result, intermediate_elements
end