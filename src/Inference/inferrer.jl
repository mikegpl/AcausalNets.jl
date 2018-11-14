#=
huang_inference:
- Julia version: 1.0
- Author: marcin
- Date: 2018-09-11
=#

import AcausalNets.Common:
    Variable,
    ncategories

import AcausalNets.Systems:
    DiscreteSystem,
    variables,
    distribution,
    permute_system

import AcausalNets.Structures:
    DiscreteBayesNet

import AcausalNets.Inference:
    JoinTree,
    normalize,
    ParentCliquesDict,
    parent_cliques_dict,
    moral_graph,
    enforce_clique,
    triangulate,
    apply_observations,
    global_propagation


struct Inferrer{S <: DiscreteSystem}
    bayes_net                   ::DiscreteBayesNet{S}

    function Inferrer{S}(dbn::DiscreteBayesNet{S}) where S
        dbn = deepcopy(dbn)
        new{S}(dbn)
    end
end

Inferrer(dbn::DiscreteBayesNet{S}) where S = Inferrer{S}(dbn)
variables(inferrer::Inferrer) = variables(inferrer.bayes_net)

function infer(
        inferrer::Inferrer{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[],
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
    result, _ = infer_debug(inferrer, vars_to_infer, observations)
    return result
end

function infer_debug(
        inferrer::Inferrer{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[]
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }

    length(vars_to_infer) > 0 || error("At least one variable to infer must be specified!")
    dbn = inferrer.bayes_net
    mg = moral_graph(dbn)
    enforced_mg = enforce_clique(dbn, mg, vars_to_infer)
    tri_mg, cliques = triangulate(enforced_mg, dbn)
    parent_cliques = parent_cliques_dict(cliques, dbn)
    initialized_jt = JoinTree(cliques, dbn)

    observations_jt = apply_observations(
                        initialized_jt,
                        parent_cliques,
                        observations
                    )
    propagated_jt = global_propagation(observations_jt)
    jt = normalize(propagated_jt)
    inferred_cluster = first([
            sys
            for (i, sys) in jt.vertex_to_cluster
            if all([
                    v in variables(sys)
                    for v in vars_to_infer
                    ])
        ])

    inference_result = sub_system(inferred_cluster, vars_to_infer)
    intermediate_elements = (
        dbn,
        mg,
        enforced_mg,
        tri_mg,
        cliques,
        parent_cliques,
        initialized_jt,
        observations_jt,
        propagated_jt,
        jt,
    )

    inference_result, intermediate_elements
end