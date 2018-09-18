#=
huang_inference:
- Julia version: 
- Author: marcin
- Date: 2018-09-11
=#

import AcausalNets.Common:
    Variable

import AcausalNets.Systems:
    DiscreteSystem,
    variables

import AcausalNets.Structures:
    DiscreteBayesNet

import AcausalNets.Inference:
    JoinTree,
    normalize,
    ParentCliquesDict,
    parent_cliques_dict,
    moral_graph,
    triangulate,
    apply_observations,
    global_propagation


struct Inferrer{S <: DiscreteSystem}
    bayes_net                   ::DiscreteBayesNet{S}
    initialized_join_tree       ::JoinTree{S}
    parent_cliques              ::ParentCliquesDict{S}

    function Inferrer{S}(dbn::DiscreteBayesNet{S}) where S
        dbn = deepcopy(dbn)
        mg = moral_graph(dbn)
        tri_mg, cliques = triangulate(mg, dbn)
        parent_cliques = parent_cliques_dict(cliques, dbn)
        join_tree = JoinTree(cliques, dbn)
        new{S}(dbn, join_tree, parent_cliques)
    end
end

Inferrer(dbn::DiscreteBayesNet{S}) where S = Inferrer{S}(dbn)

 # TODO return a system of only the specified variables
function infer(
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
    parent_cliques = inferrer.parent_cliques
    jt = normalize(
            global_propagation(
                apply_observations(
                    inferrer.initialized_join_tree,
                    parent_cliques,
                    observations
                )
            )
        )
    first([
            sys
            for (i, sys) in jt.vertex_to_cluster
            if all([
                    v in variables(sys)
                    for v in vars_to_infer
                    ])
        ])
end