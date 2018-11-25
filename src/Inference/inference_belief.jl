#=
belief_inferrer:
- Julia version: 
- Author: marcin
- Date: 2018-11-13
=#

 using LightGraphs

 import AcausalNets.Systems:
     ncategories,
     reduce_distribution,
     multiply_star,
     divide_star,
     identity_distribution,
     multiply_kron,
     permute_distribution,
     sub_system

 import AcausalNets.Structures:
     DiscreteBayesNet

 import AcausalNets.Inference:
     JoinTree,
     shallowcopy,
     enforce_clique,
     triangulate,
     apply_observations,
     global_propagation

import AcausalNets.Common:
    Variable,
    ncategories


function infer_belief(
        dbn::DiscreteBayesNet{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[]
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
    length(vars_to_infer) > 0 || error("At least one variable to infer must be specified!")
    observations_jt = unpropagated_join_tree(
                        dbn,
                        vars_to_infer,
                        observations
                    )
    inferred_cluster_ind = first([
            i
            for (i, sys) in observations_jt.vertex_to_cluster
            if all([
                v in variables(sys)
                for v in vars_to_infer
            ])
        ])
    messages_no = 2 * diameter(observations_jt.graph)
    inferred_cluster = belief(observations_jt, inferred_cluster_ind, messages_no)
    inference_result = sub_system(inferred_cluster, vars_to_infer)
    intermediate_elements = (
        observations_jt,
        messages_no,
        inferred_cluster_ind,
        inferred_cluster
    )
    return inference_result, intermediate_elements
end

function message(jt::JoinTree{S}, from::Int, to::Int, t::Int) where {D, S <: DiscreteSystem{D}}
    from_neighbors = neighbors(jt.graph, from)
    to in from_neighbors || error("$from and $to are not neighbors!")
    if t == 0
        return identity_system(jt.vertex_to_cluster[to])
    end
    from_system = jt.vertex_to_cluster[from] # (u_u)a
    to_system = jt.vertex_to_cluster[to]
    # equation 103 craziness
    previous_messages =  vcat(
        [identity_system(from_system)],
        S[message(jt, n, from, t-1) for n in from_neighbors if n !=to]
    )
    msg_distribution = multiply_star(
                    distribution(from_system),
                    multiply_star(
                        D(prod([distribution(m) for m in previous_messages])),
                        distribution(mutual_system(from_system, to_system))
                    )
    )
    msg_from = S(variables(from_system), msg_distribution)
    msg_to = sub_system(msg_from, variables(to_system))
    normalized_dist = distribution(msg_to) / tr(distribution(msg_to))
    S(variables(msg_to), D(normalized_dist))
end

function belief(jt::JoinTree{S}, cluster_ind::Int, t::Int)::S where  {D, S <: DiscreteSystem{D}}
    sys = jt.vertex_to_cluster[cluster_ind]
    messages = vcat(
        [identity_system(sys)],
        [message(jt, n, cluster_ind, t) for n in neighbors(jt.graph, cluster_ind)]
    )
    unnormalized_dist = multiply_star(
        distribution(sys),
        D(prod([distribution(m) for m in messages]))
    )
    normalized_dist = unnormalized_dist / tr(unnormalized_dist)
    S(variables(sys), D(normalized_dist))
end

function mutual_system(from::S, to::S) where {D, S <: DiscreteSystem{D}}
    # TODO implement, move to Systems module
    identity_system(from)
end