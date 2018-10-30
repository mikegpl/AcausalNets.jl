#=
join_tree_messaging:
- Julia version: 
- Author: marcin
- Date: 2018-08-19
=#

using LightGraphs
using QI

import AcausalNets.Systems:
    ncategories,
    reduce_distribution,
    multiply_star,
    divide_star,
    identity_distribution,
    multiply_kron,
    permute_distribution,
    sub_system

import AcausalNets.Algebra:
    star, unstar, event

import AcausalNets.Structures:
    DiscreteBayesNet

import AcausalNets.Inference:
    JoinTree,
    shallowcopy

function single_message_pass(from_ind::Int, to_ind::Int, jt::JoinTree{S}) where S
    if (from_ind, to_ind) in edges(jt.graph)
        println("message2 from $from_ind to $to_ind")
        jt = shallowcopy(jt)
        cluster_from = jt.vertex_to_cluster[from_ind]
        cluster_to = jt.vertex_to_cluster[to_ind]
        edge_set = Set([from_ind, to_ind])
        sepset = jt.edge_to_sepset[edge_set]

        old_sepset = sub_system(sepset, variables(cluster_to))
        new_sepset = sub_system(cluster_from, variables(cluster_to))
        jt.edge_to_sepset[edge_set] = sub_system(new_sepset, variables(sepset))

        to_distribution = distribution(cluster_to)

        new_to_distribution = multiply_star(
            divide_star(to_distribution, distribution(old_sepset)),
            distribution(new_sepset),
        )


        new_cluster_to = S(variables(cluster_to), new_to_distribution)
#         cluster_to = S(variables(cluster_to), event(to_distribution, message_ordered))
#         cluster_to = S(variables(cluster_to), multiply_star(to_distribution, message_ordered))
        jt.vertex_to_cluster[to_ind] = new_cluster_to
    end
    return jt
end


function collect_evidence(cluster_ind::Int, cluster_marks::Vector{Bool}, jt::JoinTree)
    jt = shallowcopy(jt)
    cluster_marks[cluster_ind] = false
    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt, cluster_marks = collect_evidence(neighbor, cluster_marks, jt)
            jt = single_message_pass(neighbor, cluster_ind, jt)
        end

    end
    jt, cluster_marks

end

function distribute_evidence(cluster_ind::Int, cluster_marks::Vector{Bool}, jt::JoinTree)
    jt = shallowcopy(jt)
    cluster_marks[cluster_ind] = false
    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt = single_message_pass(cluster_ind, neighbor, jt)
            # pass a message from cluster_ind to neighbor
        end
    end
    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt, cluster_mars = distribute_evidence(neighbor, cluster_marks, jt)
        end
    end
    jt, cluster_marks
end

function global_propagation(jt::JoinTree, start_ind=1)
    jt = shallowcopy(jt)
    cluster_marks = [true for k in keys(jt.vertex_to_cluster)]
    jt, cluster_marks = collect_evidence(start_ind, cluster_marks, jt)
    cluster_marks = [true for k in keys(jt.vertex_to_cluster)]
    jt, cluster_marks = distribute_evidence(start_ind, cluster_marks, jt)
    return jt
end

function message(jt::JoinTree{S}, from::Int, to::Int, t::Int) where S
#     println("message from $from to $to t=$t")
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
                        prod([distribution(m) for m in previous_messages]),
                        distribution(identity_system(from_system))
#                         distribution(
#                                 sub_system(
#                                     jt.edge_to_sepset[Set([from, to])],
#                                     variables(from_system)
#                                     )
#                             )
                    )
    )
    msg_from = S(variables(from_system), msg_distribution)
    msg_to = sub_system(msg_from, variables(to_system))
    normalized_dist = distribution(msg_to) / tr(distribution(msg_to))
    S(variables(msg_to), normalized_dist)
end


function belief(jt::JoinTree{S}, cluster_ind::Int, t::Int)::S where S
    sys = jt.vertex_to_cluster[cluster_ind]
    messages = vcat(
        [identity_system(sys)],
        [message(jt, n, cluster_ind, t) for n in neighbors(jt.graph, cluster_ind)]
    )
    unnormalized_dist = multiply_star(
        distribution(sys),
        prod([distribution(m) for m in messages])
    )
    normalized_dist = unnormalized_dist / tr(unnormalized_dist)
    S(variables(sys), normalized_dist)
end