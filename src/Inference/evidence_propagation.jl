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
