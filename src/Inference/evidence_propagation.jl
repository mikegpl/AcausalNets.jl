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
    multiply_star
    divide_star,
    identity_distribution,
    multiply_kron,
    permute_distribution

import AcausalNets.Algebra:
    star, unstar

import AcausalNets.Structures:
    DiscreteBayesNet

import AcausalNets.Inference:
    JoinTree,
    shallowcopy

function single_message_pass(from_ind::Int, to_ind::Int, jt::JoinTree{S}) where S
    if (from_ind, to_ind) in edges(jt.graph)
        jt = shallowcopy(jt)
        cluster_from = jt.vertex_to_cluster[from_ind]
        cluster_to = jt.vertex_to_cluster[to_ind]
        edge_set = Set([from_ind, to_ind])
        sepset = jt.edge_to_sepset[edge_set]
        to_trace_out_vars = setdiff(variables(cluster_from), variables(sepset))
        to_trace_out_ind = Int64[
            findfirst([v==var for var in variables(cluster_from)]) for v in to_trace_out_vars
            ]
        from_variables_sizes = [ncategories(v) for v in variables(cluster_from)]
        from_distribution = distribution(cluster_from)
        old_sepset_distribution = distribution(sepset)
        new_sepset_distribution = reduce_distribution(
            distribution(cluster_from), from_variables_sizes, to_trace_out_ind
        )
        sepset = S(variables(sepset), new_sepset_distribution)
        jt.edge_to_sepset[edge_set] = sepset

        to_distribution = distribution(cluster_to)

        message = divide_star(new_sepset_distribution, old_sepset_distribution)

        message_vars = variables(sepset)
        non_message_vars = setdiff(variables(cluster_to), message_vars)
        non_message_size = prod([ncategories(v) for v in non_message_vars])
        message = multiply_kron(
            message, identity_distribution(typeof(message), non_message_size)
        )
        message_all_vars = vcat(message_vars, non_message_vars)

        message_sorted_indices = Int64[
            findfirst([v==var for var in message_all_vars])
            for v in variables(cluster_to)
            ]

        message_dims = [ncategories(v) for v in message_all_vars]
        message_ordered = permute_distribution(message, message_dims, message_sorted_indices)
        cluster_to = S(variables(cluster_to), multiply_star(to_distribution, message_ordered))
        jt.vertex_to_cluster[to_ind] = cluster_to
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

function global_propagation(jt::JoinTree)
    jt = shallowcopy(jt)
    cluster_marks = [true for k in keys(jt.vertex_to_cluster)]
    arbitrary_cluster_ind = 1
    jt, cluster_marks = collect_evidence(arbitrary_cluster_ind, cluster_marks, jt)
    cluster_marks = [true for k in keys(jt.vertex_to_cluster)]
    jt, cluster_marks = distribute_evidence(arbitrary_cluster_ind, cluster_marks, jt)
    return jt
end
