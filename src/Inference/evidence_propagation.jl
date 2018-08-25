#=
join_tree_messaging:
- Julia version: 
- Author: marcin
- Date: 2018-08-19
=#

using LightGraphs
using QI

import AcausalNets.Systems:
    ncategories

import AcausalNets.Algebra:
    star, unstar

import AcausalNets.Structures:
    DiscreteBayesNet

import AcausalNets.Inference:
    JoinTree,
    shallowcopy

function single_message_pass(from_ind::Int, to_ind::Int, jt::JoinTree{S}, dbn::DiscreteBayesNet{S})::JoinTree{S} where S
    jt = shallowcopy(jt)
    if (from_ind, to_ind) in edges(jt.graph)
        cluster_from = jt.clusters[from_ind]
        cluster_to = jt.clusters[to_ind]
        sepset = intersect(cluster_from, cluster_to)
        to_trace_out_sym = setdiff(cluster_from, sepset)
        to_trace_out_ind = [findfirst([s==sym for sym in cluster_from]) for s in to_trace_out_sym]
        from_variables_sizes = [ncategories(v) for v in cluster_from]
        cluster_from_num = jt.vertex_to_num[from_ind]
        old_sepset_num = jt.edge_to_num[Set([from_ind, to_ind])]
        new_sepset_num = ptrace(cluster_from_num, from_variables_sizes, to_trace_out_ind)

        jt.edge_to_num[Set([from_ind, to_ind])] = new_sepset_num

        cluster_to_num = jt.vertex_to_num[to_ind]

        message = unstar(new_sepset_num, old_sepset_num)

        message_sym = Vector(sepset)
        for v in cluster_to
            if !in(v, message_sym)
                push!(message_sym, v)
                message = kron(message, eye(ncategories(v)))
            end
        end
        message_sorted_inds = [findfirst([s==sym for sym in message_sym]) for s in systems(dbn) if s in message_sym]
        message_dims = [ncategories(s) for s in message_sym]
        message_sorted = permute_systems(message, message_dims, message_sorted_inds)
        jt.vertex_to_num[to_ind] = star(cluster_to_num, message_sorted)
    end
    return jt
end

function collect_evidence(cluster_ind::Int, cluster_marks::Vector{Bool}, jt::JoinTree{S}, dbn::DiscreteBayesNet{S})::Tuple{JoinTree{S}, Vector{Bool}} where S
    jt = shallowcopy(jt)

    cluster_marks[cluster_ind] = false
    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt, cluster_marks = collect_evidence(neighbor, cluster_marks, jt, dbn)
            jt = single_message_pass(neighbor, cluster_ind, jt, dbn)
        end

    end
    jt, cluster_marks

end

function distribute_evidence(cluster_ind::Int, cluster_marks::Vector{Bool}, jt::JoinTree{S}, dbn::DiscreteBayesNet{S})::Tuple{JoinTree{S}, Vector{Bool}} where S
    cluster_marks[cluster_ind] = false
    jt = shallowcopy(jt)

    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt = single_message_pass(cluster_ind, neighbor, jt, dbn)
        end
    end
    for neighbor in neighbors(jt.graph, cluster_ind)
        if cluster_marks[neighbor]
            jt, cluster_mars = distribute_evidence(neighbor, cluster_marks, jt, dbn)
        end
    end
    jt, cluster_marks
end

function global_propagation(jt::JoinTree{S}, dbn::DiscreteBayesNet{S})::JoinTree{S} where S
    jt = shallowcopy(jt)
    cluster_marks = [true for c in jt.clusters]
    arbitrary_cluster_ind = 1
    jt, cluster_marks = collect_evidence(arbitrary_cluster_ind, cluster_marks, jt, dbn)
    cluster_marks = [true for c in jt.clusters]
    jt, cluster_marks = distribute_evidence(arbitrary_cluster_ind, cluster_marks, jt, dbn)
    return jt
end

function normalization(jt::JoinTree{S})::JoinTree{S} where S
    jt = shallowcopy(jt)
    for key in keys(jt.vertex_to_num)
        jt.vertex_to_num[key] = normalize(jt.vertex_to_num[key])
    end
    for key in keys(jt.edge_to_num)
        jt.edge_to_num[key] = normalize(jt.edge_to_num[key])
    end
end