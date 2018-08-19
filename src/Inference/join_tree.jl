#=
join_tree:
- Julia version: 
- Author: marcin
- Date: 2018-08-19
=#
using LightGraphs
using LinearAlgebra
using QI

import AcausalNets.Systems:
    DiscreteSystem,
    ncategories,
    variables

import AcausalNets.Structures:
    DiscreteBayesNet,
    parent_systems


import AcausalNets.Inference:
    SystemCluster

import AcausalNets.Algebra:
    eye,
    event, star, unstar


struct JoinTree
    graph::Graph
    clusters::Vector{SystemCluster}
    vertex_to_num::Dict{Int64, AbstractMatrix}
    edge_to_num::Dict{Set{Int}, AbstractMatrix}
end

sepset_comparator(c1::SystemCluster, c2::SystemCluster)::Int64 = -length(intersect(c1, c2))
cluster_size(cluster::SystemCluster)::Int64 = prod([ncategories(s) for s in cluster])

function JoinTree(clusters::Vector{SystemCluster}, dbn::DiscreteBayesNet)::JoinTree
    candidate_sepsets = []
    trees = Dict([c => c for c in clusters])
    chosen_sepsets = Set()

    jt = JoinTree(
        Graph(length(clusters)),
        clusters,
        Dict([v => eye(cluster_size(clusters[v])) for v=1:length(clusters)]),
        Dict()
    )

    for c1 = 1:length(clusters)
        for c2= 1:length(clusters)
            if c1 != c2
                push!(candidate_sepsets,(c1, c2))
            end
        end
    end
    candidate_sepsets = sort(candidate_sepsets, by=c -> sepset_comparator(clusters[c[1]], clusters[c[2]]))
    i = 1
    n = length(clusters)

    while length(chosen_sepsets) < n-1
        i1, i2 = candidate_sepsets[i]
        c1, c2 = clusters[i1], clusters[i2]
        sepset = intersect(c1, c2)
        if (trees[c1] != trees[c2]) && !any([sepset==s for s in chosen_sepsets])
            push!(chosen_sepsets, sepset)
            trees[c1] = trees[c2] = union(c1, c2)
            add_edge!(jt.graph, i1, i2)
            push!(jt.edge_to_num, Set([i1, i2]) => eye(cluster_size(sepset)))
        end
        i += 1
    end
    jt
end

shallowcopy(jt::JoinTree)::JoinTree = JoinTree(jt.graph, jt.clusters, deepcopy(jt.vertex_to_num), deepcopy(jt.edge_to_num))

family(ds::S, dbn::DiscreteBayesNet{S}) where {S} = union(Set([ds]), Set(parent_systems(ds, dbn)))


function initialize(jt::JoinTree, dbn::DiscreteBayesNet, assignments::Dict=Dict())::JoinTree
    jt = shallowcopy(jt)
    for v1 in systems(dbn)
        parent_cluster_ind = [c for c=1:length(jt.clusters) if is_subset(family(v1, dbn), Set(jt.clusters[c]))][1]
        mul_elem = eye(1)
        multiplied_indices = []
        multiplied_dimensions = []

        for v2 in [v for v in systems(dbn) if v in jt.clusters[parent_cluster_ind]]
            factor = eye(1)
            if v1 == v2
                assignment_factor = eye(1)
                for p in parent_systems(v2, dbn)
                    p_index = system_to_node(p, dbn)
                    push!(multiplied_indices, p_index)
                    push!(multiplied_dimensions, ncategories(p))
                    assignment_factor = kron(assignment_factor, eye(ncategories(p)))
                end
                for v in variables(v2)
                    assignment = get(assignments, v, eye(ncategories(v)))
                    assignment_factor = kron(assignment_factor, assignment)
                end
                factor = event(v2.distribution, assignment_factor) * tr(v2.distribution)
                push!(multiplied_indices, system_to_node(v2, dbn))
                push!(multiplied_dimensions, ncategories(v2))
            elseif !in(v2, parent_systems(v1, dbn))
                factor = eye(ncategories(v2))
                push!(multiplied_indices, system_to_node(v2, dbn))
                push!(multiplied_dimensions, ncategories(v2))
            end
            mul_elem = kron(mul_elem, factor)
        end
        sorted_multiplied_indices = sort(multiplied_indices)
        right_order = [findfirst([m == mul for mul in multiplied_indices]) for m in sorted_multiplied_indices]

        multiplied_dimensions = [multiplied_dimensions[r] for r in right_order]

        previous_init = jt.vertex_to_num[parent_cluster_ind]
        jt.vertex_to_num[parent_cluster_ind] = star(mul_elem, previous_init)
    end
    jt
end


