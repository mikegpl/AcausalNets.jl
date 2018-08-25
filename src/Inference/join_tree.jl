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
    variables,
    identity_distribution

import AcausalNets.Structures:
    DiscreteBayesNet,
    parent_systems


import AcausalNets.Algebra:
    eye,
    event, star, unstar


struct JoinTree{S <: DiscreteSystem}
    graph::Graph
    vertex_to_cluster::Dict{Int64, S}
    edge_to_sepset::Dict{Set{Int}, S}
end

function sepset_comparator(c1::Vector{S}, c2::Vector{S})::Int64 where S
    -length(intersect(c1, c2))
end

function clique_size(clique::Vector{S})::Int64 where S
    prod([ncategories(s) for s in clique])
end


function JoinTree(cliques::Vector{Vector{S}}, dbn::DiscreteBayesNet{S})::JoinTree{S} where S
    candidate_sepsets = []
    trees = Dict([c => c for c in cliques])
    chosen_sepsets = Set()

    jt = JoinTree{S}(
        Graph(length(cliques)),
        Dict([
            v = # systems made of merged cliques - effectively initialized in this step (without assignments yet)#> S(,identity_distribution(S, cluster_size(clusters[v])))
            for v=1:length(cliques)
            ]),
        Dict()
    )

    for c1 = 1:length(cliques)
        for c2= 1:length(cliques)
            if c1 != c2
                push!(candidate_sepsets,(c1, c2))
            end
        end
    end

    candidate_sepsets = sort(candidate_sepsets, by=c -> sepset_comparator(cliques[c[1]], cliques[c[2]]))
    i = 1
    n = length(cliques)

    while length(chosen_sepsets) < n-1
        i1, i2 = candidate_sepsets[i]
        c1, c2 = cliques[i1], cliques[i2]
        sepset = intersect(c1, c2)
        if (trees[c1] != trees[c2]) && !any([sepset==s for s in chosen_sepsets])
            push!(chosen_sepsets, sepset)
            trees[c1] = trees[c2] = union(c1, c2)
            add_edge!(jt.graph, i1, i2)
            push!(jt.edge_to_num, Set([i1, i2]) => eye(clique_size(sepset)))
        end
        i += 1
    end
    jt
end
#
# shallowcopy(jt::JoinTree)::JoinTree = JoinTree(jt.graph, jt.clusters, deepcopy(jt.vertex_to_num), deepcopy(jt.edge_to_num))
#
# family(ds::S, dbn::DiscreteBayesNet{S}) where {S} = union(Set([ds]), Set(parent_systems(ds, dbn)))
#
#
# function initialize(jt::JoinTree{S}, dbn::DiscreteBayesNet{S}, assignments::Vector{S}=S[])::JoinTree{S} where S
#     length(assignments) == 0 || error("Assignments unimplemented!")
#     jt = shallowcopy(jt)
#     for v1 in systems(dbn)
#         parent_cluster_ind = [c for c=1:length(jt.clusters) if is_subset(family(v1, dbn), Set(jt.clusters[c]))][1]
#         mul_elem = eye(1)
#         multiplied_indices = []
#         multiplied_dimensions = []
#
#         for v2 in [v for v in jt.clusters[parent_cluster_ind]]
#             factor = eye(1)
#             if v1 == v2
#                 assignment_factor = eye(1)
#                 for p in parent_systems(v2, dbn)
#                     p_index = system_to_node(p, dbn)
#                     push!(multiplied_indices, p_index)
#                     push!(multiplied_dimensions, ncategories(p))
#                     assignment_factor = kron(assignment_factor, eye(ncategories(p)))
#                 end
#                 #TODO Assignments
#                 for v in variables(v2)
#                     assignment = eye(ncategories(v))
#                     assignment_factor = kron(assignment_factor, assignment)
#                 end
#                 factor = event(v2.distribution, assignment_factor) * tr(v2.distribution)
#                 push!(multiplied_indices, system_to_node(v2, dbn))
#                 push!(multiplied_dimensions, ncategories(v2))
#             elseif !in(v2, parent_systems(v1, dbn))
#                 factor = eye(ncategories(v2))
#                 push!(multiplied_indices, system_to_node(v2, dbn))
#                 push!(multiplied_dimensions, ncategories(v2))
#             end
#             mul_elem = kron(mul_elem, factor)
#         end
#         sorted_multiplied_indices = sort(multiplied_indices)
#         right_order = [findfirst([m == mul for mul in multiplied_indices]) for m in sorted_multiplied_indices]
#
#         multiplied_dimensions = [multiplied_dimensions[r] for r in right_order]
#
#         previous_init = jt.vertex_to_num[parent_cluster_ind]
#         jt.vertex_to_num[parent_cluster_ind] = star(mul_elem, previous_init)
#     end
#     jt
# end


