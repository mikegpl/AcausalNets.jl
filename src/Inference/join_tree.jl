#=
join_tree:
- Julia version: 
- Author: marcin
- Date: 2018-08-19
=#
using LightGraphs
using LinearAlgebra
using QI

import AcausalNets.Common:
    Variable

import AcausalNets.Systems:
    DiscreteSystem,
    ncategories,
    variables,
    parents,
    distribution,
    shallowcopy,
    merge_systems,
    identity_system,
    sub_system,
    identity_distribution, permute_distribution, reduce_distribution,
    multiply_star, divide_star, multiply_kron # those methods should be system-agnostic

import AcausalNets.Structures:
    DiscreteBayesNet,
    parent_systems,
    family


function sepset_comparator(c1::Vector{S}, c2::Vector{S})::Int64 where S
    -length(intersect(c1, c2))
end

const ParentCliquesDict{S} = Dict{S, Vector{S}}

function parent_cliques_dict(cliques::Vector{Vector{S}}, dbn::DiscreteBayesNet{S})::ParentCliquesDict{S} where S
    Dict([
        sys => first([
            c for c in cliques if is_subset(family(sys, dbn), Set(c))
            ])
        for sys in systems(dbn)
    ])
end

function sys_or_id(system::S, clique::Vector{S}, parent_cliques::ParentCliquesDict{S})::S where S
    if parent_cliques[system] == clique
        return system
    else
        return identity_system(system)
    end
end

function Cluster(systems::Vector{S})::S where S
    all_variables = Variable[]
    all_parents = Variable[]
    for s in systems
        for v in variables(s)
            !(v in all_variables) || error("$v duplicated!")
            push!(all_variables, v)
        end
    end
    for s in systems
        all([p in all_variables for p in parents(s)]) ||
            error("parents of system $(variables(s)) outside cluster!")
    end
    merge_systems(systems)
#     merge_systems([
#         sub_system(
#             s,
#             [v for v in relevant_variables(s) if v in all_variables]
#         )
#         for s in systems
#     ])
end

struct JoinTree{S <: DiscreteSystem}
    graph::Graph
    vertex_to_cluster::Dict{Int64, S}
    edge_to_sepset::Dict{Set{Int}, S}
end

function JoinTree(cliques::Vector{Vector{S}}, dbn::DiscreteBayesNet{S})::JoinTree{S} where S
    candidate_sepsets = []
    trees = Dict([c => c for c in cliques])
    chosen_sepsets = Set()

    parent_cliques = parent_cliques_dict(cliques, dbn)
    jt = JoinTree{S}(
        Graph(length(cliques)),
        Dict([
            i => Cluster([sys_or_id(sys, cliques[i], parent_cliques) for sys in cliques[i]])
#             i => Cluster([sys for sys in cliques[i]])
            for i in 1:length(cliques)
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
            sepset_variables = reduce(vcat, [variables(s) for s in sepset])
#
#             # mutual probability of c1:c2 as defined in Leifer
#             c1_and_c2_sorted = [s for s in systems(dbn) if s in union(c1, c2)]
#
#             c1Uc2 = Cluster(c1_and_c2_sorted)
#
#             println(string(c1Uc2), " ", sepset_variables)
#             c1_kronned = sub_system(Cluster([s for s in systems(dbn) if s in c1]), variables(c1Uc2))
#             c2_kronned = sub_system(Cluster([s for s in systems(dbn) if s in c2]), variables(c1Uc2))
#
#             c1_mut_c2_dist = divide_star(
#                 distribution(c1Uc2),
#                 multiply_star(distribution(c1_kronned), distribution(c2_kronned))
#             )

            push!(jt.edge_to_sepset, Set([i1, i2]) => identity_system(S, sepset_variables))
#             push!(jt.edge_to_sepset, Set([i1, i2]) => S(variables(c1Uc2), c1_mut_c2_dist))
        end
        i += 1
    end
    jt
end

shallowcopy(jt::JoinTree) = JoinTree(
        jt.graph,
        Dict([
            k => shallowcopy(jt.vertex_to_cluster[k]) for k in keys(jt.vertex_to_cluster)
            ]),
        Dict([
            k => shallowcopy(jt.edge_to_sepset[k]) for k in keys(jt.edge_to_sepset)
            ])
    )

function normalize(jt::JoinTree{S}) where S
    JoinTree(
        jt.graph,
        Dict([
                v => S(
                    c.parents,
                    c.variables,
                    c.distribution / sum_distribution(c.distribution)
                    )
                for (v, c) in jt.vertex_to_cluster
            ]),
        Dict([
                e => S(
                    s.parents,
                    s.variables,
                    s.distribution / sum_distribution(s.distribution)
                    )
                for (e, s) in jt.edge_to_sepset
            ])
    )
end