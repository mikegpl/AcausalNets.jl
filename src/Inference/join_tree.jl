#=
join_tree:
- Julia version: 1.0
- Author: marcin
- Date: 2018-08-19
=#
using LightGraphs
using LinearAlgebra
using QuantumInformation

import AcausalNets.Common:
    Variable,
    is_subset

import AcausalNets.Systems:
    DiscreteSystem,
    ncategories,
    variables,
    parents,
    distribution,
    shallowcopy,
    merge_systems,
    identity_system,
    identity_distribution, permute_distribution, reduce_distribution,
    multiply_star, divide_star, multiply_kron # those methods should be system-agnostic

import AcausalNets.Structures:
    DiscreteBayesNet,
    parent_systems,
    family,
    systems,
    system_to_node,
    variable_to_node

import AcausalNets.Inference:
    Evidence,
    apply_evidence


"""
The Join Tree structure, as described in http://pages.cs.wisc.edu/~dpage/ijar95.pdf
both vertices and edges ("sepsets") are mapped to a DiscreteSystem, but we only use the systems
in vertices for belief propagation.
Perhaps systems with mutual proabililites should go to the edges?
"""
struct JoinTree{S <: DiscreteSystem}
    graph::Graph
    vertex_to_cluster::Dict{Int64, S}
    edge_to_sepset::Dict{Set{Int}, S}
end

"""
performs steps described in sections 4.1 - 5.2 with 6.3
of http://pages.cs.wisc.edu/~dpage/ijar95.pdf
"""
function unpropagated_join_tree(dbn::DiscreteBayesNet{S},
        vars_to_infer::Vector{Variable},
        observations::Vector{E} = E[]
        ) where {
            D1,
            D2 <: D1,
            S <: DiscreteSystem{D1},
            E <: Evidence{D2}
        }
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
    observations_jt
end

"""
Builds a Join Tree
http://pages.cs.wisc.edu/~dpage/ijar95.pdf, section 4.4
"""
function JoinTree(cliques::Vector{Vector{S}}, dbn::DiscreteBayesNet{S})::JoinTree{S} where S
    candidate_sepsets = []
    trees = Dict([c => Set([c]) for c in cliques])
    chosen_sepsets = Set()

    parent_cliques = parent_cliques_dict(cliques, dbn)
    jt = JoinTree{S}(
        Graph(length(cliques)),
        Dict([
            i => Cluster([sys_or_id(sys, cliques[i], parent_cliques) for sys in cliques[i]])
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
    created_edges = 0
    while created_edges < n-1
        i1, i2 = candidate_sepsets[i]
        c1, c2 = cliques[i1], cliques[i2]
        sepset = [s for s in systems(dbn) if s in intersect(c1, c2)]
        if isempty(intersect(trees[c1], trees[c2]))
            push!(chosen_sepsets, sepset)
            union_tree = union(trees[c1], trees[c2])
            [trees[tree_elem] = union_tree for  tree_elem in union_tree]
            add_edge!(jt.graph, i1, i2)
            sepset_variables = reduce(vcat, [variables(s) for s in sepset])
            push!(jt.edge_to_sepset, Set([i1, i2]) => identity_system(S, sepset_variables))
            created_edges += 1
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

function sepset_comparator(c1::Vector{S}, c2::Vector{S})::Int64 where S
    -length(intersect(c1, c2))
end

const ParentCliquesDict{S} = Dict{S, Vector{S}}

"""Chooses a parent clique in the join tree for each system in Bayes net
http://pages.cs.wisc.edu/~dpage/ijar95.pdf, section 5.2
"""
function parent_cliques_dict(cliques::Vector{Vector{S}}, dbn::DiscreteBayesNet{S})::ParentCliquesDict{S} where S
    Dict([
        sys => first([
            c for c in cliques if is_subset(family(sys, dbn), Set(c))
            ])
        for sys in systems(dbn)
    ])
end

"""A uitility function which returns the system if the clique has been designated as its parent clique
or a system with the same variables, but with an identity matrix as density matrix
if the clique is not the parent clique"""
function sys_or_id(system::S, clique::Vector{S}, parent_cliques::ParentCliquesDict{S})::S where S
    if parent_cliques[system] == clique
        return system
    else
        return identity_system(system)
    end
end

"""Builds a Join Tree vertex from a cluster (clique) of variables"""
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
            error("parents of system $variables(sys) outside cluster!")
    end
    merge_systems(systems)
end

"""Normalizes density matrices of vertices and edges in the JoinTree, so that their traces are 1"""
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

const MoralGraph = Graph


"""
moralization of the graph
http://pages.cs.wisc.edu/~dpage/ijar95.pdf
(section 4.1`)
"""
function moral_graph(dbn::DiscreteBayesNet)::MoralGraph
    result = MoralGraph(deepcopy(dbn.dag))
    for sys in systems(dbn)
        for p1 in parent_systems(sys, dbn)
            for p2 in parent_systems(sys, dbn)
                p1_ind = system_to_node(p1, dbn)
                p2_ind = system_to_node(p2, dbn)
                if(p1_ind!=p2_ind)
                    add_edge!(result, p1_ind, p2_ind)
                end
            end
        end
    end
    return result
end


function enforce_clique_for_quantum_systes(dbn::DiscreteBayesNet, mg::MoralGraph)
    for sys in systems(dbn)
        if(!isdiag(sys.distribution))
            quantum_node = system_to_node(sys, dbn)
            clique = union(all_neighbors(dbn.dag, system_to_node(sys, dbn)), quantum_node)

            for c1 in clique
                for c2 in clique
                    if c1 > c2
                        add_edge!(mg, c1, c2)
                    end
                end
            end
        end
    end
end


"""
step which is not specified in the article an introduced by us
after moralization we make sure there is a clique of systems we want to infer,
as this clique will later become one of the vertices in the junction tree
since we can only infer a system of variables which are in the same vertex of the join tree
we need to make sure all of the variables which interest us will end up in the same vertex
note: this can result in a sub-optimal join tree
"""
function enforce_clique(dbn::DiscreteBayesNet, mg::MoralGraph, vars_to_infer::Vector{Variable}=Variable[])::MoralGraph
    mg_enforced = deepcopy(mg)
    inferred_nodes = [variable_to_node(v, dbn) for v in vars_to_infer]
    for n1 in inferred_nodes
        for n2 in inferred_nodes
            if n1 != n2
                add_edge!(mg_enforced, n1, n2)
            end
        end
    end

    enforce_clique_for_quantum_systes(dbn, mg_enforced)

    return mg_enforced
end

const TriangulatedGraph = Graph

"""
Triangulation of the graph and identification of cliques in triangulted graph
http://pages.cs.wisc.edu/~dpage/ijar95.pdf
(section 4.2 and 4.3)
"""
function triangulate(mg::MoralGraph, dbn::DiscreteBayesNet{S})::Tuple{TriangulatedGraph, Vector{Vector{S}}} where S
    mg_copy = [false for _ in vertices(mg)]
    mg = deepcopy(mg)
    nl = systems(dbn)
    cliques = Set[]
    while(!all(mg_copy))
        least_edges_to_be_added = Inf
        chosen_vertex = 0
        chosen_clique = Set()
        for v=1:length(mg_copy)
            if mg_copy[v]
                continue
            else
                clique = Set()
                for e in edges(mg)
                    if (v==src(e) || v==dst(e)) && !mg_copy[src(e)] && !mg_copy[dst(e)]
                        push!(clique, src(e))
                        push!(clique, dst(e))
                    end
                end
                edges_todo = 0
                for v1 in clique
                    for v2 in clique
                        if v1 != v2 && !in((v1, v2), edges(mg))
                            edges_todo +=1
                        end
                    end
                end
                edges_todo /= 2

                if edges_todo < least_edges_to_be_added ||
                    ((edges_todo == least_edges_to_be_added) &&
                        (prod([ncategories(systems(dbn)[n]) for n in clique]) <= prod([ncategories(systems(dbn)[n]) for n in chosen_clique])))
                    least_edges_to_be_added = edges_todo
                    chosen_vertex = v
                    chosen_clique = clique
                end
            end
        end
        chosen_nodes = Set([nl[n] for n in chosen_clique])
        if !any([is_subset(chosen_nodes, clique) for clique in cliques])
            push!(cliques, chosen_nodes)
        end
        mg_copy[chosen_vertex] = true
        for v1 in chosen_clique
            for v2 in chosen_clique
                if v1 != v2 && !in((v1, v2), edges(mg))
                    add_edge!(mg, v1, v2)
                end
            end
        end
    end
    cliques = [sort([c for c in clique], by=c -> system_to_node(c, dbn)) for clique in cliques]
    return mg, cliques
end


"""
Applies observations to a Join Tree, whose vertices have been initialized
Observations about variables are applied only to the vertices which are initialized with probabilities
of those variables (so evidence of each variable is applied only to the vertex which had been designated as its "parent cliuqe"
"""
function apply_observations(
        jt::JoinTree{S},
        parent_cliques_dict::ParentCliquesDict{S},
        observations::Vector{E}
        ) where {D1, D2 <: D1, S <: DiscreteSystem{D1}, E <: Evidence{D2} }
    observations_dict = Dict{Int, E}()
    for v in keys(jt.vertex_to_cluster)
        cluster = jt.vertex_to_cluster[v]
        child_systems = [
            sys for (sys, par_cliq) in parent_cliques_dict
            if Set(variables(cluster)) == Set(vcat([variables(s) for s in par_cliq]...))
            ]
        relevant_observations = E[
            o for o in observations if
            any([
                    is_subset(Set(variables(o)), Set(variables(sys)))
                    for sys in child_systems
                ])
            ]
        observations_dict[v] = merge_systems(relevant_observations)
    end

    new_vertex_to_cluster = Dict([
            v => apply_evidence(jt.vertex_to_cluster[v], observations_dict[v])
            for v in keys(jt.vertex_to_cluster)
                            ])
    JoinTree{S}(jt.graph, new_vertex_to_cluster, jt.edge_to_sepset)
end
