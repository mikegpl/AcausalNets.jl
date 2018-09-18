#=
assignments:
- Julia version: 
- Author: marcin
- Date: 2018-09-01
=#

import AcausalNets.Algebra:
    event

import AcausalNets.Common:
    is_subset,
    ncategories


import AcausalNets.Systems:
    DiscreteSystem,
    relevant_variables,
    distribution,
    multiply_kron,
    permute_distribution,
    parents,
    variables,
    merge_systems,
    sum_distribution

import AcausalNets.Inference:
    JoinTree,
    ParentCliquesDict,
    shallowcopy




const Evidence{D} = DiscreteSystem{D}

function Evidence(variables::Vector{Variable}, distribution::D) where D
    isapprox(sum_distribution(distribution), 1) ||
        error("Sum of the distribution must be one, but is $(sum_distribution(Distribution))!")
    DiscreteSystem{D}(variables, distribution)
end


function apply_evidence(system::DiscreteSystem{D1}, evidence::Evidence{D2}) where {D1, D2 <: D1}
    system = shallowcopy(system)
    sys_vars = relevant_variables(system)
    ev_vars = relevant_variables(evidence)
    is_subset(Set(ev_vars), Set(sys_vars)) ||
        error("variables from outside the system in evidence!")

    ev_dimensions = [ncategories(v) for v in ev_vars]
    ev_dist = distribution(evidence)

    non_ev_vars = [v for v in sys_vars if !(v in ev_vars)]
    non_ev_dimensions = [ncategories(v) for v in non_ev_vars]
    non_ev_size = prod(non_ev_dimensions)

    ev_dist = multiply_kron(ev_dist, identity_distribution(D1, non_ev_size))
    ev_dimensions = vcat(ev_dimensions, non_ev_dimensions)
    ev_vars = vcat(ev_vars, non_ev_vars)
    ev_dist = permute_distribution(
                ev_dist,
                ev_dimensions,
                [findfirst(var -> var==v, ev_vars) for v in sys_vars]
            )
    DiscreteSystem{D1}(
        parents(system),
        variables(system),
        event(distribution(system), ev_dist)
    )
end

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