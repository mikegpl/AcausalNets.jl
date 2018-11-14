#=
assignments:
- Julia version: 1.0
- Author: marcin
- Date: 2018-09-01
=#

import AcausalNets.Algebra:
    event

import AcausalNets.Common:
    Variable,
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
    sum_distribution,
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

    evidence_system = sub_system(evidence, sys_vars)
#     ev_dimensions = [ncategories(v) for v in ev_vars]
#     ev_dist = distribution(evidence)
#
#     non_ev_vars = [v for v in sys_vars if !(v in ev_vars)]
#     non_ev_dimensions = [ncategories(v) for v in non_ev_vars]
#     non_ev_size = prod(non_ev_dimensions)
#
#     ev_dist = multiply_kron(ev_dist, identity_distribution(D1, non_ev_size))
#     ev_dimensions = vcat(ev_dimensions, non_ev_dimensions)
#     ev_vars = vcat(ev_vars, non_ev_vars)
#     ev_dist = permute_distribution(
#                 ev_dist,
#                 ev_dimensions,
#                 [findfirst(var -> var==v, ev_vars) for v in sys_vars]
#             )
    DiscreteSystem{D1}(
        parents(system),
        variables(system),
        event(distribution(system), distribution(evidence_system))
    )
end