import AcausalNets.Common:
    is_subset

using AcausalNets.Common

struct DiscreteSystem{D}
    parents         ::Vector{Variable}
    variables       ::Vector{Variable}
    distribution    ::D

    function DiscreteSystem{D}(
        parents::Vector{Variable},
        variables::Vector{Variable},
        distribution::D
        ) where D
        check_distribution(distribution, parents, variables) ||
            error("Dimensions of parents/variables/distribution not matching!")
        new(parents, variables, distribution)
    end
end

# Base.:(==)(ds1::DiscreteSystem{D}, ds2::DiscreteSystem{D}) where D =
#     parents(ds1) == parents(ds2) &&
#     variables(ds1) == variables(ds2)

DiscreteSystem{D}(variables::Vector{Variable}, distribution::D) where D = DiscreteSystem{D}(Variable[], variables, distribution)


parents(ds::DiscreteSystem) = ds.parents
variables(ds::DiscreteSystem) = ds.variables
distribution(ds::DiscreteSystem) = ds.distribution

parents_names(ds::DiscreteSystem) = [p.name for p in parents(ds)]
variables_names(ds::DiscreteSystem) = [v.name for v in variables(ds)]
function Base.string(ds::DiscreteSystem, verbose::Bool=true)
    result_str = join([string(n) for n in variables_names(ds)], ",")
    if verbose
        parents_str = join([String(p.name) for p in parents(ds)], ",")
        result_str = string(result_str, "|", parents_str)
    end
    return result_str
end

relevant_variables(system::DiscreteSystem) = vcat(system.parents, system.variables)


Common.ncategories(ds::DiscreteSystem) = ncategories(variables(ds))

# function enforce_parents_order(ds::DiscreteSystem{D}, existing_variables::Vector{Variable}) where D
#     new_parents_order = [v for v in existing_variables if v in parents(ds)]
#     old_parents_order = parents(ds)
#
#     new_parents_indexing = Vector{Int64}([findfirst([p == o for o in old_parents_order]) for p in new_parents_order])
#     new_variables_indexing = Vector(1:length(variables(ds)))
#     permute_system(ds, new_parents_indexing, new_variables_indexing)
# end

function expand_parents(ds::S, existing_systems::Vector{S})::S where {D, S <: DiscreteSystem{D}}
    all_parents = Vector{Variable}(vcat([
        variables(sys) for sys in existing_systems
        if !isempty(intersect(parents(ds), variables(sys)))
    ]...))

    S(
        all_parents,
        variables(ds),
        distribution(
            sub_system(ds, vcat(all_parents, variables(ds)))
        )
    )
end

function is_parent(potential_parent::DiscreteSystem{D}, potential_child::DiscreteSystem{D}) where D
    par_variables = Set(variables(potential_parent))
    child_parents = Set(parents(potential_child))

    is_subset(par_variables, child_parents) ||
        length(intersect(par_variables, child_parents)) == 0 ||
        error("Only some of the variables in parent are listed as child's parents.
                Inputting the systems into an AcausalNet may fix that problem.")
end



function merge_systems(systems::Vector{S}, verbose::Bool = true)::S where {D, S <: DiscreteSystem{D}}
    """
    We assume that systems are sorted topologically (parents first)
    We adopt the ordering defined in
    https://arxiv.org/pdf/0708.1337.pdf
    (equation 100)
    """
    all_variables = Variable[]
    all_parents = Variable[]
    for s in systems
        for v in variables(s)
            !(v in all_variables) || error("$v duplicated!")
            push!(all_variables, v)
        end
    end
    for s in systems
        for p in parents(s)
            if !(p in all_variables) && (p in all_parents)
                push!(all_parents, p)
            end
        end
    end
    all_relevant_variables = vcat(all_parents, all_variables)
    target_size = prod([ncategories(v) for v in all_relevant_variables])
    result_distribution = identity_distribution(D, target_size)

    # order of multiplication defined in
    # https://arxiv.org/pdf/0708.1337.pdf
    # (100)
    # A1 * A2 * A3 = ((A1 * A2) * A3)
    debug_str = "Id($target_size)"
    for sys in systems
        ordered_distribution = distribution(sub_system(sys, all_relevant_variables))
        result_distribution = multiply_star(result_distribution, ordered_distribution)
        debug_str = string(
            "( ",
            debug_str,
            " * ro",
            string(sys,true),
            " )"
        )
    end

    if verbose
        println(debug_str)
    end
    S(all_parents, all_variables, result_distribution)
end

function identity_system(system::DiscreteSystem{D})::DiscreteSystem{D} where D
    return DiscreteSystem{D}(variables(system), identity_distribution(D, ncategories(system)))
end

function identity_system(::Type{DiscreteSystem{D}}, variables::Vector{Variable})::DiscreteSystem{D} where D
    target_size = prod([ncategories(v) for v in variables])
    DiscreteSystem{D}(variables, identity_distribution(D, target_size))
end

shallowcopy(ds::DiscreteSystem{D}) where D = DiscreteSystem{D}(
    parents(ds),
    variables(ds),
    deepcopy(distribution(ds))
    )

# not to confuse with QI's permute_systems - this is a higher-level implementation
function permute_system(ds::DiscreteSystem{D}, new_parent_indexing::Vector{Int64}, new_variable_indexing::Vector{Int64}) where D
    new_indexing = vcat(new_parent_indexing, new_variable_indexing .+ length(parents(ds)))
    dimensions = [ncategories(v) for v in relevant_variables(ds)]
    new_distribution = permute_distribution(distribution(ds), dimensions, new_indexing)
    new_parents = Variable[parents(ds)[i] for i in new_parent_indexing]
    new_variables = Variable[variables(ds)[i] for i in new_variable_indexing]
    DiscreteSystem{D}(new_parents, new_variables, new_distribution)
end

function permute_system(ds::DiscreteSystem{D}, new_variable_indexing::Vector{Int64}) where D
    length(parents(ds)) == 0 || error("Parents permutation must be specified for this system!")
    permute_system(ds, Int64[], new_variable_indexing)
end

function sub_system(ds::DiscreteSystem{D}, desired_variables::Vector{Variable}) where D
    sys_vars = Variable[v for v in relevant_variables(ds) if v in desired_variables]
    non_sys_vars = Variable[v for v in desired_variables if !(v in sys_vars)]
    unordered_vars = vcat(sys_vars, non_sys_vars)

    redundant_indices = Vector{Int}(findall(
        (v) -> !(v in unordered_vars),
        relevant_variables(ds)
    ))

    sys_dist = reduce_distribution(
                    distribution(ds),
                    [ncategories(v) for v in relevant_variables(ds)],
                    redundant_indices
            )
    non_sys_dist = identity_distribution(D, prod([ncategories(v) for v in non_sys_vars]))
    unordered_dist = multiply_kron(sys_dist, non_sys_dist)

    ordered_indices = Int[
            findfirst(
                (v) -> v == var,
                unordered_vars
            )
            for var in desired_variables
        ]
    ordered_dist = permute_distribution(
                    unordered_dist,
                    Int[ncategories(v) for v in unordered_vars],
                    ordered_indices
            )
    DiscreteSystem{D}(desired_variables, ordered_dist)
end