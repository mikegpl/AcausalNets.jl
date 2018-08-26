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
Base.string(ds::DiscreteSystem) = join([String(n) for n in variables_names(ds)], ", ")

relevant_variables(system::DiscreteSystem) = vcat(system.parents, system.variables)


Common.ncategories(ds::DiscreteSystem) = ncategories(variables(ds))

function enforce_parents_order(ds::DiscreteSystem{D}, existing_variables::Vector{Variable}) where D
    new_parents_order = [v for v in existing_variables if v in parents(ds)]
    old_parents_order = parents(ds)

    new_parents_indexing = Vector{Int64}([findfirst([p == o for o in old_parents_order]) for p in new_parents_order])
    new_variables_indexing = Vector(1:length(variables(ds)))
    permute_system(ds, new_parents_indexing, new_variables_indexing)
end

function expand_parents(ds::DiscreteSystem{D}, existing_systems::Vector{DiscreteSystem{D}}) where D
    parent_systems = [sys for sys in existing_systems if any([v in parents(ds) for v in variables(sys)])]
    for sys in parent_systems
        for var in sys.variables
            if !(var in parents(ds))
                ds = prepend_parent(ds, var)
            end
        end
    end
    ds
end

function is_parent(potential_parent::DiscreteSystem{D}, potential_child::DiscreteSystem{D}) where D
    par_variables = Set(variables(potential_parent))
    child_parents = Set(parents(potential_child))

    is_subset(par_variables, child_parents) ||
        length(intersect(par_variables, child_parents)) == 0 ||
        error("Only some of the variables in parent are listed as child's parents.
                Inputting the systems into an AcausalNet may fix that problem.")
end



function merge_systems(systems::Vector{S})::S where {D, S <: DiscreteSystem{D}}
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

    for sys in systems
        sys_distribution = distribution(sys)
        sys_variables = relevant_variables(sys)
        sys_indices = Int64[
            findfirst([v == var for var in all_relevant_variables])
            for v in sys_variables
            ]
        sys_dimensions = [ncategories(v) for v in sys_variables]

        other_variables = [v for v in all_relevant_variables if !(v in sys_variables)]
        other_indices = Int64[
            findfirst([v == var for var in all_relevant_variables])
            for v in other_variables
            ]
        other_dimensions = [ncategories(v) for v in other_variables]
        other_distribution = identity_distribution(D, prod(other_dimensions))

        factor_distribution = multiply_kron(sys_distribution, other_distribution)
        factor_indices = vcat(sys_indices, other_indices)
        factor_dimensions = vcat(sys_dimensions, other_dimensions)

        right_order = invperm(factor_indices)
        ordered_distribution = permute_distribution(factor_distribution, factor_dimensions, right_order)
        result_distribution = multiply_star(ordered_distribution, result_distribution)
    end
    DiscreteSystem{D}(all_parents, all_variables, result_distribution)
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