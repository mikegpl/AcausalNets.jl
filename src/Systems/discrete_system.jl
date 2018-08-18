using AcausalNets.Common
import Base: ==

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