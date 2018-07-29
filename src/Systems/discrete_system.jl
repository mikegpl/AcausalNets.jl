using AcausalNets.Common

struct DiscreteSystem{D}
    parents         ::Vector{Variable}
    variables       ::Vector{Variable}
    distribution    ::D

end

parents(ds::DiscreteSystem) = ds.parents
variables(ds::DiscreteSystem) = ds.variables
parents_names(ds::DiscreteSystem) = [p.name for p in parents(ds)]
variables_names(ds::DiscreteSystem) = [v.name for v in variables(ds)]
Base.string(ds::DiscreteSystem) = join([String(n) for n in variables_names(ds)], ", ")

relevant_variables(system::DiscreteSystem) = vcat(system.parents, system.variables)

