using AcausalNets.Common

struct DiscreteSystem{D}
    parents         ::Vector{Variable}
    variables       ::Vector{Variable}
    distribution    ::D

end


variables_names(ds::DiscreteSystem) = [v.name for v in ds.variables]
Base.string(ds::DiscreteSystem) = join([String(n) for n in variables_names(ds)], ", ")

relevant_variables(system::DiscreteSystem) = system.parents + system.variables
