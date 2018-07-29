using QI
import AcausalNets.Common: Variable

const QuantumDistribution = Matrix
const DiscreteQuantumSystem = DiscreteSystem{QuantumDistribution}

function check_distribution(
    distribution::QuantumDistribution,
    parents::Vector{Variable},
    variables::Vector{Variable}
    )
    dimensions = size(distribution)
    total_ncategories = prod([v.ncategories for v in vcat(parents, variables)])
    dimensions[1] == dimensions[2] == total_ncategories
end

function prepend_parent!(dqs::DiscreteQuantumSystem, var::Variable)
    if !(var in parents(dqs))
        new_parents = vcat([var], parents(dqs))
        new_distribution = kron(eye(var.ncategories), dqs.distribution)
        return DiscreteQuantumSystem(new_parents, variables(dqs), new_distribution)
    else
        return dqs
    end
end

# not to confuse with QI's permute_systems - this is a higher-level implementation
function permute_system(dqs::DiscreteQuantumSystem, new_parent_indexing, new_variable_indexing)
    new_indexing = vcat(new_parent_indexing, new_variable_indexing + length(parents(dqs)))
    dimensions = [v.ncategories for v in relevant_variables(dqs)]
    new_distribution = permute_systems(dqs.distribution, dimensions, new_indexing)
    new_parents = [parents(dqs)[i] for i in new_parent_indexing]
    new_variables = [variables(dqs)[i] for i in new_variable_indexing]

    DiscreteQuantumSystem(new_parents, new_variables, new_distribution)
end