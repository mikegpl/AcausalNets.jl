using QI
import AcausalNets.Common: Variable
const DiscreteQuantumSystem{D <: Matrix} = DiscreteSystem{D}


function DiscreteQuantumSystem(
    parents::Vector{Variable},
    variables::Vector{Variable},
    distribution::Matrix,
)
    dimensions = size(distribution)
    total_ncategories = prod([v.ncategories for v in vcat(parents, variables)])
    dimensions[1] == dimensions[2] == total_ncategories || error("Dimensions ($(dimensions[1]), $(dimensions[2]), $total_ncategories) not matching!")

    DiscreteQuantumSystem{Matrix}(parents, variables, distribution)
end

DiscreteQuantumSystem(variables::Vector{Variable}, distribution::Matrix) = DiscreteQuantumSystem(Variable[], variables, distribution)

function enforce_parents_order(dqs::DiscreteQuantumSystem, existing_variables::Vector{Variable})
    new_parents_order = [v for v in existing_variables if v in dqs.parents]
    new_variables_order = dqs.variables
    new_vars_order = vcat(new_parents_order, new_variables_order)
    new_indexing = [findfirst(new_vars_order, p) for p in relevant_variables(dqs)]
    dimensions = [v.ncategories for v in relevant_variables(dqs)]
    new_distribution = permute_systems(dqs.distribution, dimensions, new_indexing)
    DiscreteQuantumSystem(new_parents_order, new_variables_order, new_distribution)
end