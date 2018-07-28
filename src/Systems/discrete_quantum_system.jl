using AcausalNets.Common: Variable
const DiscreteQuantumSystem{D <: Matrix} = DiscreteSystem{D}


function DiscreteQuantumSystem(
    variables::Vector{Variable},
    distribution::Matrix,
    parents::Vector{Variable} = Variable[]
)
    dimensions = size(distribution)
    total_ncategories = prod([v.ncategories for v in vcat(parents, variables)])
    if dimensions[1] == dimensions[2] == total_ncategories
        return DiscreteQuantumSystem{Matrix}(parents, variables, distribution)
     else
        throw(ArgumentError("Dimensions ($(dimensions[1]), $(dimensions[2]), $total_ncategories) not matching!"))
    end
end
