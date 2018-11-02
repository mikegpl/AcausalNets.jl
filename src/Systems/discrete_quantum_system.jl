using QI
using LinearAlgebra
import AcausalNets.Common: Variable
import AcausalNets.Algebra: eye, star, unstar

const QuantumDistribution = AbstractMatrix
const DiscreteQuantumSystem = DiscreteSystem{QuantumDistribution}

multiply_star(d1::QuantumDistribution, d2::QuantumDistribution) = star(d1, d2)
divide_star(d1::QuantumDistribution, d2::QuantumDistribution) = unstar(d1, d2)
multiply_kron(d1::QuantumDistribution, d2::QuantumDistribution) = kron(d1, d2)
sum_distribution(d::QuantumDistribution) = tr(d)

function check_distribution(
    distribution::QuantumDistribution,
    parents::Vector{Variable},
    variables::Vector{Variable}
    )
    dimensions = size(distribution)
    total_ncategories = prod([v.ncategories for v in vcat(parents, variables)])
    dimensions[1] == dimensions[2] == total_ncategories
end

permute_distribution(d::QuantumDistribution, dimensions::Vector{Int64}, order::Vector{Int64}) = permute_systems(d, dimensions, order)

function reduce_distribution(d::QuantumDistribution, dimensions::Vector{Int64}, reduce_ind::Vector{Int64})
    if length(reduce_ind) > 0
        ptrace(d, dimensions, reduce_ind)
    else
        d
    end
end



function identity_distribution(::Type{D}, size::Int64)::D where {D <: QuantumDistribution}
    eye(size)
end