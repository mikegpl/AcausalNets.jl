using QI
using LinearAlgebra
import AcausalNets.Common: Variable
import AcausalNets.Algebra: eye, star, unstar

const QuantumDistribution = Matrix{Complex{Float64}} # a bit quicker, allocates more memory
# const QuantumDistribution = Hermitian # a bit slower (due to copying matrices), less memory
# TODO to be changed to Hermitian when https://github.com/Jutho/TensorOperations.jl/issues/57 has been resolved


const QD = QuantumDistribution
const DiscreteQuantumSystem = DiscreteSystem{QuantumDistribution}

multiply_star(d1::QD, d2::QD)::QD = QD(star(d1, d2))

divide_star(d1::QD, d2::QD)::QD = QD(unstar(d1, d2))

multiply_kron(d1::QD, d2::QD)::QD = QD(kron(d1, d2))

sum_distribution(d::QD) = tr(d)

event(system::QD, e::QD)::QD = QD((e * system * e) / tr(e * system))

function check_distribution(distribution::QD, parents::Vector{Variable}, variables::Vector{Variable})::Bool
    dimensions = size(distribution)
    total_ncategories = prod([v.ncategories for v in vcat(parents, variables)])
    dimensions[1] == dimensions[2] == total_ncategories
end

permute_distribution(d::QD, dimensions::Vector{Int64}, order::Vector{Int64} )::QD = QD(permutesystems(d, dimensions, order))

function reduce_distribution(d::QD, dimensions::Vector{Int64}, reduce_ind::Vector{Int64})::QD
    if length(reduce_ind) > 0
        QD(ptrace(d, dimensions, reduce_ind))
    else
        d
    end
end

identity_distribution(::Type{QD}, size::Int64)::QD = QD(eye(size))