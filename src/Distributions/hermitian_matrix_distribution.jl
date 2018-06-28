#=
hermitian_matrix_distribution:
- Author: mprzewie, mikegpl
- Date: 2018-06-28
=#

struct HermitianMatrix{T <: Union{Complex{Float64}, Float64}}
    matrix_size::Int
    matrix::Matrix{T}

    # Replaces default constructor. Must be called with T specialization:
    # Correct: HermitianMatrix{Complex{Int64}}(size, matrix)
    # Incorrect: HermitianMatrix(size, matrix)
    function HermitianMatrix{T}(matrix_size::Int, matrix::Matrix{T}) where T
        x, y = size(matrix)
        @check_args(HermitianMatrix, matrix_size == x == y)
        new(matrix_size, matrix)
    end
end

event(system::Matrix, e::Matrix) = (e * system * e) / trace(e * system)
event(system::HermitianMatrix, e::Matrix) = HermitianMatrix(event(system.matrix, e))

HermitianMatrix{T <: Union{Complex{Float64}, Float64}}(p::Matrix{T}) = HermitianMatrix{T}(length(diag(p)),p)
HermitianMatrix{T <: Union{Complex{Float64}, Float64}}(v::Vector{T}) = HermitianMatrix{T}(length(v), diagm(v))
# TODO - construtor from integer (K => diagonal of 1/k)

Base.convert(::Type{HermitianMatrix}, m::Matrix) = HermitianMatrix(m)
Base.convert(::Type{Matrix}, hm::HermitianMatrix) = hm.matrix
