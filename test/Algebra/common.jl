#=
common:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-16
=#
import AcausalNets.Algebra: eye
using LinearAlgebra

@testset "eye" begin
    matrix_size = 5
    matrix = eye(matrix_size)
    test_diagonal = [1 for _=1:5]
    result_diagonal = diag(matrix)

    @test result_diagonal == test_diagonal
    @test isdiag(matrix)
    @test matrix isa AbstractArray
end