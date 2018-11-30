#=
distribution_operations:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-21
=#

import AcausalNets.Algebra: eye
import AcausalNets.Algebra: star, unstar, event
using LinearAlgebra

@testset "star" begin
    "return correct results"
    a = eye(5)
    b = eye(5)
    correct_diff = zeros(5, 5)
    result = star(a, b)
    @test (result - eye(5)) ≈ correct_diff

    "return correct results","return an AbstractMatrix"
    a = Diagonal([0.1, 0.8, 0.1])
    b = Diagonal([0.33, 0.66, 0.0])
    correct_result = Diagonal([0.033, 0.528, 0.0])
    result = star(a, b)
    @test result == correct_result
    @test typeof(result) <: AbstractMatrix

    "throw error if applied to matrices of unequal size"
    smaller = ones(6, 6)
    larger = ones(2, 2)
    @test_throws DimensionMismatch star(smaller, larger)
end

@testset "unstar" begin
    "revert the star operator","return an AbstractMatrix"
    a = Diagonal([0.3, 0.4, 0.5, 0.6, 0.7])
    b = Diagonal([0.1, 0.2, 0.3, 0.4, 0.5])
    c = star(a, b)
    d = unstar(c, b)
    @test d ≈ a
    @test typeof(d) <: AbstractMatrix
end

@testset "event" begin
    "throw error if applied to matrices of unequal size"
    smaller = ones(6, 6)
    larger = ones(2, 2)
    @test_throws DimensionMismatch event(smaller, larger)

    "correctly multiply evidence matrix into system matrix","return an AbstractMatrix"
    a = Diagonal([0.3, 0.3, 0.3])
    b = Diagonal([0.4, 0.3, 0.3])
    correct_result = Diagonal([0.16, 0.09, 0.09])
    result = event(a, b)
    @test result ≈ correct_result
    @test typeof(result) <: AbstractMatrix
end