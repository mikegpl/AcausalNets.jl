#=
math:
- Julia version: 
- Author: marcin
- Date: 2018-08-11
=#
using LinearAlgebra

# eye(size::Int64) = Diagonal([1 for i in 1:size])

star(A::AbstractMatrix, B::AbstractMatrix) = sqrtm(B) * A * sqrtm(B)

unstar(C, B) = star(C, pinv(B)) # odwrotność star (za #6)