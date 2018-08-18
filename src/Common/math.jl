#=
math:
- Julia version: 
- Author: marcin
- Date: 2018-08-11
=#
using LinearAlgebra


star(A::AbstractMatrix, B::AbstractMatrix) = sqrt(B) * A * sqrt(B)

unstar(C, B) = star(C, pinv(B)) # odwrotność star (za #6)