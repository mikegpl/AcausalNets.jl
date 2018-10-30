#=
distribution_operations:
- Julia version: 0.7
- Author: marcin
- Date: 2018-08-18
=#

using LinearAlgebra


star(A::AbstractMatrix, B::AbstractMatrix)::AbstractMatrix = sqrt(B) * A * sqrt(B)
# 'star' product as defined in Quantum inferring acausal structures and the Monty Hall problem,
# equation (2)

unstar(C, B)::AbstractMatrix = star(C, pinv(B)) # odwrotność star (za #6)

event(system::AbstractMatrix, e::AbstractMatrix) = (e * system * e) / tr(e * system)