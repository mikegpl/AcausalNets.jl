#=
math:
- Julia version: 
- Author: marcin
- Date: 2018-08-11
=#
using LinearAlgebra

eye(size::Int64) = Diagonal([1 for i in 1:size])