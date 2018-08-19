#=
common:
- Julia version: 
- Author: marcin
- Date: 2018-08-18
=#

using LinearAlgebra

eye(size::Int64)::Diagonal = Diagonal(ones(size))