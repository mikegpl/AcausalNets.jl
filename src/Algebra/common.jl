#=
common:
- Julia version: 1.0
- Author: marcin
- Date: 2018-08-18
=#

using LinearAlgebra

eye(size::Int64)::Diagonal = Diagonal(ones(size))