#=
Algebra:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-15
=#

module Algebra
    for file in ["common.jl", "distribution_operations"]
        include(file)
    end
end