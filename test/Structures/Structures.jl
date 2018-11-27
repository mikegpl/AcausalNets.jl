#=
Structures:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-16
=#

module Structures
    for file in ["acausal_net.jl", "bayes_net.jl"]
        include(file)
    end
end
