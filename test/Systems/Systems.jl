#=
Systems:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-17
=#

module Systems
    for file in ["discrete_quantum_system.jl", "discrete_system.jl"]
        include(file)
    end
end