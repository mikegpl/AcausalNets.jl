module Systems

    include("discrete_system.jl")
    include("discrete_quantum_system.jl")

    export
        relevant_variables,
        DiscreteQuantumSystem

end #module
