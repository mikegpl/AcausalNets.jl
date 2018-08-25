module Systems

    include("discrete_system.jl")
    include("discrete_quantum_system.jl")

    export
        parents,
        variables,
        parents_names,
        variables_names,
        relevant_variables,
        is_parent,
        DiscreteQuantumSystem

end #module
