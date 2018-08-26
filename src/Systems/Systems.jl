module Systems

    include("discrete_system.jl")
    include("discrete_quantum_system.jl")

    export
        parents,
        variables,
        distribution,
        parents_names,
        variables_names,
        relevant_variables,
        is_parent,
        DiscreteQuantumSystem

end #module
