#=
quantum_bayesian:
- Author: marcin
- Date: 2018-06-02
=#

using LightGraphs
import AcausalNets.Systems:
    DiscreteSystem,
    expand_parents,
    enforce_parents_order,
    parents,
    variables,
    variables_names
import AcausalNets.Common: VariableName, Variable
import Compat.Iterators: flatten

const DAG = DiGraph

struct DiscreteBayesNet{T <: DiscreteSystem}
    dag             ::DAG
    systems         ::Vector{T}

    DiscreteBayesNet{S}(::Type{S}) where {S <: DiscreteSystem} = new(DAG(0), S[])
end


# Returns the index of dbn's DAG node which represents the variable
variable_to_node(dbn::DiscreteBayesNet, v::Variable) = findfirst((sys -> v in variables(sys)), dbn.systems)
system_to_node(dbn::DiscreteBayesNet{S}, s::S) where S <: DiscreteSystem = findfirst(systems(dbn), s)

systems(dbn::DiscreteBayesNet) = dbn.systems
variables(dbn::DiscreteBayesNet) = Vector{Variable}(vcat([variables(s) for s in systems(dbn)]...))
variables_names(dbn::DiscreteBayesNet) = [v.name for v in variables(dbn)]


check_parents(dbn, system) = all([p in variables(dbn) for p in parents(system)])
check_variables(dbn, system) = !any([v in variables(dbn) for v in variables(system)])



function Base.push!(dbn::DiscreteBayesNet{S}, system::S) where S <: DiscreteSystem
    check_parents(dbn, system) || error("Parents of the system missing from the structure!")
    check_variables(dbn, system) || error("Variables of the system already present in the structure!")

    expanded_system = expand_parents(system, dbn.systems)
    perm_system = enforce_parents_order(expanded_system, variables(dbn))

    !add_vertex!(dbn.dag)
    push!(dbn.systems, perm_system)
    for p in perm_system.parents
        add_edge!(dbn.dag, variable_to_node(dbn, p), length(dbn.systems))
    end

    return dbn

end

