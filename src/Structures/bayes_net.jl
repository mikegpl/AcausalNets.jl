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
variable_to_node(v::Variable, dbn::DiscreteBayesNet) = findfirst((sys -> v in variables(sys)), dbn.systems)
system_to_node(s::S, dbn::DiscreteBayesNet{S}) where S <: DiscreteSystem = findfirst([s==sys for sys in systems(dbn)])

systems(dbn::DiscreteBayesNet) = dbn.systems
variables(dbn::DiscreteBayesNet) = Vector{Variable}(vcat([variables(s) for s in systems(dbn)]...))
variables_names(dbn::DiscreteBayesNet) = [v.name for v in variables(dbn)]


check_parents(system, dbn) = all([p in variables(dbn) for p in parents(system)])
check_variables(system, dbn) = !any([v in variables(dbn) for v in variables(system)])

family(ds::S, dbn::DiscreteBayesNet{S}) where {S} = union(Set([ds]), Set(parent_systems(ds, dbn)))


function parent_systems(ds::S, dbn::DiscreteBayesNet) where S
    result = S[]
    for p in parents(ds)
        parent_system = systems(dbn)[variable_to_node(p, dbn)]
        if !(parent_system in result)
            push!(result, parent_system)
        end
    end
    return result
end

function Base.push!(dbn::DiscreteBayesNet{S}, system::S) where S <: DiscreteSystem
    check_parents(system, dbn) || error("Parents of the system missing from the structure!")
    check_variables(system, dbn) || error("Variables of the system already present in the structure!")

    expanded_system = expand_parents(system, dbn.systems)
    !add_vertex!(dbn.dag)
    push!(dbn.systems, expanded_system)
    for p in expanded_system.parents
        add_edge!(dbn.dag, variable_to_node(p, dbn), length(dbn.systems))
    end

    return dbn

end

