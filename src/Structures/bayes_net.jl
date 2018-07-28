#=
quantum_bayesian:
- Author: marcin
- Date: 2018-06-02
=#

import AcausalNets.Systems: DiscreteSystem, variables_names
import AcausalNets.Common: VariableName
import LightGraphs: DiGraph

const DAG = DiGraph

struct DiscreteBayesNet{T <: DiscreteSystem}
    dag     ::DAG
    systems ::Vector{T}
end

DiscreteBayesNet() = DiscreteBayesNet(DAG(0), DiscreteSystem[])
DiscreteBayesNet{T <: DiscreteSystem}(::Type{T}) = DiscreteBayesNet(DAG(0), T[])
