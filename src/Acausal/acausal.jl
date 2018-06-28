#=
acausal_structures:
- Author: marcin
- Date: 2018-05-10
=#

using BayesNets.CPDs.DiscreteCPD
using BayesNets.DAG
using BayesNets.DiscreteBayesNet

include("discrete_mcpd.jl")
include("acausal_structure.jl")


export
    AcausalStructure,
    DiscreteQCPD





