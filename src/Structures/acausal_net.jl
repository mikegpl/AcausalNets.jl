#=
acausal_net:
- Julia version: 1.0
- Author: marcin
- Date: 2018-07-28
=#

using AcausalNets.Systems

const AcausalNet = DiscreteBayesNet{DiscreteQuantumSystem}
AcausalNet() = DiscreteBayesNet{DiscreteQuantumSystem}(DiscreteQuantumSystem)
