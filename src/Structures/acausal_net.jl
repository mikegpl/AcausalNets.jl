#=
acausal_net:
- Julia version: 
- Author: marcin
- Date: 2018-07-28
=#

using AcausalNets.Systems

const AcausalNet = DiscreteBayesNet{DiscreteQuantumSystem}
AcausalNet() = DiscreteBayesNet(DiscreteQuantumSystem)
