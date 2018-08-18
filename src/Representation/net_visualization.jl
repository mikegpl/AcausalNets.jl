#=
visualization:
- Julia version: 
- Author: marcin
- Date: 2018-07-28
=#

# using TikzGraphs
#
# import AcausalNets.Structures: DiscreteBayesNet
#
# function Base.show(dbn::DiscreteBayesNet)
#
#     node_names = [string(sys) for sys in dbn.systems]
#     return TikzGraphs.plot(dbn.dag) #, node_names)
# end