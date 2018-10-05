#=
visualization:
- Julia version: 
- Author: marcin
- Date: 2018-07-28
=#

using GraphPlot

import AcausalNets.Structures: DiscreteBayesNet

function Base.show(dbn::DiscreteBayesNet)
    node_names = [string(sys) for sys in dbn.systems]
    return gplot(dbn.dag, nodelabel = node_names)
end

