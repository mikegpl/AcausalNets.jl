#=
visualization:
- Julia version: 
- Author: marcin
- Date: 2018-07-28
=#

using GraphPlot

import AcausalNets.Structures: DiscreteBayesNet

import AcausalNets.Inference: JoinTree


function Base.show(dbn::DiscreteBayesNet)
    node_names = [string(sys) for sys in dbn.systems]
    return gplot(dbn.dag, nodelabel = node_names)
end

function Base.show(jt::JoinTree)
    node_names = ["$i:" * string(jt.vertex_to_cluster[i]) for i in 1:length(jt.vertex_to_cluster)]
    return gplot(jt.graph, nodelabel = node_names)
end

