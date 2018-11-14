#=
visualization:
- Julia version: 1.0
- Author: marcin
- Date: 2018-07-28
=#

using GraphPlot
using LightGraphs
import AcausalNets.Structures: DiscreteBayesNet

import AcausalNets.Inference: JoinTree

function Base.show(dbn::DiscreteBayesNet, verbose::Bool=false)
    node_names = [string(sys, verbose) for sys in dbn.systems]
    return gplot(dbn.dag, nodelabel = node_names)
end

function Base.show(jt::JoinTree, verbose::Bool=true)
    node_names = ["$i:" * string(jt.vertex_to_cluster[i], verbose) for i in 1:length(jt.vertex_to_cluster)]
    edge_names = [string(jt.edge_to_sepset[Set([src(e), dst(e)])], verbose) for e in edges(jt.graph)]

    return gplot(jt.graph, nodelabel = node_names, edgelabel = edge_names)
end
