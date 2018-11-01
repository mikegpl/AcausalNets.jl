#=
visualization:
- Julia version: 0.7
- Author: marcin
- Date: 2018-07-28
=#

using GraphPlot

import AcausalNets.Structures: DiscreteBayesNet

import AcausalNets.Inference: JoinTree, Inferrer


function Base.show(dbn::DiscreteBayesNet, verbose::Bool=false)
    node_names = [string(sys, verbose) for sys in dbn.systems]
    return gplot(dbn.dag, nodelabel = node_names)
end

function Base.show(jt::JoinTree, verbose::Bool=true)
    node_names = ["$i:" * string(jt.vertex_to_cluster[i], verbose) for i in 1:length(jt.vertex_to_cluster)]
    return gplot(jt.graph, nodelabel = node_names)
end

Base.show(inferrer::Inferrer) = show(inferrer.bayes_net)