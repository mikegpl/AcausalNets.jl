module AcausalNets
    __precompile__()

    include_module(module_name::String) = include(joinpath(module_name, string(module_name, ".jl")))
    # package code goes here
    include("utils.jl")

    include_module("Common")
    include_module("Algebra")
    include_module("Systems")
    include_module("Structures")
    include_module("Inference")
    include_module("Representation")
#
    using Reexport
    @reexport using AcausalNets.Algebra
    @reexport using AcausalNets.Common
    @reexport using AcausalNets.Systems
    @reexport using AcausalNets.Structures
    @reexport using AcausalNets.Inference
    @reexport using AcausalNets.Representation

end # module
