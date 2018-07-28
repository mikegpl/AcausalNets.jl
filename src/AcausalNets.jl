module AcausalNets
    __precompile__()

    # package code goes here
    include("utils.jl")

    include(joinpath("Common", "Common.jl"))
    include(joinpath("Systems", "Systems.jl"))
    include(joinpath("Structures", "Structures.jl"))

    using Reexport
    @reexport using AcausalNets.Common
    @reexport using AcausalNets.Systems
    @reexport using AcausalNets.Structures

end # module
