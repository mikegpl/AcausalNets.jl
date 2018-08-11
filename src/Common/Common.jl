module Common

    include("names.jl")
    include("variables.jl")
    include("math.jl")

    export
        VariableName,
        VariableNames,
        Variable,
        eye

end
