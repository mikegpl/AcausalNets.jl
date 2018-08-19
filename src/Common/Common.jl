module Common

    include("names.jl")
    include("variables.jl")
    include("sets.jl")

    export
        VariableName,
        VariableNames,
        Variable,
        ncategories

end
