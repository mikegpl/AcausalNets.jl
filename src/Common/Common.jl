module Common

    include("names.jl")
    include("variables.jl")

    println("exporting common")
    export
        VariableName,
        VariableNames,
        Variable

end
