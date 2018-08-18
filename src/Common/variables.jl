struct Variable
    name        ::VariableName
    ncategories ::Int
end

Base.convert(::Type{VariableName}, var::Variable) = var.name

ncategories(variables::Vector{Variable}) = prod([v.ncategories for v in variables])
ncategories(variable::Variable) = ncategories([variable])