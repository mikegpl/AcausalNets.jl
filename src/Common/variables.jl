struct Variable
    name        ::VariableName
    ncategories ::Int
end

Base.convert(::Type{VariableName}, var::Variable) = var.name
