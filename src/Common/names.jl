const VariableName = Symbol
const VariableNames = AbstractVector{VariableName}
Base.convert(::Type{VariableNames}, name::VariableName) = [name]
