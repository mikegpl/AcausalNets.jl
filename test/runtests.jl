using AcausalNets
using QuantumInformation
using Test

all_tests = vcat("Algebra", "Common", "Structures")
include_module(module_name::String) = include(joinpath(module_name, string(module_name, ".jl")))

for test_suite in all_tests
    include_module(test_suite)
end
