using AcausalNets
using Test
``
all_tests = vcat("Algebra")#, "Common", "Inference", "Structures", "Systems")
include_module(module_name::String) = include(joinpath(module_name, string(module_name, ".jl")))

for test_suite in all_tests
    include_module(test_suite)
    println(test_suite)
end
