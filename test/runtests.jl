using AcausalNets
using Test
``
all_tests = vcat("Algebra", "Common", "Inference", "Structures", "Systems")

for test_suite in all_tests
    println(test_suite)
end
