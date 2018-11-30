#=
sets:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-21
=#
import AcausalNets.Common: is_subset

@testset "is_subset" begin
    "return correct values"
    main = Set([1, 2, 3, 4])
    subset = Set([1, 2])
    not_subset = Set([4, 5, 6])
    @test is_subset(subset, main)
    @test !is_subset(not_subset, main)
end