#=
Variableiables:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-21
=#
import AcausalNets.Common: ncategories, Variable

@testset "ncategories" begin
    "return correct result"
    variables = [Variable(:a, 1), Variable(:b, 3), Variable(:c, 10)]
    correct_result = 30
    @test correct_result == ncategories(variables)
    @test 5 == ncategories(Variable(:a, 5))
end