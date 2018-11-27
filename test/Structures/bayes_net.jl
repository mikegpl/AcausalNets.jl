#=
bayes_net:
- Julia version: 1.0
- Author: mikegpl
- Date: 2018-11-26
=#
import AcausalNets.Structures: AcausalNet
import AcausalNets.Systems: DiscreteQuantumSystem
import AcausalNets.Common: Variable

import AcausalNets.Structures:
    variables_names,
    variables,
    check_parents,
    check_variables,
    variable_to_node,
    system_to_node,
    parent_systems,
    family


function get_me_net()
    roA = diagm(0 => [
            .5, #A0
            .5  #A1
            ])

    roBwA = diagm(0 => [
            .6, #A0, B0
            .4, #A0, B1
            .5, #A1, B0
            .5, #A1, B1
            ])

    roCwA = diagm(0 => [
            .8, #A0, C0
            .2, #A0, C1
            .3, #A1, C0
            .7  #A1, C1
            ])

    roDwB = diagm(0 => [
            .5, #B0, D0
            .5, #B0, D1
            .1, #B1, D0
            .9  #B1, D1
            ])

    roEwC = diagm(0 => [
            .4, #C0, E0
            .6, #C0, E1
            .7, #C1, E0
            .3  #C1, E1
            ])

    roFwDE = diagm(0 => [
            .01, #D0, E0, F0
            .99, #D0, E0, F1
            .99, #D0, E1, F0
            .01, #D0, E1, F1
            .99, #D1, E0, F0
            .01, #D1, E0, F1
            .99, #D1, E1, F0
            .01  #D1, E1, F1
            ])

    roGwC = diagm(0 => [
            .9, #C0, G0
            .1, #C0, G1
            .2, #C1, G0
            .8  #C1, G1
            ])
    roHwEG = diagm(0 => [
            .05, #E0, G0, H0
            .95, #E0, G0, H1
            .05, #E0, G1, H0
            .95, #E0, G1, H1
            .05, #E1, G0, H0
            .95, #E1, G0, H1
            .95, #E1, G1, H0
            .05  #E1, G1, H1
            ])

    var_a = Variable(:a, 2)
    var_b = Variable(:b, 2)
    var_c = Variable(:c, 2)
    var_d = Variable(:d, 2)
    var_e = Variable(:e, 2)
    var_f = Variable(:f, 2)
    var_g = Variable(:g, 2)
    var_h = Variable(:h, 2)

    sys_a = DiscreteQuantumSystem([var_a], roA)
    sys_b = DiscreteQuantumSystem([var_a], [var_b], roBwA)
    sys_c = DiscreteQuantumSystem([var_a], [var_c], roCwA)
    sys_d = DiscreteQuantumSystem([var_b], [var_d], roDwB)
    sys_e = DiscreteQuantumSystem([var_c], [var_e], roEwC)
    sys_f = DiscreteQuantumSystem([var_d, var_e], [var_f], roFwDE)
    sys_g = DiscreteQuantumSystem([var_c], [var_g], roGwC)
    sys_h = DiscreteQuantumSystem([var_e, var_g], [var_h], roHwEG)


    example_an = AcausalNet()
    push!(example_an, sys_a)
    push!(example_an, sys_b)
    push!(example_an, sys_c)
    push!(example_an, sys_d)
    push!(example_an, sys_e)
    push!(example_an, sys_f)
    push!(example_an, sys_g)
    push!(example_an, sys_h)
    variables = [var_a, var_b, var_c, var_d, var_e, var_f, var_g, var_h]
    systems = [sys_a, sys_b, sys_c, sys_d, sys_e, sys_f, sys_g, sys_h]
    return variables, systems, example_an
end

get_me_system() = DiscreteQuantumSystem([Variable(:xd, 2)], [Variable(:pdk, 2)], diagm(0 => [.28, .11, .19,.96]))

@testset "variable_to_node" begin
    "return index of networks's DAG node which represents the variable", "return Int64"
    vars, systems, net = get_me_net()
    test_var = vars[1]
    result = variable_to_node(test_var, net)
    @test 1 == result
    @test Int64 == typeof(result)

    "return Nothing if requested node doesn't exist'"
    test_var = Variable(:xd, 2)
    result = variable_to_node(test_var, net)
    @test Nothing == typeof(result)
end


@testset "system_to_node" begin
    "return index of system in network's DAG", "return Int64"
    _, sys, net = get_me_net()
    test_system = sys[1]
    result = system_to_node(test_system, net)
    @test 1 == result
    @test Int64 == typeof(result)

    "return Nothing if requested node doesn't exist"
    test_system = get_me_system()
    result = system_to_node(test_system, net)
    @test Nothing == typeof(result)
end

@testset "variables" begin
    "return correct results", "return vector of variables"
    vars, _, net = get_me_net()
    result = variables(net)
    @test vars == result
    @test typeof(result) <: Vector{Variable}
end

@testset "variables_names" begin
    _, _, net = get_me_net()
    "return correct results", "return array of symbols"
    names = Symbol[:a, :b, :c, :d, :e, :f, :g, :h]
    @test names == variables_names(net)
    @test typeof(names) <: Array{Symbol}
end

@testset "check_parents" begin
    _, systems, net = get_me_net()
    "return true if system variables have all of parent nodes in network"
    test_system = systems[1]
    @test check_parents(test_system, net)

    "return false if system has no parent nodes in network"
    test_system = get_me_system()
    @test false == check_parents(test_system, net)
end

@testset "check_variables" begin
    _, systems, net = get_me_net()
    "return false if any of system variables already exists in network"
    test_system = systems[1]
    @test false == check_variables(test_system, net)

    "return true if none of system variables exist in network"
    test_system = get_me_system()
    @test check_variables(test_system, net)
end

@testset "family" begin
    "return set only with system itself if it doesn't have any parents"
    _, systems, net = get_me_net()
        test_system = systems[1]
    result = family(test_system, net)
    @test 1 == length(result)
    @test in(test_system, result)

    "return set with its parents if it has any"
    test_system = systems[2]
    test_parent = systems[1]
    result = family(test_system, net)
    @test in(test_parent, result)
    @test in(test_system, result)
    @test 2 == length(result)
end

@testset "parent_systems" begin
    "return empty set if system's node has no parents"
    _, systems, net = get_me_net()
    test_system = systems[1]
    result = parent_systems(test_system, net)
    @test 0 == size(result)[1]

    "return correct set of parents if system's node has any parents'"
    test_system = systems[2]
    test_parent = systems[1]
    result = parent_systems(test_system, net)
    # Set(Array -> problem
    @test in(test_parent, result)
    @test 1 == size(result)[1]
end

