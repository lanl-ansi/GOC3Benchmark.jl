using Test
import JSON
import JuMP
import GOC3Benchmark as goc3

filedir = dirname(@__FILE__)
datadir = joinpath(filedir, "data")

@testset "Test build and solve models" begin

    # Division 1
    problem_file = joinpath(datadir, "C3E4N00073D1_scenario_303.json")
    problem_data = Dict{String, Any}()
    open(problem_file) do io
        problem_data = JSON.parse(io)
    end
    input_data = goc3.process_input_data(problem_data)
    model, schedule_data = goc3.schedule_power_copperplate(input_data)
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(model) == JuMP.OPTIMAL

    opf_model, opf_data = goc3.compute_optimal_power_flow_at_interval(
        input_data, schedule_data, 1
    )
    @test JuMP.primal_status(opf_model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(opf_model) == JuMP.LOCALLY_SOLVED

    # Division 2
    problem_file = joinpath(datadir, "C3E4N00073D2_scenario_303.json")
    problem_data = Dict{String, Any}()
    open(problem_file) do io
        problem_data = JSON.parse(io)
    end
    input_data = goc3.process_input_data(problem_data)
    model, schedule_data = goc3.schedule_power_copperplate(input_data)
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(model) == JuMP.OPTIMAL

    opf_model, opf_data = goc3.compute_optimal_power_flow_at_interval(
        input_data, schedule_data, 1
    )
    @test JuMP.primal_status(opf_model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(opf_model) == JuMP.LOCALLY_SOLVED

    ## Division 3
    problem_file = joinpath(datadir, "C3E4N00073D3_scenario_303.json")
    problem_data = Dict{String, Any}()
    open(problem_file) do io
        problem_data = JSON.parse(io)
    end
    input_data = goc3.process_input_data(problem_data)
    model, schedule_data = goc3.schedule_power_copperplate(input_data)
    @test JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(model) == JuMP.OPTIMAL

    opf_model, opf_data = goc3.compute_optimal_power_flow_at_interval(
        input_data, schedule_data, 1
    )
    @test JuMP.primal_status(opf_model) == JuMP.FEASIBLE_POINT
    @test JuMP.termination_status(opf_model) == JuMP.LOCALLY_SOLVED

end
