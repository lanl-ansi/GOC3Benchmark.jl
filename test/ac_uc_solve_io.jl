using Test
import JSON
import GOC3Benchmark as goc3

filedir = dirname(@__FILE__)
datadir = joinpath(filedir, "data")

"""
Check the solution file for correct structure and that some
bounds are satisfied.
"""
function check_solution_file(problem_file, solution_file)
    problem_data = Dict{String, Any}()
    open(problem_file) do io
        problem_data = JSON.parse(io)
    end
    input_data = goc3.process_input_data(problem_data)
    solution_data = Dict{String, Any}()
    open(solution_file) do io
        solution_data = JSON.parse(io)
    end
    @test "time_series_output" in keys(solution_data)
    ts_out = solution_data["time_series_output"]
    @test "bus" in keys(ts_out)
    @test "shunt" in keys(ts_out)
    @test "simple_dispatchable_device" in keys(ts_out)
    @test "ac_line" in keys(ts_out)
    @test "two_winding_transformer" in keys(ts_out)
    @test "dc_line" in keys(ts_out)

    @test length(ts_out["bus"]) == length(input_data.bus_ids)
    @test length(ts_out["shunt"]) == length(input_data.shunt_ids)
    @test length(ts_out["simple_dispatchable_device"]) == length(input_data.sdd_ids)
    @test length(ts_out["ac_line"]) == length(input_data.ac_line_ids)
    @test length(ts_out["two_winding_transformer"]) == length(input_data.twt_ids)
    @test length(ts_out["dc_line"]) == length(input_data.dc_line_ids)

    n_periods = length(input_data.periods)

    # Make sure bus voltage magnitudes are within bounds
    for bus in ts_out["bus"]
        @test "uid" in keys(bus)
        @test "vm" in keys(bus)
        uid = bus["uid"]
        ub = input_data.bus_lookup[uid]["vm_ub"]
        lb = input_data.bus_lookup[uid]["vm_lb"]
        @test length(bus["vm"]) == n_periods
        for vm in bus["vm"]
            @test lb .- 1e-8 <= vm && vm <= ub .+ 1e-8
        end
    end

    # Make sure generator levels are within bounds
    for sdd in ts_out["simple_dispatchable_device"]
        @test "uid" in keys(sdd)
        @test "on_status" in keys(sdd)
        @test "p_on" in keys(sdd)
        @test "q" in keys(sdd)
        uid = sdd["uid"]
        p_lb_array = input_data.sdd_ts_lookup[uid]["p_lb"]
        p_ub_array = input_data.sdd_ts_lookup[uid]["p_ub"]
        q_lb_array = input_data.sdd_ts_lookup[uid]["q_lb"]
        q_ub_array = input_data.sdd_ts_lookup[uid]["q_ub"]

        on_status_array = sdd["on_status"]
        p_on_array = sdd["p_on"]
        q_array = sdd["q"]
        @test length(sdd["on_status"]) == n_periods
        @test length(sdd["p_on"]) == n_periods
        @test length(sdd["q"]) == n_periods
        # Note that these test "<=" for all coordinates
        @test on_status_array .* p_lb_array .- 1e-8 <= p_on_array
        @test p_on_array <= on_status_array .* p_ub_array .+ 1e-8
        @test q_lb_array .- 1e-8 <= q_array && q_array <= q_ub_array .+ 1e-8
        for u_on in on_status_array
            @test u_on == 0 || u_on == 1
        end
    end

end

@testset "Test run_ac_uc_solver function" begin

    # Division 1
    problem_file = joinpath(datadir, "C3E4N00073D1_scenario_303.json")
    solution_file = joinpath(tempdir(), "C3E4N00073D1_scenario_303_solution.json")
    args = Dict(
        "case" => problem_file,
        "solution_file" => solution_file,
        "parallel_opf" => false
    )
    # This runs the solver and puts the solution file in the specified directory
    goc3.run_ac_uc_solver(args)
    check_solution_file(problem_file, solution_file)

    # Division 2
    problem_file = joinpath(datadir, "C3E4N00073D2_scenario_303.json")
    solution_file = joinpath(tempdir(), "C3E4N00073D2_scenario_303_solution.json")
    args = Dict(
        "case" => problem_file,
        "solution_file" => solution_file,
        "parallel_opf" => false
    )
    # This runs the solver and puts the solution file in the specified directory
    goc3.run_ac_uc_solver(args)
    check_solution_file(problem_file, solution_file)

    # Division 3
    problem_file = joinpath(datadir, "C3E4N00073D3_scenario_303.json")
    solution_file = joinpath(tempdir(), "C3E4N00073D3_scenario_303_solution.json")
    args = Dict(
        "case" => problem_file,
        "solution_file" => solution_file,
        "parallel_opf" => false
    )
    # This runs the solver and puts the solution file in the specified directory
    goc3.run_ac_uc_solver(args)
    check_solution_file(problem_file, solution_file)

end
