#!/usr/bin/env -S julia -t16

using Gurobi

include("common.jl")
# TODO: eval.jl should be a "conditional dependency". I think this can be
# accomplished by moving it into some condition branch at the bottom of this
# file.
include("eval.jl")
#include("ac-uc-solver-000.jl")

import GOC3Benchmark as benchmark

"""A prototype of the solver we plan to run in Event 4.

Simplifications:

- If the network is > 6500 buses and we are running D1, we:
    - Skip the second OPF solve with rounded shunts
    - Schedule by only penalizing deviation from initial status

"""

function code1(case_file::String, time_limit::Int, division::Int, model::String, switching_allowed::Int, output_file_path::String; use_hsl::Bool=true)
    time_start = time()

    verbose = true

    # The character after the first "D" in the file name should be the
    # division number
    division_loc = findfirst("D", case_file)[1] + 1
    division = case_file[division_loc:division_loc]
    division = parse(Int, division)

    # The five characters after the first "N" in the file name should be
    # the number of buses
    network_loc = findfirst("N", case_file)[1] + 1
    network = case_file[network_loc:network_loc+4]
    nbus = parse(Int, network)

    println("Running solver for division $division on network with $nbus buses")

    if time_limit > 14400
        # If time_limit is greater than the largest time limit we expect, relax
        # the MIP time limit on D1 and D2
        if nbus > 20000
            division_to_mip_time_limit = Dict(1 => 3600, 2 => 14400, 3 => 14400)
        else
            division_to_mip_time_limit = Dict(1 => 3600, 2 => 7200, 3 => 7200)
        end
        mip_time_limit = division_to_mip_time_limit[division]
    else
        mip_time_limit = Int(round(time_limit / 2))
    end

    if nbus > 6500 && division == 1
        resolve_rounded_shunts = false
        schedule_close_to_initial = true
    else
        resolve_rounded_shunts = true
        schedule_close_to_initial = false
    end

    relax_copperplate_balances = (division == 2 && nbus > 20000)

    args = Dict(
        "case" => case_file,
        "solution_file" => output_file_path,
        "print_io_info" => verbose,
        "print_solution_info" => verbose,
        "print_program_info" => verbose,
        "print_solver_info" => verbose,
        "scheduler_time_limit" => mip_time_limit,
        "mip_optimizer" => Gurobi.Optimizer,
        "linear_solver" => "ma27",
        "allow_switching" => Bool(switching_allowed),
        "warmstart_mip" => true,

        # This runs OPF solves at individual time points in parallel, rounds
        # the result to respect ramp limits, and writes a solution file.
        # Note that the sequential solve is performed after this and a new
        # solution file is written if time allows.
        "sequential_opf" => false,
        "parallel_opf" => true,

        # Apply simplifications for large networks
        "schedule_close_to_initial" => schedule_close_to_initial,
        "resolve_rounded_shunts" => resolve_rounded_shunts,

        "relax_copperplate_balances" => relax_copperplate_balances,

        # If we are only penalizing deviation-from-initial status, including
        # reserves doesn't change the solution.
        #
        # We may want to restrict when reserves are included if this takes too
        # long.
        "include_reserves_in_schedule" => !schedule_close_to_initial,
    )

    println("Running AC-UC solver with arguments:")
    display(args)

    benchmark.run_ac_uc_solver(args)

    solve_time = time() - time_start

    # TODO: evalution_summary should be run when running this script directly,
    # but not when running via MyJulia1.jl
    #evaluation_summary(case_file, output_file_path, 1, solve_time)
end

# TODO: Check whether this file is the one being executed, rather than whether
# we are running interactively
if isinteractive() == false
    # NOTE: This script depends on arg parsing functionality in common.jl
    main(parse_goc_c3_args())
else
    # TODO: There is no settings.jl. This branch should not exist. This code is
    # only intended to be used as part of command line-executable scripts, e.g.
    # by running this file directoy or running it through MyJulia1
    include("settings.jl")
    main(settings)
end
