#!/usr/bin/env julia

using ArgParse
using JSON

using JuMP
#import MathOptInterface as MOI
using HiGHS
using Ipopt

#include("scheduler.jl")
# : schedule_power_copperplate

#include("opf_model.jl")
# : get_ac_opf_model, extract_data_from_model

#include("opf.jl")
# This already includes opf_model.jl and process_bounds.jl
# : compute_optimal_power_flow_at_interval

#include("process_data.jl")
# : process_input_data

#include("process_bounds.jl")
# : tighten_bounds_using_ramp_limits

#include("reserves.jl")

#include("solution_data.jl")
# : construct_solution_dict


"""
AC-unit commitment solver

args must contain "case", which specifies the file to use for input data.

To change the linear solver used by Ipopt (default MA57), use e.g.:
    args["linear_solver"] = "mumps"
To fix real powers to scheduled values, use:
    args["fix_real_power"] = true
To require strict satisfaction (within tolerance) of power balances, use:
    args["relax_power_balance"] = false 
To penalize deviation of real powers from scheduled values, use:
    args["penalize_power_deviation"] = true

"""
function run_ac_uc_solver(args::Dict)
    # TODO: I think optional arguments would be less error-prone than a dict
    # of args here.
    start_time = time()

    timing_data = Dict{String, Float64}()

    #
    # Process arguments
    #
    if !("case" in keys(args))
        throw(Exception)
    end

    # scheduler time limit setting
    scheduler_time_limit = get(args, "scheduler_time_limit", -1.0)
    # Note that time_limit in MyJulia1 is an int. Here we convert to float
    # for consistency with JuMP.set_time_limit_sec
    scheduler_time_limit = Float64(scheduler_time_limit)

    print_io_info = get(args, "print_io_info", false)
    print_solution_info = get(args, "print_solution_info", false)
    print_program_info = get(args, "print_program_info", false)
    print_solver_info = get(args, "print_solver_info", false)
    print_projected_devices = get(args, "print_projected_devices", true)
    print_violated_qbalance = get(args, "print_violated_qbalance", false)

    simultaneous_acuc = get(args, "simultaneous_acuc", false)
    include_reserves_in_acuc = get(args, "include_reserves_in_acuc", false)
    minlp_optimizer = get(args, "minlp_optimizer", nothing)
    if simultaneous_acuc && minlp_optimizer === nothing
        throw(ArgumentError("minlp_optimizer must be provided if simultaneous_acuc is used"))
    end

    #
    # Options for MIP solve
    #
    mip_optimizer = get(args, "mip_optimizer", HiGHS.Optimizer)
    relax_integrality = get(args, "relax_integrality", false)
    halt_on_infeasible_schedule = get(args, "halt_on_infeasible_schedule", false)
    include_reserves_in_schedule = get(args, "include_reserves_in_schedule", true)
    relax_copperplate_balances = get(args, "relax_copperplate_balances", false)
    warmstart_mip = get(args, "warmstart_mip", false)
    overcommitment_factor = get(args, "overcommitment_factor", 1.05)

    # The default is to first attempt the copper-plate scheduling model, then
    # the relaxed copper-plate scheduling model, then the independent scheduling
    # model. This option overrides these to go straight to the independent
    # scheduling model.
    schedule_independently = get(args, "schedule_independently", false)

    # This option overrides the "copper-plate" and "independent" scheduling
    # problems and instead solves a problem where deviation from the initial
    # status is penalized.
    schedule_close_to_initial = get(args, "schedule_close_to_initial", false)
    if schedule_independently && schedule_close_to_initial
        println(
            """schedule_independently and schedule_to_initial_point override
            the defualt scheduling procedure and cannot both be set
            """
        )
        throw(Exception)
    end

    case_json_file = args["case"]
    if print_io_info
        println("loading: $(case_json_file)")
    end

    tighten_bounds_pre_opf = get(args, "tighten_bounds_pre_opf", true)

    #
    # Options for NLP solves
    #
    nlp_optimizer = get(args, "nlp_optimizer", nothing)
    linear_solver = get(args, "linear_solver", nothing)
    fix_ac_real_power = get(args, "fix_ac_real_power", false)
    sequential_opf = get(args, "sequential_opf", true)
    parallel_opf = get(args, "parallel_opf", false)
    multiperiod_opf = get(args, "multiperiod_opf", false)
    # This is only used if parallel_opf is true
    post_parallel_sequential = get(args, "post_parallel_sequential", true)
    if sequential_opf && parallel_opf
        println("\"sequential_opf\" and \"parallel_opf\" cannot both be true")
        throw(Exception)
    end
    allow_switching = get(args, "allow_switching", true)
    resolve_rounded_shunts = get(args, "resolve_rounded_shunts", true)
    relax_thermal_limits = get(args, "relax_thermal_limits", true)
    sdd_fraction_to_constrain = get(args, "sdd_fraction_to_constrain", 0.0)

    postprocess_final_solution = get(args, "postprocess_final_solution", true)

    write_duals = get(args, "write_duals", false)
    ###

    #
    # Read input data from file
    #
    case_data = Dict{String,Any}()

    open(case_json_file, "r") do io
        case_data = JSON.parse(io)
    end
    ###

    #
    # Process input data into convenient format for analyzing
    # solution from model
    #
    input_data = process_input_data(case_data)
    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
     azr_lookup, azr_ts_lookup, azr_ids, rzr_lookup, rzr_ts_lookup, rzr_ids,
    ) = input_data
    ###

    #
    # Initialize data structure for results of solve
    #
    metadata = Dict{String, Any}(
        "input_file" => case_json_file,
        "output_file" => get(args, "solution_file", nothing),
    )
    solve_data = Dict{String, Any}(
        "metadata" => metadata,
        "feasible" => nothing,
        "run_time" => nothing,
        "objective_value" => nothing,
    )
    ###

    # TODO
    # Before scheduling, attempt to tighten time-dependent bounds based
    # on limits implied by adjacent bounds and ramping constraints.
    # processed_data = tighten_bounds_using_ramp_limits(input_data)

    timing_data["load_initial_data"] = time() - start_time

    if simultaneous_acuc
        if fix_ac_real_power
            throw(ArgumentError("simultaneous_acuc cannot be used with fix_ac_real_power"))
        end
        if resolve_rounded_shunts
            throw(ArgumentError("simultaneous_acuc does not support re-solving shunts yet"))
        end
        if relax_thermal_limits
            throw(ArgumentError("simultaneous_acuc does not support relax_thermal_limits"))
        end
        acuc_model = get_ac_uc_model(
            input_data,
            include_reserves=include_reserves_in_acuc,
        )
        JuMP.set_optimizer(acuc_model, minlp_optimizer)
        JuMP.optimize!(acuc_model)

        # We can extract scheduling data easily as this model uses the same
        # variable names as the full model
        schedule_data = extract_data_from_scheduling_model(
            input_data,
            acuc_model,
            include_reserves=include_reserves_in_acuc,
        )
        # We should be able to extract OPF data by simply overriding the default
        # variable names used for p and q
        mpacopf_solution = extract_data_from_multiperiod_model(
            acuc_model,
            input_data,
            schedule_data.on_status,
            p_sdd=acuc_model[:p],
            q_sdd=acuc_model[:q],
        )
        # Reshape the solution data for consistency with other subroutines
        acopf_solutions = [
            (nothing, Dict(
                key => Dict(
                    uid => Dict(
                        attr => device[attr][i] for attr in keys(device)
                    ) for (uid, device) in devicedict
                ) for (key, devicedict) in mpacopf_solution
            )) for i in periods
        ]
    else
        time_mip_start = time()

        if warmstart_mip
            # Note that this simply copies initial status and does not solve
            # a MIP
            naive_schedule = schedule_to_initial_point(input_data)
        else
            naive_schedule = nothing
        end

        if schedule_independently
            model, schedule_data = schedule_power_independently(
                input_data;
                optimizer = mip_optimizer,
                relax_integrality = relax_integrality,
                time_limit = scheduler_time_limit,
                include_reserves = include_reserves_in_schedule,
                warmstart_data = naive_schedule,
            )
        elseif schedule_close_to_initial
            model, schedule_data = schedule_close_to_initial_point(
                input_data;
                optimizer = mip_optimizer,
                relax_integrality = relax_integrality,
                time_limit = scheduler_time_limit,
            )
        else
            # If we were not told to schedule independently, we schedule
            # with a 2-step procedure, where we first try with a strict
            # copper-plate balance, then with a relaxed copperplate balance
            # with penalized violation.
            #
            # Determine a valid generation/consumption schedule assuming a
            # copper-plate grid
            #
            model, schedule_data = schedule_power_copperplate(
                input_data;
                optimizer = mip_optimizer,
                relax_integrality = relax_integrality,
                time_limit = scheduler_time_limit,
                relax_balances = relax_copperplate_balances,
                include_reserves = include_reserves_in_schedule,
                warmstart_data = naive_schedule,
                overcommitment_factor = overcommitment_factor,
            )
            if schedule_data === nothing && !relax_copperplate_balances
                # If we can't find a feasible solution within the time limit, AND
                # we did not already relax the copperplate balances, relax the
                # copperplate balances and try again.
                # We may need to provide or determine a different time limit for this
                # solve.
                model, schedule_data = schedule_power_copperplate(
                    input_data;
                    optimizer = mip_optimizer,
                    relax_integrality = relax_integrality,
                    time_limit = scheduler_time_limit,
                    relax_balances = true,
                    include_reserves = include_reserves_in_schedule,
                    warmstart_data = naive_schedule,
                    overcommitment_factor = overcommitment_factor,
                )
            end
        end
        if schedule_data === nothing
            # If the scheduling procedure above (independent or otherwise) fails,
            # we simply copy scheduling decisions from the initial point to try and
            # write *some* solution file.
            schedule_data = schedule_to_initial_point(input_data)
        end
        pr_status = primal_status(model)
        feasible_schedule = (pr_status == FEASIBLE_POINT)

        if feasible_schedule
            solve_data["feasible"] = true
            solve_data["objective_value"] = objective_value(model)
        else
            solve_data["feasible"] = false
            if halt_on_infeasible_schedule
                # TODO: Find a place for this useful code
                #compute_conflict!(model)
                #println("Constraints in conflict:")
                #for con in all_constraints(model, include_variable_in_set_constraints=true)
                #    status = MOI.get(model, MOI.ConstraintConflictStatus(), con)
                #    if status == MOI.IN_CONFLICT
                #        println(con)
                #    end
                #end
                return solve_data
            end
        end

        # All the conditions required for schedule_data to have reserves
        if (
            include_reserves_in_schedule
            && feasible_schedule
            && !schedule_close_to_initial
            && sdd_fraction_to_constrain > 0.0
        )
            pos_reserve_value, neg_reserve_value = compute_device_reserve_values(input_data, schedule_data)

            topo_data = preprocess_topology_data(input_data)
            sdd2bus = Dict()
            for bus_id in input_data.bus_ids
                for sdd_id in topo_data.bus_sdd_ids[bus_id]
                    sdd2bus[sdd_id] = bus_id
                end
            end
            has_shunt = Dict(uid => (length(topo_data.bus_shunt_ids[sdd2bus[uid]]) >= 1) for uid in input_data.sdd_ids)

            println("Computing costs for fixing devices")
            # Estimate of (local) power balance violation penalty we will take by
            # *preventing the device from decreasing*
            balance_costs = compute_device_balance_costs(input_data, schedule_data, allow_switching = allow_switching)
            println("Done computing costs")

            net_value_for_lb = Dict(uid => (neg_reserve_value[uid] - balance_costs[uid]) for uid in input_data.sdd_ids)
            net_value_for_ub = Dict(uid => (pos_reserve_value[uid] - balance_costs[uid]) for uid in input_data.sdd_ids)

            sdd_by_ub_value = sort(input_data.sdd_ids, by = uid -> net_value_for_ub[uid], rev = true)
            sdd_by_lb_value = sort(input_data.sdd_ids, by = uid -> net_value_for_lb[uid], rev = true)
            sdd_by_ub_pos = [uid for uid in sdd_by_ub_value if net_value_for_ub[uid] > 0 && balance_costs[uid] == 0 && !has_shunt[uid]]
            sdd_by_lb_pos = [uid for uid in sdd_by_lb_value if net_value_for_lb[uid] > 0 && balance_costs[uid] == 0 && !has_shunt[uid]]

            n_to_ub = Int64(round(sdd_fraction_to_constrain * length(sdd_by_ub_pos)))
            n_to_lb = Int64(round(sdd_fraction_to_constrain * length(sdd_by_lb_pos)))

            sdd_to_ub = sdd_by_ub_pos[1:n_to_ub]
            sdd_to_lb = sdd_by_lb_pos[1:n_to_lb]

            println("SDD TO BE UPPER BOUNDED THROUGHOUT OPF SOLVES")
            println("---------------------------------------------")
            for (i, uid) in enumerate(sdd_to_ub)
                println("$i: uid=$uid, reserve value = $(pos_reserve_value[uid]), balance cost = $(balance_costs[uid]), net value = $(net_value_for_ub[uid])")
            end
            println("---------------------------------------------")
            println("SDD TO BE LOWER BOUNDED THROUGHOUT OPF SOLVES")
            println("---------------------------------------------")
            for (i, uid) in enumerate(sdd_to_lb)
                println("$i: uid=$uid, reserve value = $(neg_reserve_value[uid]), balance cost = $(balance_costs[uid]), net value = $(net_value_for_lb[uid])")
            end
            println("---------------------------------------------")
        else
            sdd_to_ub = []
            sdd_to_lb = []
        end

        timing_data["schedule_total"] = time() - time_mip_start

        time_opf_start = time()

        on_status = schedule_data.on_status
        real_power = schedule_data.real_power

        if tighten_bounds_pre_opf && feasible_schedule
            # Only tighten bounds if we have a feasible schedule. Otherwise bounds
            # can cross, and we don't handle that case here.
            processed_data = tighten_bounds_using_ramp_limits(
                input_data, on_status, real_power
            )
        else
            processed_data = input_data
        end

        #
        # Solve ACOPF problem for every time period
        #
        if parallel_opf
            # First, compute in parallel
            acopf_solutions = compute_opf_in_parallel(
                processed_data,
                on_status,
                real_power;
                optimizer = nlp_optimizer,
                ipopt_linear_solver = linear_solver,
                fix_real_power = fix_ac_real_power,
                penalize_power_deviation = true,
                allow_switching = allow_switching,
                resolve_rounded_shunts = resolve_rounded_shunts,
                fix_shunt_steps = !resolve_rounded_shunts,
                relax_thermal_limits = relax_thermal_limits,
                sdd_to_lb = sdd_to_lb,
                sdd_to_ub = sdd_to_ub,
            )

            acopf_solve_data = [data for (_, data) in acopf_solutions]
            #
            # Construct serializable solution data structure
            #
            solfile_dict = construct_solution_dict(
                input_data,
                schedule_data;
                opf_data = acopf_solve_data,
                write_duals = write_duals,
                include_reserves = true,
                # We need postprocessing after the parallel solves to make sure that
                # ramping constraints are satisfied.
                #
                # TODO: Why is postprocessing combined with construction of the
                # solution data structure? Just to make this function simpler?
                postprocess = true,
                print_projected_devices = print_projected_devices,
            )

            #
            # Write out solution file
            #
            if "solution_file" in keys(args)
                solution_fpath = args["solution_file"]
                # Defer making the solution file directory until now so that we hopefully
                # only make this directory if we are going to put something in it.
                mkpath(dirname(solution_fpath))
                open(io -> JSON.print(io, solfile_dict, 4), solution_fpath, "w")
                println("Wrote solution file at $(time() - start_time) s")
                println("(after start of run_ac_uc_solver)")
            end

            #
            # Now, if we used a parallel, fully-independent solve previously, we
            # re-solve with sequential OPF solves to try to achieve a better solution
            # than the solve-independently-and-project solution gives us.
            #
            if post_parallel_sequential
                acopf_solutions = compute_optimal_power_flows(
                    processed_data,
                    on_status,
                    real_power;
                    # In this post-parallel solve, we always solve sequentially.
                    sequential = true,
                    optimizer = nlp_optimizer,
                    ipopt_linear_solver = linear_solver,
                    fix_real_power = fix_ac_real_power,
                    penalize_power_deviation = true,
                    allow_switching = allow_switching,
                    resolve_rounded_shunts = resolve_rounded_shunts,
                    fix_shunt_steps = !resolve_rounded_shunts,
                    relax_thermal_limits = relax_thermal_limits,
                    sdd_to_lb = sdd_to_lb,
                    sdd_to_ub = sdd_to_ub,
                )
            end
        elseif multiperiod_opf
            if fix_ac_real_power
                throw(ArgumentError("multiperiod_opf cannot be used with fix_ac_real_power"))
            end
            if resolve_rounded_shunts
                throw(ArgumentError("multiperiod_opf does not support re-solving shunts yet"))
            end
            if relax_thermal_limits
                throw(ArgumentError("multiperiod_opf does not support relax_thermal_limits"))
            end
            _, mpacopf_solution = compute_multiperiod_opf(
                processed_data,
                on_status;
                optimizer = nlp_optimizer,
                ipopt_linear_solver = linear_solver,
                return_model = false,
            )
            # Reshape the solution data for consistency with other subroutines
            acopf_solutions = [
                (nothing, Dict(
                    key => Dict(
                        uid => Dict(
                            attr => device[attr][i] for attr in keys(device)
                        ) for (uid, device) in devicedict
                    ) for (key, devicedict) in mpacopf_solution
                )) for i in periods
            ]
        else
            # This is a non-parallel default. It can solve OPF subproblems
            # independently or sequentially. This is controlled by the sequential_opf
            # option.
            acopf_solutions = compute_optimal_power_flows(
                processed_data,
                on_status,
                real_power;
                sequential = sequential_opf,
                optimizer = nlp_optimizer,
                ipopt_linear_solver = linear_solver,
                fix_real_power = fix_ac_real_power,
                penalize_power_deviation = true,
                allow_switching = allow_switching,
                resolve_rounded_shunts = resolve_rounded_shunts,
                fix_shunt_steps = !resolve_rounded_shunts,
                relax_thermal_limits = relax_thermal_limits,
                sdd_to_lb = sdd_to_lb,
                sdd_to_ub = sdd_to_ub,
            )
        end

        timing_data["opf_total"] = time() - time_opf_start
    end # End "!simultaneous_acuc"" clause.

    acopf_solve_data = [data for (_, data) in acopf_solutions]

    # TODO: Should this printing happen in the OPF function?
    for (acopf_model, _) in acopf_solutions
        if acopf_model !== nothing
            obj_val = objective_value(acopf_model)
            ac_pr_status = primal_status(acopf_model)
            println("  $(ac_pr_status), $(obj_val)")
        end
    end

    #
    # Display results
    #
    # TODO: Should put this in a separate function, e.g.
    # display_scheduler_results
    #if pr_status == FEASIBLE_POINT && print_solution_info
    #    #
    #    # Generate local references to model variables
    #    #
    #    cost_by_period = model[:cost_by_period]
    #    p = model[:p]
    #    p_on = model[:p_on]
    #    p_on_status = model[:p_on_status]
    #    device_cost = model[:device_cost]
    #    ###

    #    println("results")

    #    println("cost: $(objective_value(model))")
    #    println("  $(value.(cost_by_period))")

    #    for uid in sdd_ids
    #        println("$(uid)")
    #        println("  p    - $(value.([p[uid, i] for i in periods]))")
    #        println("  z    - $(value.([p_on_status[uid, i] for i in periods]))")
    #        println("  cost - $(value.([device_cost[uid, i] for i in periods]))")
    #    end
    #elseif print_solution_info
    #    println("$(case_json_file): $(pr_status)")
    #end
    ###

    time_solution_start = time()

    #
    # Construct serializable dict with solution data
    #
    solfile_dict = construct_solution_dict(
        input_data,
        schedule_data;
        opf_data = acopf_solve_data,
        write_duals = write_duals,
        include_reserves = true,
        postprocess = postprocess_final_solution,
        print_projected_devices = print_projected_devices,
    )

    #
    # Write out solution file
    #
    if "solution_file" in keys(args)
        solution_fpath = args["solution_file"]
        # Defer making the solution file directory until now so that we hopefully
        # only make this directory if we are going to put something in it.
        mkpath(dirname(solution_fpath))
        open(io -> JSON.print(io, solfile_dict, 4), solution_fpath, "w")
        println("Wrote solution file at $(time() - start_time) s")
        println("(after start of run_ac_uc_solver)")
    end

    timing_data["write_solution"] = time() - time_solution_start

    # I don't think solve_data is being used for anything. Can probably
    # be removed?
    solve_data["problem"] = case_data
    solve_data["solution"] = solfile_dict
    solve_data["timing"] = timing_data

    display(timing_data)

    # TODO: Populate "run_time" in solve_data
    return solve_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--case", "-c"
            help = "the case file"
            required = true
    end

    run_ac_uc_solver(parse_args(s))
end
