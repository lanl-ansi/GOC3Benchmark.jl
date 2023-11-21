using JuMP
using Ipopt
using Printf
import MathOptSymbolicAD


#include("opf_model.jl")
# : get_ac_opf_model

#include("process_bounds.jl")
# : check_bounds_not_crossing


"""
A callback to help us terminate at an acceptable point if Ipopt is taking
a long time.
"""
function ipopt_acceptable_callback(
    alg_mod,
    iter_count,
    obj_value,
    inf_pr,
    inf_du,
    mu,
    d_norm,
    regularization_size,
    alpha_du,
    alpha_pr,
    ls_trials,
)
    # If Ipopt is taking a long time, just stop at the first feasible
    # solution that is somewhat close to optimal.
    if (
        iter_count > 200
        # These appear to be unscaled infeasibilities...
        # TODO: Use GetIpoptCurrentViolations to get scaled infeasibilities?
        && inf_pr <= 1e-3
        && inf_du <= 1e0
    )
        terminate = true
    else
        terminate = false
    end
    return !terminate
end


"""
Compute optimal power flow at the specified time interval from the given input
(problem) data.

`input_data` is the `NamedTuple` returned by `process_input_data`.

This method accepts on_status and real_power in the form expected by
get_ac_opf_model.

"""
function compute_optimal_power_flow_at_interval(
    input_data::NamedTuple,
    interval::Int64;
    on_status = nothing,
    real_power = nothing,
    penalize_power_deviation::Bool = false,
    relax_power_balance::Bool = true,
    relax_p_balance::Bool = relax_power_balance,
    relax_q_balance::Bool = relax_power_balance,
    max_balance_violation = nothing,
    fix_real_power::Bool = false,
    set_silent = false,
    optimizer = nothing,
    ipopt_linear_solver = nothing,
    allow_switching = true,
    fix_shunt_steps = false,
    resolve_rounded_shunts = true,
    relax_thermal_limits = false,
    sdd_to_lb = Vector(),
    sdd_to_ub = Vector(),
)
    args = Dict{String, Any}(
        "on_status" => on_status,
        "real_power" => real_power,
        "penalize_power_deviation" => penalize_power_deviation,
        "relax_power_balance" => relax_power_balance,
        "relax_p_balance" => relax_p_balance,
        "relax_q_balance" => relax_q_balance,
        "max_balance_violation" => max_balance_violation,
        "fix_real_power" => fix_real_power,
        "allow_switching" => allow_switching,
        "fix_shunt_steps" => fix_shunt_steps,
        "relax_thermal_limits" => relax_thermal_limits,
        "sdd_to_lb" => sdd_to_lb,
        "sdd_to_ub" => sdd_to_ub,
    )
    time_model_start = time()
    model = get_ac_opf_model(
        input_data,
        interval;
        args = args,
    )

    if optimizer === nothing
        if ipopt_linear_solver === nothing
            # linear solver was not specified. We defer selection
            # of the default to Ipopt
            optimizer = JuMP.optimizer_with_attributes(
                Ipopt.Optimizer,
                "honor_original_bounds" => "yes",
                "tol" => 1e-6,
                "acceptable_tol" => 1e-4,
            )
        else
            # Use the linear solver specified
            optimizer = JuMP.optimizer_with_attributes(
                Ipopt.Optimizer,
                "honor_original_bounds" => "yes",
                "acceptable_tol" => 1e-4,
                "linear_solver" => ipopt_linear_solver,
            )
        end
    end

    JuMP.set_optimizer(model, optimizer)

    if set_silent
        JuMP.set_silent(model)
    end

    #
    # Display variables that have crossed bounds. This can occur due
    # to bugs in the bound tightening that occurs prior to OPF solves.
    #
    vars_w_crossed_bounds = Vector{JuMP.VariableRef}()
    for var in JuMP.all_variables(model)
        if JuMP.has_upper_bound(var) && JuMP.has_lower_bound(var)
            if JuMP.upper_bound(var) - JuMP.lower_bound(var) <= -1e-8
                push!(vars_w_crossed_bounds, var)
            end
        end
    end
    if !isempty(vars_w_crossed_bounds)
        println("Variables with lb > ub")
        for var in vars_w_crossed_bounds
            println(
                "  $var: ($(JuMP.lower_bound(var)), $(JuMP.upper_bound(var)))"
            )
        end
        println()
    end

    time_solve_start = time()

    # NOTE: This only makes sense if we are using Ipopt (which we are, for now)
    MOI.set(model, Ipopt.CallbackFunction(), ipopt_acceptable_callback)

    JuMP.optimize!(
        model,
        _differentiation_backend = MathOptSymbolicAD.DefaultBackend()
    )

    # Round shunt steps to integers and re-solve
    if resolve_rounded_shunts
        shunt_int_values = Dict(
            uid => round(value(model[:shunt_step][uid]))
            for uid in input_data.shunt_ids
        )
        for uid in input_data.shunt_ids
            step = model[:shunt_step][uid]
            JuMP.fix(step, shunt_int_values[uid], force = true)
            # Interestingly, once I fix this, I can no longer query the model
            # for values
        end
        JuMP.optimize!(
            model,
            _differentiation_backend = MathOptSymbolicAD.DefaultBackend()
        )
    end

    opf_solution_data = extract_data_from_model(
        model,
        input_data,
        on_status,
        real_power;
        # TODO: Find a good way to specify tolerances consistently. (I.e.
        # should this be the NLP tolerance, the integrality tolerance, or
        # something else?)
        tolerance = 1e-6,
        allow_switching = allow_switching,
    )

    return model, opf_solution_data
end


"""
Method that directly accepts input_data, schedule_data, and an integer
interval.
"""
function compute_optimal_power_flow_at_interval(
    input_data::NamedTuple,
    schedule_data::NamedTuple,
    interval::Int64;
    penalize_power_deviation::Bool = false,
    relax_power_balance::Bool = true,
    relax_p_balance::Bool = relax_power_balance,
    relax_q_balance::Bool = relax_power_balance,
    max_balance_violation = nothing,
    fix_real_power::Bool = false,
    set_silent = false,
    optimizer = nothing,
    ipopt_linear_solver = nothing,
    allow_switching = true,
    fix_shunt_steps = false,
    resolve_rounded_shunts = true,
    relax_thermal_limits = false,
    sdd_to_lb = Vector(),
    sdd_to_ub = Vector(),
)
    sdd_ids = input_data.sdd_ids
    on_status = schedule_data.on_status
    real_power = schedule_data.real_power
    interval_on_status = Dict(uid => on_status[uid][interval] for uid in sdd_ids)
    interval_real_power = Dict(uid => real_power[uid][interval] for uid in sdd_ids)
    return compute_optimal_power_flow_at_interval(
        input_data,
        interval;
        on_status = interval_on_status,
        real_power = interval_real_power,
        penalize_power_deviation = penalize_power_deviation,
        relax_power_balance = relax_power_balance,
        relax_p_balance = relax_p_balance,
        relax_q_balance = relax_q_balance,
        max_balance_violation = max_balance_violation,
        fix_real_power = fix_real_power,
        set_silent = set_silent,
        optimizer = optimizer,
        ipopt_linear_solver = ipopt_linear_solver,
        allow_switching = allow_switching,
        fix_shunt_steps = fix_shunt_steps,
        resolve_rounded_shunts = resolve_rounded_shunts,
        relax_thermal_limits = relax_thermal_limits,
        sdd_to_lb = sdd_to_lb,
        sdd_to_ub = sdd_to_ub,
    )
end


"""
Solve a square power flow model at the specified interval.
"""
function compute_power_flow_at_interval(
    input_data::NamedTuple,
    schedule_data::NamedTuple,
    interval::Int64;
    set_silent = false,
    optimizer = nothing,
    ipopt_linear_solver = nothing,
    allow_switching = true,
)
    model = get_power_flow_model_at_interval(
        input_data,
        schedule_data,
        interval;
        allow_switching = allow_switching,
    )
    # TODO: Set Optimizer attributes and solve

    sdd_ids = input_data.sdd_ids
    # Need these for extract_data_from_model
    on_status = schedule_data.on_status
    real_power = schedule_data.real_power
    interval_on_status = Dict(uid => on_status[uid][interval] for uid in sdd_ids)
    interval_on_status = Dict(uid => on_status[uid][interval] for uid in sdd_ids)
    interval_real_power = Dict(uid => real_power[uid][interval] for uid in sdd_ids)
    power_flow_data = extract_data_from_model(
        model,
        input_data,
        on_status,
        real_power;
        # TODO: Find a good way to specify tolerances consistently. (I.e.
        # should this be the NLP tolerance, the integrality tolerance, or
        # something else?)
        tolerance = 1e-6,
        allow_switching = allow_switching,
    )

    return model, power_flow_data
end


function tighten_bounds_at_interval_using_ramp_limits!(
    data,
    interval,
    previous_on_status,
    previous_real_power,
    current_on_status;
    tolerance = 1e-8,
)
    dt = data.dt[interval]
    i = interval
    for uid in data.sdd_ids
        u_prev = previous_on_status[uid]
        p_prev = previous_real_power[uid]
        u_current = current_on_status[uid]
        sdd = data.sdd_lookup[uid]
        sdd_ts = data.sdd_ts_lookup[uid]
        if round(u_prev) == 1 && round(u_current) == 1
            # Only attempt to tighten bounds if we were on in the previous
            # interval and are on in the current interval. Otherwise the ramp
            # limits do not apply.
            ramp_ub = p_prev + dt*sdd["p_ramp_up_ub"]
            ramp_lb = p_prev - dt*sdd["p_ramp_down_ub"]

            if ramp_ub < sdd_ts["p_ub"][interval] - tolerance
                sdd_ts["p_ub"][interval] = ramp_ub
            end
            if ramp_lb > sdd_ts["p_lb"][interval] + tolerance
                sdd_ts["p_lb"][interval] = ramp_lb
            end
            # Is this an appropriate way to make sure these variables exist
            # in the desired scope?
            local lb, ub
            try
                lb, ub = check_bounds_not_crossing(
                    sdd_ts["p_lb"][interval],
                    sdd_ts["p_ub"][interval],
                    tolerance = tolerance,
                )
            catch err
                println("Interval: $i")
                println("UID: $uid")
                println("p_prev: $p_prev")
                throw(err)
            end
            sdd_ts["p_lb"][interval] = lb
            sdd_ts["p_ub"][interval] = ub
            @assert ramp_ub >= ramp_lb
        end
    end
end


function compute_optimal_power_flows(
    input_data::NamedTuple,
    on_status,
    real_power;
    sequential = false,
    penalize_power_deviation::Bool = false,
    relax_power_balance::Bool = true,
    relax_p_balance::Bool = relax_power_balance,
    relax_q_balance::Bool = relax_power_balance,
    max_balance_violation = nothing,
    fix_real_power::Bool = false,
    set_silent = false,
    optimizer = nothing,
    ipopt_linear_solver = nothing,
    return_models = false,
    allow_switching = true,
    fix_shunt_steps = false,
    resolve_rounded_shunts = true,
    relax_thermal_limits = false,
    sdd_to_lb = Vector(),
    sdd_to_ub = Vector(),
)
    # TODO: Performance implications of deepcopy?
    if sequential
        # We will update this data (bounds) as we solve sequentially.
        # Copy so we don't alter caller's data structure.
        input_data = deepcopy(input_data)
    end
    sdd_key = "simple_dispatchable_device"
    i1 = first(input_data.periods)
    acopf_solutions = Vector{Tuple}()
    sdd_ids = input_data.sdd_ids
    println("BEGINNING OPF SOLVES: SEQUENTIAL = $sequential")
    for i in input_data.periods
        interval_on_status = Dict(uid => on_status[uid][i] for uid in sdd_ids)
        interval_real_power = Dict(uid => real_power[uid][i] for uid in sdd_ids)
        if i != i1 && sequential
            #
            # Adjust p_lb and p_ub in input_data, in-place. This step is what
            # makes the solve "sequential". Otherwise the solves at individual
            # time intervals are completely independent.
            #
            # TODO: This functionality will need to removed from a "subproblem
            # loop" intended to run in parallel

            # Previous on status just comes from the scheduler solution
            previous_on_status = Dict(uid => on_status[uid][i-1] for uid in sdd_ids)
            _, prev_sol_data = acopf_solutions[i-1]
            # Previous real power comes from the previous ACOPF solution
            previous_real_power = Dict(
                uid => prev_sol_data[sdd_key][uid]["p_on"] for uid in sdd_ids
            )
            tighten_bounds_at_interval_using_ramp_limits!(
                input_data,
                i,
                previous_on_status,
                previous_real_power,
                interval_on_status,
            )
        end

        # Solve ACOPF at i and extract data
        if !set_silent
            println("Starting the OPF solve at interval $i")
        end
        model, sol_data = compute_optimal_power_flow_at_interval(
            input_data,
            i;
            on_status = interval_on_status,
            real_power = interval_real_power,
            penalize_power_deviation = penalize_power_deviation,
            relax_power_balance = relax_power_balance,
            relax_p_balance = relax_p_balance,
            relax_q_balance = relax_q_balance,
            max_balance_violation = max_balance_violation,
            fix_real_power = fix_real_power,
            set_silent = set_silent,
            optimizer = optimizer,
            ipopt_linear_solver = ipopt_linear_solver,
            allow_switching = allow_switching,
            fix_shunt_steps = fix_shunt_steps,
            resolve_rounded_shunts = resolve_rounded_shunts,
            relax_thermal_limits = relax_thermal_limits,
            sdd_to_lb = sdd_to_lb,
            sdd_to_ub = sdd_to_ub,
        )

        # Add to vector of solutions
        #
        # TODO: If this loop is to run in parallel, we will need to pre-allocate
        # the acopf_solutions vector
        if return_models
            push!(acopf_solutions, (model, sol_data))
        else
            push!(acopf_solutions, (nothing, sol_data))
        end

        # TODO: OpfData NamedTuple with data, primal_status, model, etc. fields

        # Update previous status for next interval
        #previous_on_status = interval_on_status
        #previous_real_power = interval_real_power
    end

    # Return vector of (model, data)
    return acopf_solutions
end

function compute_opf_in_parallel(
    input_data::NamedTuple,
    on_status,
    real_power;
    penalize_power_deviation::Bool = false,
    relax_power_balance::Bool = true,
    relax_p_balance::Bool = relax_power_balance,
    relax_q_balance::Bool = relax_power_balance,
    max_balance_violation = nothing,
    fix_real_power::Bool = false,
    optimizer = nothing,
    ipopt_linear_solver = nothing,
    return_models = false,
    allow_switching = true,
    fix_shunt_steps = false,
    resolve_rounded_shunts = true,
    relax_thermal_limits = false,
    sdd_to_lb = Vector(),
    sdd_to_ub = Vector(),
)
    start_time = time()
    # To allow parallel solves, we pre-allocate a vector that will store the
    # result of each solve. Each item is Tuple{Nothing | JuMP.Model, NamedTuple}
    acopf_solutions = Vector{Any}([nothing for _ in input_data.periods])

    # Here, we do not need to copy the input_data NamedTuple. We do not modify
    # it at any point during this function as no bound tightening is performed.
    sdd_key = "simple_dispatchable_device"
    sdd_ids = input_data.sdd_ids
    println("BEGINNING PARALLEL OPF SOLVES")

    #
    # Parallel loop
    #
    Threads.@threads for i in input_data.periods

        # Solve ACOPF at i and extract data
        println("Starting the OPF solve at interval $i")
        interval_on_status = Dict(uid => on_status[uid][i] for uid in sdd_ids)
        interval_real_power = Dict(uid => real_power[uid][i] for uid in sdd_ids)
        model, sol_data = compute_optimal_power_flow_at_interval(
            input_data,
            i;
            on_status = interval_on_status,
            real_power = interval_real_power,
            penalize_power_deviation = penalize_power_deviation,
            relax_power_balance = relax_power_balance,
            relax_p_balance = relax_p_balance,
            relax_q_balance = relax_q_balance,
            max_balance_violation = max_balance_violation,
            fix_real_power = fix_real_power,
            # Force set_silent = true for parallel solves. Log is not meaningful
            # as messages are displayed in non-deterministic order
            set_silent = true,
            optimizer = optimizer,
            ipopt_linear_solver = ipopt_linear_solver,
            allow_switching = allow_switching,
            fix_shunt_steps = fix_shunt_steps,
            resolve_rounded_shunts = resolve_rounded_shunts,
            relax_thermal_limits = relax_thermal_limits,
            sdd_to_lb = sdd_to_lb,
            sdd_to_ub = sdd_to_ub,
        )

        # TODO: Would be nice to include these as fields in a NamedTuple in the
        # compute_opf_at_interval function.
        # TODO: Other information to include: Max p/q balance penalties, and
        # generators at which they appear.
        sol_data["primal_status"] = JuMP.primal_status(model)
        sol_data["objective_value"] = JuMP.objective_value(model)
        sol_data["solve_time"] = JuMP.solve_time(model)

        # This option exists to combat the potential memory overhead of storing all
        # models in memory simultaneously. If we do not return them from this
        # function, they are free to get GC'd. It is unclear whether this is actually
        # necessary. TODO: remove option if not necessary.
        if return_models
            acopf_solutions[i] = (model, sol_data)
        else
            acopf_solutions[i] = (nothing, sol_data)
        end

        println("Finished the OPF solve at interval $i")
    end

    println("Results of parallel OPF solves")
    println("------------------------------")
    for (i, (_, sol_data)) in enumerate(acopf_solutions)
        status = sol_data["primal_status"]
        obj = sol_data["objective_value"]
        solve_time = sol_data["solve_time"]
        println("t=$(@sprintf "%2i" i), $status, objective=$(@sprintf "%1.2e" obj), time=$(@sprintf "%3.0f" solve_time)")
    end
    println("------------------------------")

    end_time = time()
    solve_time = end_time - start_time
    println("Total time required by parallel OPF: $(@sprintf "%3.0f s" solve_time)")

    return acopf_solutions
end
