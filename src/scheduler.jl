using JuMP

include("scheduling_model.jl")

function _process_schedule_data(input_data::NamedTuple, schedule_data::NamedTuple)
    periods = input_data.periods
    sdd_ids = input_data.sdd_ids
    on_status = schedule_data.on_status
    for uid in sdd_ids
        for i in periods
            on_status[uid][i] = round(on_status[uid][i])
        end
    end
    return schedule_data
end

"""
Compute a generation/consumption schedule for simple dispatchable devices
under a copper-plate assumption.

`input_data` is the NamedTuple resulting from running `process_input_data`
on the parsed input JSON file.

This function:
- Constructs a JuMP model
- Populates the JuMP model with variables, constraints, and objective
- Solves the model
- Extracts data from the model

"""
function schedule_power_copperplate(
    input_data::NamedTuple;
    optimizer = nothing,
    relax_integrality::Bool = false,
    time_limit::Float64 = 0.0,
    set_silent::Bool = false,
    relax_balances::Bool = false,
    include_reserves::Bool = false,
    relax_reserves::Bool = true,
    penalize_only_reserve_shortfalls::Bool = false,
    warmstart_data::Union{NamedTuple,Nothing} = nothing,
    overcommitment_factor::Float64 = 1.0,
)
    time_model_start = time()
    model = get_copperplate_scheduling_model(
        input_data;
        relax_balances = relax_balances,
        include_reserves = include_reserves,
        relax_reserves = relax_reserves,
        penalize_only_reserve_shortfalls = penalize_only_reserve_shortfalls,
        overcommitment_factor = overcommitment_factor,
    )

    if optimizer === nothing
        optimizer = HiGHS.Optimizer
    end

    if relax_integrality
        JuMP.relax_integrality(model)
    end

    if time_limit > 0.0
        JuMP.set_time_limit_sec(model, time_limit)
    end

    JuMP.set_optimizer(model, optimizer)

    if set_silent
        JuMP.set_silent(model)
    end

    if warmstart_data !== nothing
        warmstart!(model, input_data, warmstart_data)
    end

    time_solve_start = time()
    optimize!(model)

    if JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
        schedule_data = extract_data_from_scheduling_model(
            input_data, model;
            include_reserves = include_reserves,
        )
        schedule_data = _process_schedule_data(input_data, schedule_data)
    else
        schedule_data = nothing
    end

    return model, schedule_data
end

function schedule_power_independently(
    input_data::NamedTuple;
    optimizer = nothing,
    relax_integrality::Bool = false,
    time_limit::Float64 = 0.0,
    set_silent::Bool = false,
    include_reserves::Bool = false,
    relax_reserves::Bool = true,
    warmstart_data::Union{NamedTuple,Nothing} = nothing,
)
    model = get_copperplate_scheduling_model(
        input_data,
        include_balances = false,
        include_reserves = include_reserves,
        relax_reserves = relax_reserves,
    )

    if optimizer === nothing
        optimizer = HiGHS.Optimizer
    end

    if relax_integrality
        JuMP.relax_integrality(model)
    end

    if time_limit > 0.0
        JuMP.set_time_limit_sec(model, time_limit)
    end

    JuMP.set_optimizer(model, optimizer)

    if set_silent
        JuMP.set_silent(model)
    end

    if warmstart_data !== nothing
        warmstart!(model, input_data, warmstart_data)
    end

    optimize!(model)

    if JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
        schedule_data = extract_data_from_scheduling_model(
            input_data, model, include_reserves = include_reserves
        )
        schedule_data = _process_schedule_data(input_data, schedule_data)
    else
        schedule_data = nothing
    end

    return model, schedule_data
end

function schedule_to_initial_point(input_data::NamedTuple)
    real_power = Dict(
        uid => [
            input_data.sdd_lookup[uid]["initial_status"]["p"]
            for i in input_data.periods
        ]
        for uid in input_data.sdd_ids
    )
    on_status = Dict(
        uid => [
            input_data.sdd_lookup[uid]["initial_status"]["on_status"]
            for i in input_data.periods
        ]
        for uid in input_data.sdd_ids
    )
    return (
        real_power = real_power,
        on_status = on_status,
    )
end

"""
Solves a scheduling MIP where the cost function penalizes deviation from the
initial point.
"""
function schedule_close_to_initial_point(
    input_data::NamedTuple;
    optimizer = nothing,
    relax_integrality::Bool = false,
    time_limit::Float64 = 0.0,
    set_silent::Bool = false,
    # It may eventually be worthwhile to support reserves in this function.
    # NOTE: If reserves are supported, they must be rigorously enforced
    # (we won't include the shortfall penalties in the objective).
    #include_reserves::Bool = false,
)
    # Scheduling model we construct here should:
    # - not include objective
    # - not include any balance constraints
    # - not include reserves
    model = get_copperplate_scheduling_model(
        input_data;
        include_reserves = false,
        include_balances = false,
        include_objective = false,
    )

    sdd_ids = input_data.sdd_ids
    sdd_lookup = input_data.sdd_lookup
    periods = input_data.periods
    init_p = Dict(uid => sdd_lookup[uid]["initial_status"]["p"] for uid in sdd_ids)
    init_q = Dict(uid => sdd_lookup[uid]["initial_status"]["q"] for uid in sdd_ids)
    init_on = Dict(uid => sdd_lookup[uid]["initial_status"]["on_status"] for uid in sdd_ids)
    p_on = model[:p_on]
    q = model[:q]

    #
    # Add deviation-from initial point objective
    #
    @variable(model, 0 <= p_on_init_slack_pos[sdd_ids, periods])
    @variable(model, 0 <= p_on_init_slack_neg[sdd_ids, periods])
    @variable(model, 0 <= q_init_slack_pos[sdd_ids, periods])
    @variable(model, 0 <= q_init_slack_neg[sdd_ids, periods])

    @constraint(model, p_on_init_slack_con[uid in sdd_ids, i in periods],
        init_p[uid] * init_on[uid]
        == p_on[uid, i] + p_on_init_slack_pos[uid, i] - p_on_init_slack_neg[uid, i]
    )
    # Not entirely obvious to me whether init_q should be multiplied by init_on.
    # Also not entirely obvious that we should even consider initial q here.
    @constraint(model, q_init_slack_con[uid in sdd_ids, i in periods],
        init_q[uid] * init_on[uid]
        == q[uid, i] + q_init_slack_pos[uid, i] - q_init_slack_neg[uid, i]
    )

    @objective(model, Min,
        sum(
            p_on_init_slack_pos[uid, i] + p_on_init_slack_neg[uid, i]
            + q_init_slack_pos[uid, i] + q_init_slack_neg[uid, i]
            for uid in sdd_ids for i in periods
        )
    )
    ###

    if optimizer === nothing
        optimizer = HiGHS.Optimizer
    end

    if relax_integrality
        JuMP.relax_integrality(model)
    end

    if time_limit > 0.0
        JuMP.set_time_limit_sec(model, time_limit)
    end

    JuMP.set_optimizer(model, optimizer)

    if set_silent
        JuMP.set_silent(model)
    end

    optimize!(model)

    if JuMP.primal_status(model) == JuMP.FEASIBLE_POINT
        schedule_data = extract_data_from_scheduling_model(input_data, model)
        schedule_data = _process_schedule_data(input_data, schedule_data)
    else
        schedule_data = nothing
    end

    return model, schedule_data
end
