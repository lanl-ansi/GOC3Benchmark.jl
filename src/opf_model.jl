using JuMP

#include("reserves.jl")

function _filter_inactive_sdds(
    data::NamedTuple, on_status::Dict, p::Dict; tolerance = 0.0
)::NamedTuple
    inactive_sdds = Set(
        uid for uid in data.sdd_ids if (
            # The device is not on
            on_status[uid] <= tolerance
            # and the device has zero real power (could have nonzero power if
            # we are in a SU/SD curve - in this case we don't want to filter).
            #
            # ... why would uid ever not be in keys(p)? ...
            #
            # Note that this is not a reliable check. p[uid] could be nonzero
            # when u[uid] is 0 due to primal feas tol in scheduler. However, we
            # have a "structural" notion of SU/SD curves, which we should use
            # for this check instead. TODO
            && (!(uid in keys(p)) || (abs(p[uid]) <= tolerance))
        )
    )
    # Note that we are creating new vectors, not modifying existing
    # vectors in-place. This is necessary for parallel implementation.
    sdd_ids = [uid for uid in data.sdd_ids if !(uid in inactive_sdds)]
    sdd_ids_producer = [
        uid for uid in data.sdd_ids_producer if !(uid in inactive_sdds)
    ]
    sdd_ids_consumer = [
        uid for uid in data.sdd_ids_consumer if !(uid in inactive_sdds)
    ]
    return (
        sdd_ids = sdd_ids,
        sdd_ids_producer = sdd_ids_producer,
        sdd_ids_consumer = sdd_ids_consumer,
    )
end


# A method where, if an interval is provided, we treat on_status and real_power
# as containing the full time-series data from the scheduler.
function _filter_inactive_sdds(
    data::NamedTuple,
    interval::Int64,
    on_status::Dict,
    real_power::Dict;
    tolerance = 0.0,
)::NamedTuple
    i = interval
    sdd_ids = data.sdd_ids
    interval_on_status = Dict(uid => on_status[uid][i] for uid in sdd_ids)
    interval_real_power = Dict(uid => real_power[uid][i] for uid in sdd_ids)
    return _filter_inactive_sdds(
        data,
        interval_on_status,
        interval_real_power,
        tolerance = tolerance,
    )
end

function preprocess_topology_data(input_data::NamedTuple)
    # Input data that is actually necessary for this function:
    # - ac_line_lookup
    # - twt_lookup
    # - dc_line_lookup
    # - bus_ids
    # - sdd_ids (consumer/producer vectors not used)
    # - sdd_lookup
    # - shunt_lookup
    ac_line_lookup = input_data.ac_line_lookup
    twt_lookup = input_data.twt_lookup
    dc_line_lookup = input_data.dc_line_lookup
    bus_ids = input_data.bus_ids
    sdd_ids = input_data.sdd_ids
    sdd_lookup = input_data.sdd_lookup
    shunt_lookup = input_data.shunt_lookup

    #
    # Preprocess data necessary to define constraints on the network topology
    #
    # TODO: Does this need to be a set?
    branch_fr_keys = Set()
    for (k,v) in ac_line_lookup
        push!(branch_fr_keys, (k, v["fr_bus"], v["to_bus"]))
    end
    for (k,v) in twt_lookup
        push!(branch_fr_keys, (k, v["fr_bus"], v["to_bus"]))
    end
    for (k,v) in dc_line_lookup
        push!(branch_fr_keys, (k, v["fr_bus"], v["to_bus"]))
    end

    branch_to_keys = Set()
    for (uid,fr,to) in branch_fr_keys
        push!(branch_to_keys, (uid,to,fr))
    end

    # TODO: Does this need to be a set?
    branch_keys = Set([branch_fr_keys..., branch_to_keys...])

    # TODO: Do the values in this dict need to be a set?
    bus_sdd_ids = Dict(uid => Set() for uid in bus_ids)
    # Why do we use sdd_lookup here and sdd_ids below?
    # Only sdd_ids has been filtered
    for (uid, sdd) in sdd_lookup
        push!(bus_sdd_ids[sdd["bus"]], uid)
    end

    bus_sdd_producer_ids = Dict(uid => Vector() for uid in bus_ids)
    bus_sdd_consumer_ids = Dict(uid => Vector() for uid in bus_ids)
    # Iterate over sdd_ids, because this is the list we will filter
    for uid in sdd_ids
        sdd = sdd_lookup[uid]
        if sdd["device_type"]=="producer"
            push!(bus_sdd_producer_ids[sdd["bus"]], uid)
        end
        if sdd["device_type"]=="consumer"
            push!(bus_sdd_consumer_ids[sdd["bus"]], uid)
        end
    end

    bus_shunt_ids = Dict(uid => Set() for uid in bus_ids)
    for (uid, shunt) in shunt_lookup
        push!(bus_shunt_ids[shunt["bus"]], uid)
    end

    bus_branch_keys = Dict(uid => Set() for uid in bus_ids)
    for k in branch_keys
        (uid,fr,to) = k
        push!(bus_branch_keys[fr], k)
    end

    return (
        branch_fr_keys = branch_fr_keys,
        branch_to_keys = branch_to_keys,
        branch_keys = branch_keys,
        bus_sdd_ids = bus_sdd_ids,
        bus_sdd_producer_ids = bus_sdd_producer_ids,
        bus_sdd_consumer_ids = bus_sdd_consumer_ids,
        bus_shunt_ids = bus_shunt_ids,
        bus_branch_keys = bus_branch_keys,
    )
end


function add_real_reactive_linear_constraints!(model, sdd_ids, data)
    sdd_lookup = data.sdd_lookup

    q_bound = Dict(uid => sdd_lookup[uid]["q_bound_cap"] for uid in sdd_ids)
    q_linear = Dict(uid => sdd_lookup[uid]["q_linear_cap"] for uid in sdd_ids)
    pq_bound_sdds = [uid for uid in sdd_ids if q_bound[uid] == 1]
    pq_linear_sdds = [uid for uid in sdd_ids if q_linear[uid] == 1]
    q_0_ub = Dict(uid => sdd_lookup[uid]["q_0_ub"] for uid in pq_bound_sdds)
    q_0_lb = Dict(uid => sdd_lookup[uid]["q_0_lb"] for uid in pq_bound_sdds)
    beta_ub = Dict(uid => sdd_lookup[uid]["beta_ub"] for uid in pq_bound_sdds)
    beta_lb = Dict(uid => sdd_lookup[uid]["beta_lb"] for uid in pq_bound_sdds)
    q_0 = Dict(uid => sdd_lookup[uid]["q_0"] for uid in pq_linear_sdds)
    beta = Dict(uid => sdd_lookup[uid]["beta"] for uid in pq_linear_sdds)

    p = model[:p_sdd]
    q = model[:q_sdd]
    pq_ub = @constraint(
        model,
        [uid in pq_bound_sdds],
        q[uid] <= q_0_ub[uid] + beta_ub[uid]*p[uid]
    )
    pq_lb = @constraint(
        model,
        [uid in pq_bound_sdds],
        q[uid] >= q_0_lb[uid] + beta_lb[uid]*p[uid]
    )
    pq_eq = @constraint(
        model,
        [uid in pq_linear_sdds],
        q[uid] == q_0[uid] + beta[uid]*p[uid]
    )
    return pq_ub, pq_lb, pq_eq
end


function get_power_flow_model_at_interval(
    input_data,
    schedule_data,
    interval;
    allow_switching = true,
)
    vad_ub = deg2rad(30)
    vad_lb = -vad_ub

    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
     azr_lookup, azr_ts_lookup, azr_ids, rzr_lookup, rzr_ts_lookup, rzr_ids,
    ) = input_data

    #
    # Filter IDs of "inactive" SDDs
    #
    sdd_ids = input_data.sdd_ids
    on_status = schedule_data.on_status
    real_power = schedule_data.real_power
    interval_on_status = Dict(uid => on_status[uid][interval] for uid in sdd_ids)
    interval_real_power = Dict(uid => real_power[uid][interval] for uid in sdd_ids)
    sdd_ids, sdd_ids_producer, sdd_ids_consumer = _filter_inactive_sdds(
        input_data, interval_on_status, interval_real_power; tolerance = 1e-6
        # 1e-8 is the default float-int tolerance in the config file.
        # However, the more important function of this tolerance is deciding
        # when a SDD is in a SU/SD curve (as opposed to "completely off").
        # The tolerance around zero that should be "off" is determined by the
        # primal feasibility tolerance in the MIP solver. For Gurobi, the
        # default is 1e-6. (For HiGHS, the default is 1e-7.)
    )
    ###

    #
    # Filter IDs of inactive AC lines
    #
    if !allow_switching
        ac_line_ids = [
            uid for uid in ac_line_ids
            if ac_line_lookup[uid]["initial_status"]["on_status"] == 1
        ]
        twt_ids = [
            uid for uid in twt_ids
            if twt_lookup[uid]["initial_status"]["on_status"] == 1
        ]
        # We frequently iterate over these dicts, so we need to filter
        # them as well.
        ac_line_lookup = Dict(uid=>ac_line_lookup[uid] for uid in ac_line_ids)
        twt_lookup = Dict(uid=>twt_lookup[uid] for uid in twt_ids)
    end

    model = JuMP.Model()
    return model
end

"""
given data and a time period, solves an AC-OPF model for that period
presumes a unit commitment schedule has been previously determined
e.g., from the scheduling model

    data - dictionary obtained from loading the input json file
    period - time interval for which the AC-OPF model is constructed
    args - dictionary of additional arguments
        "on_status" => dictionary mapping SDD UIDs to binary on/off status
        "real_power" => dictionary mapping SDD UIDs to power produced or consumed

"""
function get_ac_opf_model(data::NamedTuple, period::Int; args=nothing)
    if args === nothing
        args = Dict()
    end
    print_program_info = get(args, "print_program_info", false)

    fix_real_power = get(args, "fix_real_power", false)
    relax_power_balance = get(args, "relax_power_balance", true)
    # relax_power_balance serves as a default for relax_p/q_balance.
    relax_p_balance = get(args, "relax_p_balance", relax_power_balance)
    relax_q_balance = get(args, "relax_q_balance", relax_power_balance)
    penalize_power_deviation = get(args, "penalize_power_deviation", false)
    max_balance_violation = get(args, "max_balance_violation", nothing)
    allow_switching = get(args, "allow_switching", true)
    fix_shunt_steps = get(args, "fix_shunt_steps", false)
    relax_thermal_limits = get(args, "relax_thermal_limits", false)
    sdd_to_lb = get(args, "sdd_to_lb", Vector())
    sdd_to_ub = get(args, "sdd_to_ub", Vector())

    ### Data Preparation ###
    if print_program_info
        println("data preparation")
        println("  time period: $(period)")
    end

    vad_ub = deg2rad(30)
    vad_lb = -vad_ub

    processed_data = data
    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
     azr_lookup, azr_ts_lookup, azr_ids, rzr_lookup, rzr_ts_lookup, rzr_ids,
    ) = processed_data

    #
    # If not provided, create default dicts that do not filter any SDDs
    #
    on_status_dict = get(
        args, "on_status", Dict(uid => 1 for uid in sdd_ids)
    )
    # By default, do not populate
    p_dict = get(args, "real_power", Dict())
    ###

    #
    # Filter IDs of "inactive" SDDs
    #
    # Note that this should not modify any of the Dicts/Vectors in
    # processed_data in-place.
    sdd_ids, sdd_ids_producer, sdd_ids_consumer = _filter_inactive_sdds(
        processed_data, on_status_dict, p_dict; tolerance = 1e-6
        # 1e-8 is the default float-int tolerance in the config file.
        # However, the more important function of this tolerance is deciding
        # when a SDD is in a SU/SD curve (as opposed to "completely off").
        # The tolerance around zero that should be "off" is determined by the
        # primal feasibility tolerance in the MIP solver. For Gurobi, the
        # default is 1e-6. (For HiGHS, the default is 1e-7.)
    )
    ###

    #
    # Filter IDs of inactive AC lines
    #
    if !allow_switching
        ac_line_ids = [
            uid for uid in ac_line_ids
            if ac_line_lookup[uid]["initial_status"]["on_status"] == 1
        ]
        twt_ids = [
            uid for uid in twt_ids
            if twt_lookup[uid]["initial_status"]["on_status"] == 1
        ]
        # We frequently iterate over these dicts, so we need to filter
        # them as well.
        # TODO: Only iterate over the arrays rather than the dicts
        ac_line_lookup = Dict(uid=>ac_line_lookup[uid] for uid in ac_line_ids)
        twt_lookup = Dict(uid=>twt_lookup[uid] for uid in twt_ids)
    end

    filtered_data = (
        sdd_ids = sdd_ids,
        sdd_ids_producer = sdd_ids_producer,
        sdd_ids_consumer = sdd_ids_consumer,
        ac_line_ids = ac_line_ids,
        ac_line_lookup = ac_line_lookup,
        twt_ids = twt_ids,
        twt_lookup = twt_lookup,
        dc_line_lookup = dc_line_lookup,
        bus_ids = bus_ids,
        sdd_lookup = sdd_lookup,
        shunt_lookup = shunt_lookup,
    )
    #
    # Preprocess data necessary to define constraints on the network topology
    #
    # Don't just use input_data here as we want to capture the filtered
    # arrays of UIDs
    topo_data = preprocess_topology_data(filtered_data)

    branch_fr_keys = topo_data.branch_fr_keys
    branch_to_keys = topo_data.branch_to_keys
    branch_keys = topo_data.branch_keys
    bus_sdd_ids = topo_data.bus_sdd_ids
    bus_sdd_producer_ids = topo_data.bus_sdd_producer_ids
    bus_sdd_consumer_ids = topo_data.bus_sdd_consumer_ids
    bus_shunt_ids = topo_data.bus_shunt_ids
    bus_branch_keys = topo_data.bus_branch_keys

    ### Mathematical Model ###
    if print_program_info
        println("model building")
    end

    model = JuMP.Model()

    t = period

    @variable(model, 
        bus_lookup[uid]["vm_lb"] <= vm[uid in bus_ids] <= bus_lookup[uid]["vm_ub"], start=1.0
    )
    @variable(model, va[uid in bus_ids])

    # TODO: Add a parameter to tighten these bounds. This will require a dict of
    # scheduled values for real power.
    @variable(model, 
        sdd_ts_lookup[uid]["p_lb"][t] <= p_sdd[uid in sdd_ids] <= sdd_ts_lookup[uid]["p_ub"][t]
    )
    @variable(model, 
        sdd_ts_lookup[uid]["q_lb"][t] <= q_sdd[uid in sdd_ids] <= sdd_ts_lookup[uid]["q_ub"][t]
    )

    pq_ub, pq_lb, pq_eq = add_real_reactive_linear_constraints!(
        model, sdd_ids, processed_data
    )

    @variable(model, p_branch[uid in branch_keys])
    @variable(model, q_branch[uid in branch_keys])

    if relax_thermal_limits
        @variable(model, ac_thermal_slack[uid in filtered_data.ac_line_ids] >= 0.0)
        @variable(model, twt_thermal_slack[uid in filtered_data.twt_ids] >= 0.0)
    end

    #
    # Continuous shunt_step variable. This will be rounded during
    # post-processing.
    #
    @variable(
        model,
        shunt_lookup[uid]["step_lb"] <= shunt_step[uid in shunt_ids] <= shunt_lookup[uid]["step_ub"],
        start = shunt_lookup[uid]["initial_status"]["step"],
    )
    if fix_shunt_steps
        for uid in shunt_ids
            init_step = shunt_lookup[uid]["initial_status"]["step"]
            JuMP.fix(shunt_step[uid], init_step, force = true)
        end
    end

    #
    # Fix variables that are specified in the "real power" argument
    # Without relaxing power balances, I expect this to be infeasible
    #
    if fix_real_power
        for uid in sdd_ids
            if uid in keys(p_dict) # Hopefully calling keys() is fast...
                # Note that this overrides the bounds of this variable
                #
                # ... why would any UID in sdd_ids (filtered already)
                # not be in p_dict ...
                JuMP.fix(p_sdd[uid], p_dict[uid]; force = true)
            end
        end
    else
        # Otherwise, we only fix "active" SDDs with on-status == 0.
        # These are the SDDs that are currently in SU/SD.
        #
        # Note that sdd_ids has already been filtered
        # TODO: Allow specification of tolerance somewhere
        on_sdds = [uid for uid in sdd_ids if abs(on_status_dict[uid] - 1.0) <= 1e-6]
        susd_sdds = [uid for uid in sdd_ids if abs(on_status_dict[uid]) <= 1e-6]
        @assert length(on_sdds) + length(susd_sdds) == length(sdd_ids)

        # Fix p[uid] for those in SU/SD curves
        for uid in susd_sdds
            JuMP.fix(p_sdd[uid], p_dict[uid]; force = true)
        end

        # Fix SDDs that we were told to fix
        to_ub_set = Set(sdd_to_ub)
        to_lb_set = Set(sdd_to_lb)
        for uid in sdd_ids
            if uid in to_lb_set
                if JuMP.has_upper_bound(p_sdd[uid])
                    ub = JuMP.upper_bound(p_sdd[uid])
                    if ub < p_dict[uid]
                        #println("WARNING: p_sdd[$uid] has an upper bound $ub below its scheduled value of $(p_dict[uid])")
                        JuMP.set_upper_bound(p_sdd[uid], p_dict[uid])
                    end
                end
                if !JuMP.is_fixed(p_sdd[uid])
                    JuMP.set_lower_bound(p_sdd[uid], p_dict[uid])
                end
            end
            if uid in to_ub_set
                if JuMP.has_lower_bound(p_sdd[uid])
                    lb = JuMP.lower_bound(p_sdd[uid])
                    if lb > p_dict[uid]
                        #println("WARNING: p_sdd[$uid] has a lower bound $lb above its scheduled value of $(p_dict[uid])")
                        JuMP.set_lower_bound(p_sdd[uid], p_dict[uid])
                    end
                end
                if !JuMP.is_fixed(p_sdd[uid])
                    JuMP.set_upper_bound(p_sdd[uid], p_dict[uid])
                end
            end
        end

        # Penalize deviation from scheduled power
        if penalize_power_deviation
            p_slack_pos = @variable(model, p_slack_pos[on_sdds] >= 0)
            p_slack_neg = @variable(model, p_slack_neg[on_sdds] >= 0)

            # Why is the penalty coefficient p-bus-vio-cost? Deviating
            # from scheduled power does not lead to a balance violation,
            # as we include the balance constraints in this model.
            # What are we trying to correct for with this penalty?
            # Violation of inter-temporal constraints and reserve constraints
            # that may have determined the power levels in the scheduling model.
            # TODO: one of these penalty coefficients should be used instead...
            #
            sched_dev_penalty = Dict(
                uid => 1e-2*violation_cost["p_bus_vio_cost"] for uid in sdd_ids
            )
            sched_dev_eqn = @constraint(
                model,
                sched_dev_eqn[uid in on_sdds],
                p_sdd[uid] - p_dict[uid] == p_slack_pos[uid] - p_slack_neg[uid],
            )
        end
    end
    ###

    #
    # Add slack variables and define penalty parameters for bus power balance
    # Note that these are only used if relax_power_balance is true
    #
    if relax_p_balance
        if max_balance_violation !== nothing
            p_balance_slack_pos = @variable(
                model,
                0 <= p_balance_slack_pos[bus_ids] <= max_balance_violation,
            )
            p_balance_slack_neg = @variable(
                model,
                0 <= p_balance_slack_neg[bus_ids] <= max_balance_violation,
            )
        else
            p_balance_slack_pos = @variable(model, 0 <= p_balance_slack_pos[bus_ids])
            p_balance_slack_neg = @variable(model, 0 <= p_balance_slack_neg[bus_ids])
        end
        p_balance_penalty = violation_cost["p_bus_vio_cost"]
    end
    if relax_q_balance
        if max_balance_violation !== nothing
            q_balance_slack_pos = @variable(
                model,
                0 <= q_balance_slack_pos[bus_ids] <= max_balance_violation,
            )
            q_balance_slack_neg = @variable(
                model,
                0 <= q_balance_slack_neg[bus_ids] <= max_balance_violation,
            )
        else
            q_balance_slack_pos = @variable(model, 0 <= q_balance_slack_pos[bus_ids])
            q_balance_slack_neg = @variable(model, 0 <= q_balance_slack_neg[bus_ids])
        end
        q_balance_penalty = violation_cost["q_bus_vio_cost"]
    end
    ###

    # angle refrence constraint
    for (uid,bus) in bus_lookup
        @constraint(model, va[uid] == 0.0)
        break
    end

    #
    # Power balance constraints
    #
    if relax_p_balance
        @NLconstraint(model, 
            p_balance[uid in bus_ids],
            sum(p_branch[k] for k in bus_branch_keys[uid], init = 0) ==
            sum(p_sdd[ssd_id] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(p_sdd[ssd_id] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) -
            sum(
                shunt_lookup[shunt_id]["gs"]*shunt_step[shunt_id]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid]^2
            #gs*vm[uid]^2
            # Add positive and negative slack variables
            + (p_balance_slack_pos[uid] - p_balance_slack_neg[uid])
        )
    else
        @NLconstraint(model, 
            p_balance[uid in bus_ids],
            sum(p_branch[k] for k in bus_branch_keys[uid], init = 0) ==
            sum(p_sdd[ssd_id] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(p_sdd[ssd_id] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) -
            sum(
                shunt_lookup[shunt_id]["gs"]*shunt_step[shunt_id]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid]^2
            #gs*vm[uid]^2
        )
    end
    if relax_q_balance
        @NLconstraint(model, 
            q_balance[uid in bus_ids],
            sum(q_branch[k] for k in bus_branch_keys[uid], init = 0) ==
            sum(q_sdd[ssd_id] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(q_sdd[ssd_id] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) +
            sum(
                shunt_lookup[shunt_id]["bs"]*shunt_step[shunt_id]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid]^2
            #bs*vm[uid]^2
            # Add positive and negative slack variables
            + (q_balance_slack_pos[uid] - q_balance_slack_neg[uid])
        )
    else
        @NLconstraint(model, 
            q_balance[uid in bus_ids],
            sum(q_branch[k] for k in bus_branch_keys[uid], init = 0) ==
            sum(q_sdd[ssd_id] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(q_sdd[ssd_id] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) +
            sum(
                shunt_lookup[shunt_id]["bs"]*shunt_step[shunt_id]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid]^2
            #bs*vm[uid]^2
        )
    end

    for (uid,ac_line) in ac_line_lookup
        branch = ac_line
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key]
        q_fr = q_branch[fr_key]
        p_to = p_branch[to_key]
        q_to = q_branch[to_key]

        if relax_thermal_limits
            s_thermal = ac_thermal_slack[uid]
        else
            s_thermal = 0.0
        end

        vm_fr = vm[branch["fr_bus"]]
        vm_to = vm[branch["to_bus"]]
        va_fr = va[branch["fr_bus"]]
        va_to = va[branch["to_bus"]]

        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        tm = 1.0
        ta = 0.0
        tr = tm * cos(ta)
        ti = tm * sin(ta)
        ttm = tm^2

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        # variable bounds
        JuMP.set_lower_bound(p_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(p_to, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_to, -10.0*branch["mva_ub_nom"])

        JuMP.set_upper_bound(p_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(p_to,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_to,  10.0*branch["mva_ub_nom"])

        # From side of the branch flow
        JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, vad_lb <= va_fr - va_to <= vad_ub)

        # Apparent power limit, from side and to side
        # TODO: What is the best way to formulate this as a soft constraint?
        # Square root? Or (s_max + slack)^2?
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= (branch["mva_ub_nom"] + s_thermal)^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= (branch["mva_ub_nom"] + s_thermal)^2)
    end


    for (uid,twt) in twt_lookup
        branch = twt
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key]
        q_fr = q_branch[fr_key]
        p_to = p_branch[to_key]
        q_to = q_branch[to_key]

        if relax_thermal_limits
            s_thermal = twt_thermal_slack[uid]
        else
            s_thermal = 0.0
        end

        vm_fr = vm[branch["fr_bus"]]
        vm_to = vm[branch["to_bus"]]
        va_fr = va[branch["fr_bus"]]
        va_to = va[branch["to_bus"]]

        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        initial_status = branch["initial_status"]

        tm = initial_status["tm"]
        ta = initial_status["ta"]
        tr = tm * cos(ta)
        ti = tm * sin(ta)
        ttm = tm^2

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        # variable bounds
        JuMP.set_lower_bound(p_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(p_to, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_to, -10.0*branch["mva_ub_nom"])

        JuMP.set_upper_bound(p_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(p_to,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_to,  10.0*branch["mva_ub_nom"])

        # From side of the branch flow
        JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, vad_lb <= va_fr - va_to <= vad_ub)

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= (branch["mva_ub_nom"] + s_thermal)^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= (branch["mva_ub_nom"] + s_thermal)^2)
    end


    for (uid,dc_line) in dc_line_lookup
        branch = dc_line
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key]
        q_fr = q_branch[fr_key]
        p_to = p_branch[to_key]
        q_to = q_branch[to_key]


        dc_line["pminf"] = -dc_line["pdc_ub"]
        dc_line["pmaxf"] = dc_line["pdc_ub"]
        dc_line["pmint"] = -dc_line["pdc_ub"]
        dc_line["pmaxt"] = dc_line["pdc_ub"]

        dc_line["qminf"] = dc_line["qdc_fr_lb"]
        dc_line["qmaxf"] = dc_line["qdc_fr_ub"]
        dc_line["qmint"] = dc_line["qdc_to_lb"]
        dc_line["qmaxt"] = dc_line["qdc_to_ub"]


        # variable bounds
        JuMP.set_lower_bound(p_fr, -branch["pdc_ub"])
        JuMP.set_lower_bound(q_fr,  branch["qdc_fr_lb"])
        JuMP.set_lower_bound(p_to, -branch["pdc_ub"])
        JuMP.set_lower_bound(q_to,  branch["qdc_to_lb"])

        JuMP.set_upper_bound(p_fr, branch["pdc_ub"])
        JuMP.set_upper_bound(q_fr, branch["qdc_fr_ub"])
        JuMP.set_upper_bound(p_to, branch["pdc_ub"])
        JuMP.set_upper_bound(q_to, branch["qdc_to_ub"])

        # From side of the branch flow
        JuMP.@constraint(model, p_fr + p_to == 0.0)
    end

    ### Mathematical Model - Objective ###

    device_cost = Dict()
    for uid in sdd_ids
        cost_blocks = sdd_ts_lookup[uid]["cost"][t]
        cost_block_p = @variable(model,
            [i in 1:length(cost_blocks)],
            lower_bound = 0.0,
            upper_bound = cost_blocks[i][2]
        )
        device_cost[uid] = @expression(model,
            sum(
                cb[1]*cost_block_p[i] for (i,cb) in enumerate(cost_blocks),
                init = 0
            )
        )
        JuMP.@constraint(model, p_sdd[uid] == sum(cost_block_p, init = 0))
    end

    # Compute a term penalizing deviation from scheduled power.
    if fix_real_power || !penalize_power_deviation
        # If we fixed real power, this term is constant, so we ignore it.
        # This overrides the penalize_power_deviation flag.
        p_sched_penalty_term = 0.0
    else
        # If we were instructed to penalize power deviation from scheduled,
        # and have not fixed real power, create a term penalizing the slacks
        # applied to p == p_sched.
        p_sched_penalty_term = @expression(model,
            p_sched_penalty_term,
            sum(
                sched_dev_penalty[uid]*(p_slack_pos[uid] + p_slack_neg[uid])
                for uid in on_sdds,
                # NOTE: "init = 0," gives a syntax error. Why is this trailing
                # comma not allowed?
                init = 0
            )
        )
    end
    if relax_p_balance
        p_balance_penalty_term = @expression(model,
            p_balance_penalty*sum(
                p_balance_slack_pos[uid] + p_balance_slack_neg[uid]
                for uid in bus_ids,
                init = 0
            )
        )
    else
        p_balance_penalty_term = 0.0
    end
    if relax_q_balance
        q_balance_penalty_term = @expression(model,
            q_balance_penalty*sum(
                q_balance_slack_pos[uid] + q_balance_slack_neg[uid]
                for uid in bus_ids,
                init = 0
            )
        )
    else
        q_balance_penalty_term = 0.0
    end
    # TODO: Option to make thermal penalties hard constraints?
    if relax_thermal_limits
        @expression(model, ac_thermal_penalty_term,
            violation_cost["s_vio_cost"]*sum(
                ac_thermal_slack[uid] for uid in filtered_data.ac_line_ids
            )
        )
        @expression(model, twt_thermal_penalty_term,
            violation_cost["s_vio_cost"]*sum(
                twt_thermal_slack[uid] for uid in filtered_data.twt_ids
            )
        )
    else
        ac_thermal_penalty_term = 0.0
        twt_thermal_penalty_term = 0.0
    end
    @objective(model, Max, dt[t]*(
            sum(device_cost[uid] for uid in sdd_ids_consumer, init = 0) -
            sum(device_cost[uid] for uid in sdd_ids_producer, init = 0)
            - p_sched_penalty_term
            - p_balance_penalty_term
            - q_balance_penalty_term
            - ac_thermal_penalty_term
            - twt_thermal_penalty_term
        )
    )
    #end

    #println(model)

    return model
end


function get_multiperiod_acopf_model(
    data::NamedTuple;
    args=nothing,
    optimizer=nothing,
)
    if args === nothing
        args = Dict{String, Any}()
    end
    # These are the options supported by get_ac_opf_model. The exact options
    # we will support here are TBD.
    print_program_info = get(args, "print_program_info", false)
    fix_real_power = get(args, "fix_real_power", false)
    relax_power_balance = get(args, "relax_power_balance", true)
    # relax_power_balance serves as a default for relax_p/q_balance.
    relax_p_balance = get(args, "relax_p_balance", relax_power_balance)
    relax_q_balance = get(args, "relax_q_balance", relax_power_balance)
    penalize_power_deviation = get(args, "penalize_power_deviation", false)
    max_balance_violation = get(args, "max_balance_violation", nothing)
    allow_switching = get(args, "allow_switching", true)
    fix_shunt_steps = get(args, "fix_shunt_steps", false)
    #relax_thermal_limits = get(args, "relax_thermal_limits", false)
    sdd_to_lb = get(args, "sdd_to_lb", Vector())
    sdd_to_ub = get(args, "sdd_to_ub", Vector())

    if (
        fix_real_power
        || penalize_power_deviation
        || !allow_switching
        || !isempty(sdd_to_lb)
        || !isempty(sdd_to_ub)
    )
        throw(ArgumentError(
            "Invalid argument. fix_real_power, penalize_power_deviation,
            allow_switching, sdd_to_lb, and sdd_to_ub should not be used
            by get_multiperiod_acopf_model"
        ))
    end

    vad_ub = deg2rad(30)
    vad_lb = -vad_ub
    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
     azr_lookup, azr_ts_lookup, azr_ids, rzr_lookup, rzr_ts_lookup, rzr_ids,
    ) = data

    # If not provided, create default dicts that do not filter any SDDs
    #
    # TODO: on_status is basically required. Maybe it should be a positional
    # arg?
    #
    # on_status is used to determine which devices should have their power levels
    # fixed (i.e. are in su/sd curves). We use a dict that is ducktype-compatible
    # with on_status variables from the UC problem
    if "on_status" in keys(args)
        # We store the original on-status, as this data structure is required for
        # a few data processing functions
        orig_on_status = args["on_status"]
        on_status = Dict((uid, i) => args["on_status"][uid][i] for uid in sdd_ids for i in periods)
    else
        orig_on_status = Dict(uid => [1 for i in periods] for uid in sdd_ids)
        on_status = Dict((uid, i) => 1 for uid in sdd_ids for i in periods)
    end

    if !allow_switching
        println("switching must be allowed in multiperiod ACOPF")
        throw(Exception)
    end
    ac_on_status_dict = Dict(uid => [1 for _ in periods] for uid in ac_line_ids)
    twt_on_status_dict = Dict(uid => [1 for _ in periods] for uid in twt_ids)

    topo_data = preprocess_topology_data(data)
    branch_fr_keys = topo_data.branch_fr_keys
    branch_to_keys = topo_data.branch_to_keys
    branch_keys = topo_data.branch_keys
    bus_sdd_ids = topo_data.bus_sdd_ids
    bus_sdd_producer_ids = topo_data.bus_sdd_producer_ids
    bus_sdd_consumer_ids = topo_data.bus_sdd_consumer_ids
    bus_shunt_ids = topo_data.bus_shunt_ids
    bus_branch_keys = topo_data.bus_branch_keys

    if optimizer === nothing
        model = JuMP.Model()
    else
        model = JuMP.Model(optimizer)
    end

    # In a loop, construct variables and constraints for ACOPF at every point
    # in time.
    @variable(model, bus_lookup[uid]["vm_lb"] <= vm[uid in bus_ids, t in periods] <= bus_lookup[uid]["vm_ub"], start=1.0)
    @variable(model, va[uid in bus_ids, t in periods])

    # Since we will have on-status implication constraints, the domain of p/q-sdd
    # variables always need to include zero.
    @variable(model,
        min(sdd_ts_lookup[uid]["p_lb"][t], 0.0)
        <= p_sdd[uid in sdd_ids, t in periods]
        <= max(sdd_ts_lookup[uid]["p_ub"][t], 0.0)
    )
    @variable(model,
        min(sdd_ts_lookup[uid]["q_lb"][t], 0.0)
        <= q_sdd[uid in sdd_ids, t in periods]
        <= max(sdd_ts_lookup[uid]["q_ub"][t], 0.0)
    )

    @variable(model, p_branch[uid in branch_keys, t in periods])
    @variable(model, q_branch[uid in branch_keys, t in periods])
    # TODO: Support thermal limit relaxation if necessary
    #if relax_thermal_limits
    #    @variable(model, ac_thermal_slack[uid in ac_line_ids, t in periods] >= 0.0)
    #    @variable(model, twt_thermal_slack[uid in twt_ids, t in periods] >= 0.0)
    #end

    T_su_pc, T_sd_pc = get_inverse_power_curve_intervals(data)
    # With fixed on status, we can compute u_su/sd
    # TODO: don't attempt to compute these if variables are provided (rather than
    # constants)
    u_su, u_sd = get_su_sd_from_on_status(data, orig_on_status)
    # Create dicts that are ducktype-compatible with JuMP variables
    u_su = Dict((uid, i) => u_su[uid][i] for uid in sdd_ids for i in periods)
    u_sd = Dict((uid, i) => u_sd[uid][i] for uid in sdd_ids for i in periods)

    add_real_reactive_linking_constraints!(
        model, data, T_su_pc, T_sd_pc;
        pq=(p_sdd, q_sdd), u_on=on_status, u_susd=(u_su, u_sd), include_reserves = false
    )

    in_supc, p_su, in_sdpc, p_sd = get_supc_sdpc_lookups(data, orig_on_status)

    # Fix power for devices in power curves
    # NOTE: This will need to change if on-status is not fixed
    for uid in sdd_ids
        for i in periods
            if Bool(in_supc[uid][i])
                JuMP.fix(p_sdd[uid, i], p_su[uid][i], force = true)
            end
            if Bool(in_sdpc[uid][i])
                JuMP.fix(p_sdd[uid, i], p_sd[uid][i], force = true)
            end
            if !(Bool(in_supc[uid][i]) || Bool(in_sdpc[uid][i]) || Bool(on_status[uid, i]))
                # If we're not online or in a power curve, fix power to zero
                JuMP.fix(p_sdd[uid, i], 0.0, force = true)
            end
        end
    end

    # This function should be used if we have distinct variables for p_on/su/sd
    #add_on_su_sd_implication_constraints!(model, data; include_reserves=false,
    #    p_on=mock_p_on, u_on=on_status, p_su=mock_p_su, p_sd=mock_p_sd,
    #)
    # We implement our own implication constraints here, as we only need a subset
    # of those from the scheduling_model
    online_device_periods = [(uid, i) for uid in sdd_ids for i in periods if Bool(on_status[uid, i])]
    p_on_implication = @constraint(model, [(uid, i) in online_device_periods],
        sdd_ts_lookup[uid]["p_lb"][i] * on_status[uid, i]
        <= p_sdd[uid, i]
        <= sdd_ts_lookup[uid]["p_ub"][i] * on_status[uid, i]
    )

    add_reactive_power_implication_constraints!(model, data, T_su_pc, T_sd_pc;
        include_reserves=false, q=q_sdd, u_on=on_status, u_su=u_su, u_sd=u_sd,
    )

    # TODO: Integer domain and explicitly relax?
    @variable(
        model,
        shunt_lookup[uid]["step_lb"]
        <= shunt_step[uid in shunt_ids, i in periods]
        <= shunt_lookup[uid]["step_ub"]
    )
    if fix_shunt_steps
        for uid in shunt_ids
            init_step = shunt_lookup[uid]["initial_status"]["step"]
            JuMP.fix.(shunt_step[uid, :], init_step, force = true)
        end
    end

    # Add slack variables for power balance relaxation (which we usually use)
    if relax_p_balance
        if max_balance_violation !== nothing
            p_balance_slack_pos = @variable(
                model,
                0 <= p_balance_slack_pos[bus_ids, periods] <= max_balance_violation,
            )
            p_balance_slack_neg = @variable(
                model,
                0 <= p_balance_slack_neg[bus_ids, periods] <= max_balance_violation,
            )
        else
            p_balance_slack_pos = @variable(model, 0 <= p_balance_slack_pos[bus_ids, periods])
            p_balance_slack_neg = @variable(model, 0 <= p_balance_slack_neg[bus_ids, periods])
        end
        p_balance_penalty = violation_cost["p_bus_vio_cost"]
    end
    if relax_q_balance
        if max_balance_violation !== nothing
            q_balance_slack_pos = @variable(
                model,
                0 <= q_balance_slack_pos[bus_ids, periods] <= max_balance_violation,
            )
            q_balance_slack_neg = @variable(
                model,
                0 <= q_balance_slack_neg[bus_ids, periods] <= max_balance_violation,
            )
        else
            q_balance_slack_pos = @variable(model, 0 <= q_balance_slack_pos[bus_ids, periods])
            q_balance_slack_neg = @variable(model, 0 <= q_balance_slack_neg[bus_ids, periods])
        end
        q_balance_penalty = violation_cost["q_bus_vio_cost"]
    end

    for (uid, bus) in bus_lookup
        @constraint(model, [i in periods], va[uid, i] == 0.0)
        break
    end

    # Implement piecewise-linear cost functions
    device_cost = Dict()
    for uid in sdd_ids
        for t in periods
            cost_blocks = sdd_ts_lookup[uid]["cost"][t]
            cost_block_p = @variable(model,
                [i in 1:length(cost_blocks)],
                lower_bound = 0.0,
                upper_bound = cost_blocks[i][2]
            )
            device_cost[uid, t] = @expression(model,
                sum(
                    cb[1]*cost_block_p[i] for (i, cb) in enumerate(cost_blocks),
                    init = 0
                )
            )
            JuMP.@constraint(model, p_sdd[uid, t] == sum(cost_block_p, init = 0))
        end
    end

    # Create penalty terms for power balance violations
    if relax_p_balance
        p_balance_penalty_term = @expression(model,
            [i in periods],
            p_balance_penalty*sum(
                p_balance_slack_pos[uid, i] + p_balance_slack_neg[uid, i]
                for uid in bus_ids,
                init = 0
            )
        )
    else
        p_balance_penalty_term = 0.0
    end
    if relax_q_balance
        q_balance_penalty_term = @expression(model,
            [i in periods],
            q_balance_penalty*sum(
                q_balance_slack_pos[uid, i] + q_balance_slack_neg[uid, i]
                for uid in bus_ids,
                init = 0
            )
        )
    else
        q_balance_penalty_term = 0.0
    end
    # TODO: Thermal limit penalty expressions?
    # (Only if I allow relaxing these constraints.)

    @objective(model, Max, sum(
        dt[t] * (
            sum(device_cost[uid, t] for uid in sdd_ids_consumer, init = 0)
            - sum(device_cost[uid, t] for uid in sdd_ids_producer, init = 0)
            - p_balance_penalty_term[t]
            - q_balance_penalty_term[t]
        ) for t in periods
    ))

    # Constraints
    # Note that we don't include max-energy-over-interval constraints/penalties,
    # although they could be relevant here.

    if relax_p_balance
        @NLconstraint(model,
            p_balance[uid in bus_ids, i in periods],
            sum(p_branch[k, i] for k in bus_branch_keys[uid], init = 0) ==
            sum(p_sdd[ssd_id, i] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(p_sdd[ssd_id, i] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) -
            sum(
                shunt_lookup[shunt_id]["gs"]*shunt_step[shunt_id, i]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid, i]^2
            #gs*vm[uid]^2
            # Add positive and negative slack variables
            + (p_balance_slack_pos[uid, i] - p_balance_slack_neg[uid, i])
        )
    else
        @NLconstraint(model,
            p_balance[uid in bus_ids, i in periods],
            sum(p_branch[k, i] for k in bus_branch_keys[uid], init = 0) ==
            sum(p_sdd[ssd_id, i] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(p_sdd[ssd_id, i] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) -
            sum(
                shunt_lookup[shunt_id]["gs"]*shunt_step[shunt_id, i]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid, i]^2
            #gs*vm[uid]^2
        )
    end
    if relax_q_balance
        @NLconstraint(model,
            q_balance[uid in bus_ids, i in periods],
            sum(q_branch[k, i] for k in bus_branch_keys[uid], init = 0) ==
            sum(q_sdd[ssd_id, i] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(q_sdd[ssd_id, i] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) +
            sum(
                shunt_lookup[shunt_id]["bs"]*shunt_step[shunt_id, i]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid, i]^2
            #bs*vm[uid]^2
            # Add positive and negative slack variables
            + (q_balance_slack_pos[uid, i] - q_balance_slack_neg[uid, i])
        )
    else
        @NLconstraint(model, 
            q_balance[uid in bus_ids, i in periods],
            sum(q_branch[k, i] for k in bus_branch_keys[uid], init = 0) ==
            sum(q_sdd[ssd_id, i] for ssd_id in bus_sdd_producer_ids[uid], init = 0) -
            sum(q_sdd[ssd_id, i] for ssd_id in bus_sdd_consumer_ids[uid], init = 0) +
            sum(
                shunt_lookup[shunt_id]["bs"]*shunt_step[shunt_id, i]
                for shunt_id in bus_shunt_ids[uid],
                init = 0
            )*vm[uid, i]^2
            #bs*vm[uid]^2
        )
    end

    for (uid,ac_line) in ac_line_lookup, i in periods
        branch = ac_line
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key, i]
        q_fr = q_branch[fr_key, i]
        p_to = p_branch[to_key, i]
        q_to = q_branch[to_key, i]

        #if relax_thermal_limits
        #    s_thermal = ac_thermal_slack[uid]
        #else
        #    s_thermal = 0.0
        #end

        vm_fr = vm[branch["fr_bus"], i]
        vm_to = vm[branch["to_bus"], i]
        va_fr = va[branch["fr_bus"], i]
        va_to = va[branch["to_bus"], i]

        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        tm = 1.0
        ta = 0.0
        tr = tm * cos(ta)
        ti = tm * sin(ta)
        ttm = tm^2

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        # variable bounds
        JuMP.set_lower_bound(p_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(p_to, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_to, -10.0*branch["mva_ub_nom"])

        JuMP.set_upper_bound(p_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(p_to,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_to,  10.0*branch["mva_ub_nom"])

        # From side of the branch flow
        JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, vad_lb <= va_fr - va_to <= vad_ub)

        # Apparent power limit, from side and to side
        # TODO: What is the best way to formulate this as a soft constraint?
        # Square root? Or (s_max + slack)^2?
        # TODO: Potentially re-add option to relax
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= (branch["mva_ub_nom"])^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= (branch["mva_ub_nom"])^2)
    end

    for (uid,twt) in twt_lookup, i in periods
        branch = twt
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key, i]
        q_fr = q_branch[fr_key, i]
        p_to = p_branch[to_key, i]
        q_to = q_branch[to_key, i]

        #if relax_thermal_limits
        #    s_thermal = twt_thermal_slack[uid]
        #else
        #    s_thermal = 0.0
        #end

        vm_fr = vm[branch["fr_bus"], i]
        vm_to = vm[branch["to_bus"], i]
        va_fr = va[branch["fr_bus"], i]
        va_to = va[branch["to_bus"], i]

        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        initial_status = branch["initial_status"]

        tm = initial_status["tm"]
        ta = initial_status["ta"]
        tr = tm * cos(ta)
        ti = tm * sin(ta)
        ttm = tm^2

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        # variable bounds
        JuMP.set_lower_bound(p_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_fr, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(p_to, -10.0*branch["mva_ub_nom"])
        JuMP.set_lower_bound(q_to, -10.0*branch["mva_ub_nom"])

        JuMP.set_upper_bound(p_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_fr,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(p_to,  10.0*branch["mva_ub_nom"])
        JuMP.set_upper_bound(q_to,  10.0*branch["mva_ub_nom"])

        # From side of the branch flow
        JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, vad_lb <= va_fr - va_to <= vad_ub)

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= (branch["mva_ub_nom"])^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= (branch["mva_ub_nom"])^2)
    end

    for (uid,dc_line) in dc_line_lookup, i in periods
        branch = dc_line
        fr_key = (uid, branch["fr_bus"], branch["to_bus"])
        to_key = (uid, branch["to_bus"], branch["fr_bus"])

        p_fr = p_branch[fr_key, i]
        q_fr = q_branch[fr_key, i]
        p_to = p_branch[to_key, i]
        q_to = q_branch[to_key, i]

        dc_line["pminf"] = -dc_line["pdc_ub"]
        dc_line["pmaxf"] = dc_line["pdc_ub"]
        dc_line["pmint"] = -dc_line["pdc_ub"]
        dc_line["pmaxt"] = dc_line["pdc_ub"]

        dc_line["qminf"] = dc_line["qdc_fr_lb"]
        dc_line["qmaxf"] = dc_line["qdc_fr_ub"]
        dc_line["qmint"] = dc_line["qdc_to_lb"]
        dc_line["qmaxt"] = dc_line["qdc_to_ub"]

        # variable bounds
        JuMP.set_lower_bound(p_fr, -branch["pdc_ub"])
        JuMP.set_lower_bound(q_fr,  branch["qdc_fr_lb"])
        JuMP.set_lower_bound(p_to, -branch["pdc_ub"])
        JuMP.set_lower_bound(q_to,  branch["qdc_to_lb"])

        JuMP.set_upper_bound(p_fr, branch["pdc_ub"])
        JuMP.set_upper_bound(q_fr, branch["qdc_fr_ub"])
        JuMP.set_upper_bound(p_to, branch["pdc_ub"])
        JuMP.set_upper_bound(q_to, branch["qdc_to_ub"])

        # From side of the branch flow
        JuMP.@constraint(model, p_fr + p_to == 0.0)
    end

    # Add ramping constraints
    add_ramp_constraints!(model, sdd_lookup, periods, sdd_ids, dt;
        p=p_sdd, u_on=on_status, u_su=u_su, u_sd=u_sd,
    )

    return model
end


function extract_data_from_model(
    model::JuMP.Model, input_data::NamedTuple, on_status::Dict, p::Dict;
    tolerance = 1e-6,
    allow_switching = true,
)::Dict
    data = input_data

    # Filter out inactive SDDs. Model variables are only defined for active
    # SDDs, but we must know about all the SDDs to populate solution files.
    # Note that unlike the model function above, this function actually
    # relies on on_status and p being provided.
    filtered_data = _filter_inactive_sdds(
        data, on_status, p; tolerance = tolerance
    )
    active_sdd_set = Set(uid for uid in filtered_data.sdd_ids)

    # Dict mapping UIDs to dicts of output parameters
    opf_solution = Dict{String, Any}()

    model_has_duals = JuMP.has_duals(model)

    #
    # Add solution data for buses
    #
    bus_key = "bus"
    bus_map = Dict{String, Dict}()
    for uid in data.bus_ids
        bus_uid_map = Dict{String, Any}()
        bus_uid_map["vm"] = JuMP.value(model[:vm][uid])
        bus_uid_map["va"] = JuMP.value(model[:va][uid])
        if model_has_duals
            # Need names for these constraints
            bus_uid_map["p_balance_dual"] = JuMP.dual(model[:p_balance][uid])
            bus_uid_map["q_balance_dual"] = JuMP.dual(model[:q_balance][uid])
        end
        bus_map[uid] = bus_uid_map
    end
    opf_solution[bus_key] = bus_map
    ###

    #
    # Add solution data for shunts
    #
    shunt_key = "shunt"
    shunt_map = Dict{String, Any}()
    for uid in data.shunt_ids
        shunt_uid_map = Dict{String, Any}(
            #"step" => data.shunt_lookup[uid]["initial_status"]["step"],
            "step" => JuMP.value(model[:shunt_step][uid]),
        )
        shunt_map[uid] = shunt_uid_map
    end
    opf_solution[shunt_key] = shunt_map
    ###

    #
    # Add solution data for SDDs
    #
    sdd_key = "simple_dispatchable_device"
    sdd_map = Dict{String, Any}()
    for uid in data.sdd_ids
        sdd_uid_map = Dict{String, Any}(
            # Just use the on status that was provided
            # Probably don't even need to include on_status in this dict.
            "on_status" => on_status[uid],
            # Just use the real power that was provided
            # TODO: If this gets relaxed in the OPF model, send the
            # power calculated by OPF.
            #
            # This is incorrect. p_on is the power decision variable, which
            # is different from the power sent to the OPF model, which includes
            # SU/SD power. These SU/SD powers should not be included in solution
            # files, and do not correspond to the field p_on.
            #"p_on" => p[uid],
            # Use reactive power calculated by OPF, if active. Else zero.
        )
        if uid in active_sdd_set
            # active_sdd_set includes SDDs in SU/SD curves. However, these
            # should not be given a "p_on" value
            #if abs(on_status[uid] - 1.0) <= tolerance
            if round(on_status[uid]) == 1
                sdd_uid_map["p_on"] = JuMP.value(model[:p_sdd][uid])
            else
                sdd_uid_map["p_on"] = 0.0
            end
            # Even if the device is not on, it can have a reactive power
            # if it is an SU/SD curve.
            sdd_uid_map["q"] = JuMP.value(model[:q_sdd][uid])
        else
            sdd_uid_map["p_on"] = 0.0
            sdd_uid_map["q"] = 0.0
        end
        sdd_map[uid] = sdd_uid_map
    end
    opf_solution[sdd_key] = sdd_map
    ###

    #
    # Add solution data for AC lines
    #
    ac_key = "ac_line"
    ac_map = Dict{String, Any}()
    for uid in data.ac_line_ids
        if !allow_switching
            # If switching is not allowed, lines with initial status of "off"
            # were excluded from the formulation. In this case, use the initial
            # status as the on status.
            init_on = data.ac_line_lookup[uid]["initial_status"]["on_status"]
            ac_uid_map = Dict{String, Any}("on_status" => init_on)
        else
            init_on = data.ac_line_lookup[uid]["initial_status"]["on_status"]
            branch = input_data.ac_line_lookup[uid]
            fr_key = (uid, branch["fr_bus"], branch["to_bus"])
            to_key = (uid, branch["to_bus"], branch["fr_bus"])
            p_fr = model[:p_branch][fr_key]
            p_to = model[:p_branch][to_key]
            q_fr = model[:q_branch][fr_key]
            q_to = model[:q_branch][to_key]
            if (JuMP.value(p_fr) <= tolerance && JuMP.value(p_to) <= tolerance
                && JuMP.value(q_fr) <= tolerance && JuMP.value(q_to) <= tolerance)
                # If we are not using the line, just use the initial status. This
                # prevents turning off the line, which incurs a shutdown cost
                # and can lead to a disconnected graph.
                line_on_status = init_on
            else
                # If we are using the line, make sure it is on.
                # By allowing OPF to determine whether the line is on, we assume
                # that the benefit gained by more flexible OPF is greater than
                # the SU/SD cost incurred by potentially line switching.
                line_on_status = 1
            end
            ac_uid_map = Dict{String, Any}("on_status" => line_on_status)
        end
        ac_map[uid] = ac_uid_map
    end
    opf_solution[ac_key] = ac_map
    ###

    #
    # Add solution data for TWT
    #
    twt_key = "two_winding_transformer"
    twt_map = Dict{String, Any}()
    for uid in data.twt_ids
        if !allow_switching
            # If switching is not allowed, we use initial status for everything
            init_on = data.twt_lookup[uid]["initial_status"]["on_status"]
            twt_uid_map = Dict(
                "tm" => data.twt_lookup[uid]["initial_status"]["tm"],
                "ta" => data.twt_lookup[uid]["initial_status"]["ta"],
                "on_status" => init_on,
            )
        else
            init_on = data.twt_lookup[uid]["initial_status"]["on_status"]
            twt = input_data.twt_lookup[uid]
            fr_key = (uid, twt["fr_bus"], twt["to_bus"])
            to_key = (uid, twt["to_bus"], twt["fr_bus"])
            p_fr = model[:p_branch][fr_key]
            p_to = model[:p_branch][to_key]
            q_fr = model[:q_branch][fr_key]
            q_to = model[:q_branch][to_key]
            if (JuMP.value(p_fr) <= tolerance && JuMP.value(p_to) <= tolerance
                && JuMP.value(q_fr) <= tolerance && JuMP.value(q_to) <= tolerance)
                # If we are not using the line, just use the initial status. This
                # prevents turning off the line, which incurs a shutdown cost.
                line_on_status = init_on
            else
                # If we are using the line, make sure it is on.
                # By allowing OPF to determine whether the line is on, we assume
                # that the benefit gained by more flexible OPF is greater than
                # the SU/SD cost incurred by potentially line switching.
                line_on_status = 1
            end
            twt_uid_map = Dict(
                "tm" => data.twt_lookup[uid]["initial_status"]["tm"],
                "ta" => data.twt_lookup[uid]["initial_status"]["ta"],
                "on_status" => line_on_status,
            )
        end
        twt_map[uid] = twt_uid_map
    end
    opf_solution[twt_key] = twt_map

    #
    # Add solution data for DC lines
    #
    dc_key = "dc_line"
    dc_map = Dict{String, Any}()
    for uid in data.dc_line_ids
        fr_bus = data.dc_line_lookup[uid]["fr_bus"]
        to_bus = data.dc_line_lookup[uid]["to_bus"]
        dc_uid_map = Dict{String, Any}()
        dc_uid_map["pdc_fr"] = JuMP.value(model[:p_branch][(uid, fr_bus, to_bus)])
        dc_uid_map["qdc_fr"] = JuMP.value(model[:q_branch][(uid, fr_bus, to_bus)])
        dc_uid_map["qdc_to"] = JuMP.value(model[:q_branch][(uid, to_bus, fr_bus)])
        dc_map[uid] = dc_uid_map
    end
    opf_solution[dc_key] = dc_map

    return opf_solution
end

"""Extract data from solution to multiperiod model

on_status is required so we can distinguish between p_on and p_su/sd.
    It should map uids to an arrays of binary values.

"""
function extract_data_from_multiperiod_model(
    model, 
    data,
    on_status;
    p_sdd=nothing,
    q_sdd=nothing,
)
    solution_data = Dict{String, Any}()

    if p_sdd === nothing
        p_sdd = model[:p_sdd]
    end
    if q_sdd === nothing
        q_sdd = model[:q_sdd]
    end

    #
    # Add solution data for buses
    #
    bus_key = "bus"
    bus_map = Dict{String, Dict}()
    for uid in data.bus_ids
        bus_uid_map = Dict{String, Any}()
        bus_uid_map["vm"] = Vector(JuMP.value.(model[:vm][uid, :]))
        bus_uid_map["va"] = Vector(JuMP.value.(model[:va][uid, :]))
        # TODO: Extract duals in solution data if requested
        #if model_has_duals
        #    # Need names for these constraints
        #    bus_uid_map["p_balance_dual"] = JuMP.dual(model[:p_balance][uid])
        #    bus_uid_map["q_balance_dual"] = JuMP.dual(model[:q_balance][uid])
        #end
        bus_map[uid] = bus_uid_map
    end
    solution_data[bus_key] = bus_map

    #
    # Add solution data for shunts
    #
    shunt_key = "shunt"
    shunt_map = Dict{String, Any}()
    for uid in data.shunt_ids
        shunt_uid_map = Dict{String, Any}(
            #"step" => data.shunt_lookup[uid]["initial_status"]["step"],
            "step" => Vector(JuMP.value.(model[:shunt_step][uid, :])),
        )
        shunt_map[uid] = shunt_uid_map
    end
    solution_data[shunt_key] = shunt_map

    #
    # Add solution data for SDDs
    #
    sdd_key = "simple_dispatchable_device"
    sdd_map = Dict{String, Any}()
    for uid in data.sdd_ids
        sdd_uid_map = Dict{String, Any}(
            # Just use the on status that was provided
            # Probably don't even need to include on_status in this dict.
            "on_status" => on_status[uid],
            "p_on" => [0.0 for _ in data.periods],
            "q" => Vector(JuMP.value.(q_sdd[uid, :])),
        )
        for i in data.periods
            if Bool(on_status[uid][i])
                sdd_uid_map["p_on"][i] = JuMP.value(p_sdd[uid, i])
            end
        end
        sdd_map[uid] = sdd_uid_map
    end
    solution_data[sdd_key] = sdd_map

    #
    # Add solution data for AC lines
    #
    ac_key = "ac_line"
    ac_map = Dict{String, Any}()
    for uid in data.ac_line_ids
        # TODO: Potentially support !allow_switching?
        ac_map[uid] = Dict{String, Any}(
            "on_status" => [1 for _ in data.periods],
        )
    end
    solution_data[ac_key] = ac_map

    #
    # Add solution data for TWT
    #
    twt_key = "two_winding_transformer"
    twt_map = Dict{String, Any}()
    for uid in data.twt_ids
        # TODO: Potentially support !allow_switching
        twt_map[uid] = Dict(
            "tm" => [data.twt_lookup[uid]["initial_status"]["tm"] for _ in data.periods],
            "ta" => [data.twt_lookup[uid]["initial_status"]["ta"] for _ in data.periods],
            "on_status" => [1 for _ in data.periods],
        )
    end
    solution_data[twt_key] = twt_map

    #
    # Add solution data for DC lines
    #
    dc_key = "dc_line"
    dc_map = Dict{String, Any}()
    for uid in data.dc_line_ids
        fr_bus = data.dc_line_lookup[uid]["fr_bus"]
        to_bus = data.dc_line_lookup[uid]["to_bus"]
        dc_uid_map = Dict{String, Any}()
        dc_uid_map["pdc_fr"] = Vector(JuMP.value.(model[:p_branch][(uid, fr_bus, to_bus), :]))
        dc_uid_map["qdc_fr"] = Vector(JuMP.value.(model[:q_branch][(uid, fr_bus, to_bus), :]))
        dc_uid_map["qdc_to"] = Vector(JuMP.value.(model[:q_branch][(uid, to_bus, fr_bus), :]))
        dc_map[uid] = dc_uid_map
    end
    solution_data[dc_key] = dc_map

    return solution_data
end

function compute_line_capacity_ub(input_data, branch_id, frbus_id, tobus_id)
    acl_set = keys(input_data.ac_line_lookup)
    twt_set = keys(input_data.twt_lookup)
    dcl_set = keys(input_data.dc_line_lookup)

    vm_fr = input_data.bus_lookup[frbus_id]["vm_ub"]
    vm_to = input_data.bus_lookup[tobus_id]["vm_ub"]

    if branch_id in acl_set
        # For AC lines, tm=1 and ta=0, so ttm=1, tr=tm=1, and ti=0
        branch = input_data.ac_line_lookup[branch_id]
        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        p_fr_ub = (
            abs((g+g_fr)*vm_fr^2)
            + abs((-g)*vm_fr*vm_to) #*cos(va_fr-va_to)
            + abs((-b)*vm_fr*vm_to) #*sin(va_fr-va_to)
        )
        p_to_ub = (
            abs((g+g_to)*vm_to^2)
            + abs((-g)*vm_to*vm_fr) #*cos(va_to-va_fr)
            + abs((-b)*vm_to*vm_fr) #*sin(va_to-va_fr)
        )

        # Just take the upper bound of from and to. Really, we should
        # choose one depending on whether the bus has producers or
        # consumers.
        line_capacity_ub = max(p_fr_ub, p_to_ub)
    elseif branch_id in twt_set
        branch = input_data.twt_lookup[branch_id]

        r = branch["r"]
        x = branch["x"]
        y = 1 / (r + im*x)
        g = real(y)
        b = imag(y)

        # tm always participates in the denominator, so use its lower bound
        tm = branch["tm_lb"]
        # Assume cos/sin of ta give us 1.0
        tr = tm
        ti = tm
        ttm = tm^2

        g_fr = 0.0
        b_fr = branch["b"]/2.0
        g_to = 0.0
        b_to = branch["b"]/2.0

        if branch["additional_shunt"] == 1
            g_fr += branch["g_fr"]
            b_fr += branch["b_fr"]
            g_to += branch["g_to"]
            b_to += branch["b_to"]
        end

        p_fr_ub = (
            abs((g+g_fr)/ttm*vm_fr^2)
            + abs((-g*tr+b*ti)/ttm*vm_fr*vm_to) #*cos(va_fr-va_to)
            + abs((-b*tr-g*ti)/ttm*vm_fr*vm_to) #*sin(va_fr-va_to)
        )
        p_to_ub = (
            abs((g+g_to)*vm_to^2)
            + abs((-g*tr-b*ti)/ttm*vm_to*vm_fr) #*cos(va_to-va_fr)
            + abs((-b*tr+g*ti)/ttm*vm_to*vm_fr) #*sin(va_to-va_fr)
        )

        # Just take the upper bound of from and to. Really, we should
        # choose one depending on whether the bus has producers or
        # consumers.
        line_capacity_ub = max(p_fr_ub, p_to_ub)
    elseif branch_id in dcl_set
        line_capacity_ub = input_data.dc_line_lookup[branch_id]["pdc_ub"]
    end
    return line_capacity_ub
end

function get_line_ub_em(input_data, branch_id)
    acl_set = keys(input_data.ac_line_lookup)
    twt_set = keys(input_data.twt_lookup)
    dcl_set = keys(input_data.dc_line_lookup)

    if branch_id in acl_set
        return input_data.ac_line_lookup[branch_id]["mva_ub_em"]
    elseif branch_id in twt_set
        return input_data.twt_lookup[branch_id]["mva_ub_em"]
    elseif branch_id in dcl_set
        line_capacity_ub = input_data.dc_line_lookup[branch_id]["pdc_ub"]
    end
    return line_capacity_ub
end

function get_line_ub_nom(input_data, branch_id)
    acl_set = keys(input_data.ac_line_lookup)
    twt_set = keys(input_data.twt_lookup)
    dcl_set = keys(input_data.dc_line_lookup)

    if branch_id in acl_set
        return input_data.ac_line_lookup[branch_id]["mva_ub_nom"]
    elseif branch_id in twt_set
        return input_data.twt_lookup[branch_id]["mva_ub_nom"]
    elseif branch_id in dcl_set
        line_capacity_ub = input_data.dc_line_lookup[branch_id]["pdc_ub"]
    end
    return line_capacity_ub
end

function get_initial_on_status(input_data, branch_id)
    acl_set = keys(input_data.ac_line_lookup)
    twt_set = keys(input_data.twt_lookup)
    dcl_set = keys(input_data.dc_line_lookup)

    if branch_id in acl_set
        return input_data.ac_line_lookup[branch_id]["initial_status"]["on_status"]
    elseif branch_id in twt_set
        return input_data.twt_lookup[branch_id]["initial_status"]["on_status"]
    elseif branch_id in dcl_set
        return 1
    end
end

# Compute (a lower bound on) the cost incurred if each device is fixed to its
# scheduled value.
function compute_device_balance_costs(input_data, schedule_data; allow_switching = true)
    topo_data = preprocess_topology_data(input_data)

    sdd2bus = Dict()
    for bus_id in input_data.bus_ids
        for sdd_id in topo_data.bus_sdd_ids[bus_id]
            sdd2bus[sdd_id] = bus_id
        end
    end

    balance_cost = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    for sdd_id in input_data.sdd_ids
        # Array of scheduled power at each time point
        p_sched = schedule_data.real_power[sdd_id]

        bus_id = sdd2bus[sdd_id]
        consumers = [
            uid for uid in topo_data.bus_sdd_ids[bus_id]
            if input_data.sdd_lookup[uid]["device_type"] == "consumer"
        ]
        producers = [
            uid for uid in topo_data.bus_sdd_ids[bus_id]
            if input_data.sdd_lookup[uid]["device_type"] == "producer"
        ]
        has_both_prod_and_cons = (!isempty(consumers) && !isempty(producers))

        # Upper bound on the in/out capacity of all adjacent lines
        line_capacity_ub = 0.0
        for (branch_id, frbus_id, tobus_id) in topo_data.bus_branch_keys[bus_id]
            if allow_switching || get_initial_on_status(input_data, branch_id) == 1
                line_capacity_ub += min(
                    compute_line_capacity_ub(input_data, branch_id, frbus_id, tobus_id),
                    get_line_ub_nom(input_data, branch_id),
                )
            end
        end

        # If a bus only has producers or consumers, fixing a device to above
        # adjacent line flow limits will cause a power balance violation.
        if !has_both_prod_and_cons
            # Assume p is positive
            overcapacity = [max((p - line_capacity_ub), 0.0) for p in p_sched]
            p_bus_vio_cost = input_data.violation_cost["p_bus_vio_cost"]
            dt = input_data.dt
            balance_penalty = [dt[i]*p_bus_vio_cost*p for (i, p) in enumerate(overcapacity)]
            total_penalty = sum(balance_penalty; init = 0.0)
            balance_cost[sdd_id] = total_penalty
        end
    end
    return balance_cost
end
