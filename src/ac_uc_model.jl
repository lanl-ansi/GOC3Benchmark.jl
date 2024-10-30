import JuMP

function get_ac_uc_model(
    data::NamedTuple;
    optimizer=nothing,
    direct=false,
    relax_balances=true, # Use slack for AC power balances
    relax_p_balance=relax_balances,
    relax_q_balance=relax_balances,
    include_reserves=false,
    fix_shunt_steps=true,
    max_balance_violation=nothing,
)
    # Various other arguments, e.g. those to allow feasibility tests, may
    # be useful eventually.
    if optimizer === nothing
        if direct
            throw(ArgumentError("Direct model requires an optimizer"))
        end
        model = JuMP.Model()
    elseif direct
        model = JuMP.direct_model(optimizer)
    else
        model = JuMP.Model(optimizer)
    end

    # Data processing for unit commitment model
    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
    ) = data
    # Compute start, stop, and midpoint of each interval
    interval_times = Vector{Any}()
    start = 0.0
    stop = 0.0
    mid = 0.0
    for (i, duration) in enumerate(dt)
        start = stop
        stop += duration
        mid = 0.5*(start + stop)
        push!(interval_times, (start=start, stop=stop, mid=mid))
    end

    # Data processing for ACOPF model
    vad_ub = deg2rad(30)
    vad_lb = -vad_ub
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

    # Add scheduling model variables and constraints
    p, q = add_power_variables!(model, sdd_ts_lookup, periods, sdd_ids)
    p_on_status = add_p_on_status_variables!(model, sdd_ts_lookup, periods, sdd_ids)
    mustrun_init, outage_init = add_init_mustrun_outage_constraints!(
        model, sdd_lookup, periods, sdd_ids, interval_times
    )
    out = add_su_sd_variables_constraints!(
        model, sdd_lookup, periods, sdd_ids
    )
    u_su, u_sd, on_status_evolution, prohibit_su_sd = out

    min_up, min_dn = add_min_up_down_constraints!(
        model, sdd_lookup, periods, sdd_ids, interval_times
    )

    max_starts_over_intervals = add_maximum_starts_over_intervals!(
        model, sdd_lookup, periods, sdd_ids, interval_times
    )
    _, p_su_ru, _, p_sd_rd = get_contiguous_power_curves(data)
    T_su_pc, T_sd_pc = get_inverse_power_curve_intervals(data)
    p_su, p_sd, su_power_curve, sd_power_curve = add_startup_shutdown_power!(
        model, sdd_ts_lookup, periods, sdd_ids, p_su_ru, T_su_pc, p_sd_rd, T_sd_pc
    )
    p_on, p_on_su_sd_aggregation = add_p_on_su_sd_aggregation!(
        model, sdd_ts_lookup, sdd_ids, periods
    )
    ramp_ub, ramp_lb = add_ramp_constraints!(
        model, sdd_lookup, periods, sdd_ids, dt
    )
    max_en_over_intervals, min_en_over_intervals = add_max_min_energy_over_intervals!(
        model, sdd_lookup, periods, sdd_ids, interval_times, dt
    )

    # Add reserves, which are necessary for implication constraints if present
    azr_ids = data.azr_ids
    rzr_ids = data.rzr_ids
    if include_reserves
        reserve_vars = add_reserve_variables!(model, data)
        (p_rgu, p_rgd, p_scr, p_nsc, p_rru_on, p_rrd_on, p_rru_off,
         p_rrd_off, q_qru, q_qrd) = reserve_vars

        res_costs = _get_reserve_costs(data)
        (c_rgu, c_rgd, c_scr, c_nsc, c_rru_on, c_rrd_on, c_rru_off,
         c_rrd_off, c_qru, c_qrd) = res_costs

        res_shortfalls = add_reserve_shortfall_variables!(model, data)
        (p_rgu_slack, p_rgd_slack, p_scr_slack, p_nsc_slack, p_rru_slack,
         p_rrd_slack, q_qru_slack, q_qrd_slack) = res_shortfalls

        res_penalties = _get_reserve_shortfall_penalties(data)
        z_rgu, z_rgd, z_scr, z_nsc, z_rru, z_rrd, z_qru, z_qrd = res_penalties

        add_reserve_balances!(
            # Note we are relaxing reserves here. This will change if we support
            # solving reserve-feasibility problems.
            model, data, model[:p]; relax = true
        )

        # Reserve cost expressions
        @expression(model,
            c_rgu_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rgu[uid][i] * p_rgu[uid, i]
        )
        @expression(model,
            c_rgd_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rgd[uid][i] * p_rgd[uid, i]
        )
        @expression(model,
            c_scr_expr[uid in sdd_ids, i in periods],
            dt[i] * c_scr[uid][i] * p_scr[uid, i]
        )
        @expression(model,
            c_nsc_expr[uid in sdd_ids, i in periods],
            dt[i] * c_nsc[uid][i] * p_nsc[uid, i]
        )
        @expression(model,
            c_rru_on_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rru_on[uid][i] * p_rru_on[uid, i]
        )
        @expression(model,
            c_rrd_on_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rrd_on[uid][i] * p_rrd_on[uid, i]
        )
        @expression(model,
            c_rru_off_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rru_off[uid][i] * p_rru_off[uid, i]
        )
        @expression(model,
            c_rrd_off_expr[uid in sdd_ids, i in periods],
            dt[i] * c_rrd_off[uid][i] * p_rrd_off[uid, i]
        )
        @expression(model,
            c_qru_expr[uid in sdd_ids, i in periods],
            dt[i] * c_qru[uid][i] * q_qru[uid, i]
        )
        @expression(model,
            c_qrd_expr[uid in sdd_ids, i in periods],
            dt[i] * c_qrd[uid][i] * q_qrd[uid, i]
        )

        # Reserve shortfall penalty expressions
        # Active zonal reserves
        @expression(model,
            rgu_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_rgu[uid] * p_rgu_slack[uid, i]
        )
        @expression(model,
            rgd_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_rgd[uid] * p_rgd_slack[uid, i]
        )
        @expression(model,
            scr_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_scr[uid] * p_scr_slack[uid, i]
        )
        @expression(model,
            nsc_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_nsc[uid] * p_nsc_slack[uid, i]
        )
        @expression(model,
            rru_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_rru[uid] * p_rru_slack[uid, i]
        )
        @expression(model,
            rrd_sf_expr[uid in azr_ids, i in periods],
            dt[i] * z_rrd[uid] * p_rrd_slack[uid, i]
        )
        # Reactive zonal reserves
        @expression(model,
            qru_sf_expr[uid in rzr_ids, i in periods],
            dt[i] * z_qru[uid] * q_qru_slack[uid, i]
        )
        @expression(model,
            qrd_sf_expr[uid in rzr_ids, i in periods],
            dt[i] * z_qrd[uid] * q_qrd_slack[uid, i]
        )

    else
        # If reserves are not included, create dummy expressions to use in the
        # objective function without more nested branching.
        @expression(model, c_rgu_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_rgd_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_scr_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_nsc_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_rru_on_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_rrd_on_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_rru_off_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_rrd_off_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_qru_expr[uid in sdd_ids, i in periods], 0.0)
        @expression(model, c_qrd_expr[uid in sdd_ids, i in periods], 0.0)

        @expression(model, rgu_sf_expr[uid in azr_ids, i in periods], 0.0)
        @expression(model, rgd_sf_expr[uid in azr_ids, i in periods], 0.0)
        @expression(model, scr_sf_expr[uid in azr_ids, i in periods], 0.0)
        @expression(model, nsc_sf_expr[uid in azr_ids, i in periods], 0.0)
        @expression(model, rru_sf_expr[uid in azr_ids, i in periods], 0.0)
        @expression(model, rrd_sf_expr[uid in azr_ids, i in periods], 0.0)

        @expression(model, qru_sf_expr[uid in rzr_ids, i in periods], 0.0)
        @expression(model, qrd_sf_expr[uid in rzr_ids, i in periods], 0.0)
    end

    add_on_su_sd_implication_constraints!(
        model, data; include_reserves = include_reserves
    )
    add_reactive_power_implication_constraints!(
        model, data, T_su_pc, T_sd_pc; include_reserves = include_reserves
    )
    add_real_reactive_linking_constraints!(
        model, data, T_su_pc, T_sd_pc; include_reserves = include_reserves
    )

    on_status = model[:p_on_status]

    # TODO: Insert AC variables and equations here
    # Will need to audit all the constraints being added to make sure they still
    # make sense in this version of the model. Major differences are:
    # - This model has distinct p_on/p_su/p_sd variables
    # - on_status is a variable
    # - u_su and u_sd are variables
    # - We are including reserves
    #
    # Any constraints I will need to add that were not relevant before?
    @variable(model, bus_lookup[uid]["vm_lb"] <= vm[uid in bus_ids, t in periods] <= bus_lookup[uid]["vm_ub"], start=1.0)
    @variable(model, va[uid in bus_ids, t in periods])

    # Note that device p/q variables are already defined

    @variable(model, p_branch[uid in branch_keys, t in periods])
    @variable(model, q_branch[uid in branch_keys, t in periods])
    # TODO: Support thermal limit relaxation if necessary
    #if relax_thermal_limits
    #    @variable(model, ac_thermal_slack[uid in ac_line_ids, t in periods] >= 0.0)
    #    @variable(model, twt_thermal_slack[uid in twt_ids, t in periods] >= 0.0)
    #end

    # Cost blocks for objective
    cost_blocks, cost_block_p, p_disagg, device_cost = add_device_power_cost!(
        model, sdd_ts_lookup, periods, sdd_ids 
    )
    cost_by_period = @expression(
        model,
        cost_by_period[i=periods],
        sum(device_cost[uid, i] for uid in sdd_ids_consumer) -
        sum(device_cost[uid, i] for uid in sdd_ids_producer)
    )

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

    # Arbitrarily fix voltage angle for a reference bus
    for (uid, bus) in bus_lookup
        @constraint(model, [i in periods], va[uid, i] == 0.0)
        break
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

    on_cost, su_cost, sd_cost = add_on_su_sd_cost!(model, data)

    # TODO: Add penalty terms for AC slack variables
    # This may be considered "base operating cost". We will add to this
    # the reserve operating cost, and the penalties we incur for not satisfying
    # zonal reserve requirements and AC power balances.
    cost_expr = @expression(model,
        cost_expr,
        sum(dt[i]*cost_by_period[i] for i in periods)
        - on_cost - su_cost - sd_cost
    )
    reserve_cost_expr = @expression(model,
        reserve_cost_expr,
        # Reserve cost expressions
        # If reserves are not included, these are zero. Note that these have
        # already been multiplied by dt.
        - sum(c_rgu_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_rgd_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_scr_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_nsc_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_rru_on_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_rrd_on_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_rru_off_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_rrd_off_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_qru_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
        - sum(c_qrd_expr[uid, i] for uid in sdd_ids for i in periods, init = 0)
    )

    @objective(model,
        Max,
        # These are zero if we are only penalizing reserve shortfall
        cost_expr + reserve_cost_expr
        # Reserve shortfall expressions
        # If reserves are not included, these are zero. Note that these have
        # already been multiplied by dt.
        # Active
        - sum(rgu_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        - sum(rgd_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        - sum(scr_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        - sum(nsc_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        - sum(rru_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        - sum(rrd_sf_expr[uid, i] for uid in azr_ids for i in periods, init = 0)
        # Reactive
        - sum(qru_sf_expr[uid, i] for uid in rzr_ids for i in periods, init = 0)
        - sum(qrd_sf_expr[uid, i] for uid in rzr_ids for i in periods, init = 0)
        - sum(
            dt[t] * (p_balance_penalty_term[t] + q_balance_penalty_term[t])
            for t in periods
        )
    )

    # Define aliases for p and q using the variable names we used in the OPF model
    p_sdd = p
    q_sdd = q

    if relax_p_balance
        @constraint(model,
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
        @constraint(model,
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
        @constraint(model,
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
        @constraint(model, 
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
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

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
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

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

    return model
end
