using JSON

include("reserves.jl")

function read_solution_data_from_file(fname::String; process = true)
    data = Dict{String, Any}()
    open(fname, "r") do io
        data = JSON.parse(io)
    end
    if process
        data = process_solution_data(data)
    end
    return data
end

function process_solution_data(data::Dict)
    sdd_key = "simple_dispatchable_device"
    sdds = data["time_series_output"][sdd_key]
    sdd_ts_lookup = Dict(sdd["uid"] => sdd for sdd in sdds)
    return (
        sdd_ts_lookup = sdd_ts_lookup,
    )
end

"""
Parameters:
- input_data
- schedule_data
- opf_data
- write_duals
- include_reserves: if opf_data is provided, we compute reserves via LP.
  otherwise, we use reserve_data
- reserve_data: if opf_data is not provided and reserves are requested,
  these are the reserves we will use. Presumably these come from the scheduling
  model.
"""
function construct_solution_dict(
    input_data::NamedTuple,
    schedule_data::NamedTuple;
    opf_data=nothing,
    write_duals=false,
    include_reserves=true,
    reserve_data=nothing,
    postprocess=true,
    print_projected_devices=true,
)::Dict
    solfile_dict = Dict{String, Any}()
    solfile_dict["time_series_output"] = Dict{String, Any}()
    tsout = solfile_dict["time_series_output"]

    sdd_key = "simple_dispatchable_device"
    ac_key = "ac_line"
    dc_key = "dc_line"
    twt_key = "two_winding_transformer"
    shunt_key = "shunt_key"
    bus_key = "bus"

    has_opf_data = !(opf_data === nothing)

    sdd_ids = input_data.sdd_ids
    periods = input_data.periods
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup

    # Populate attributes for simple dispatchable devices
    sdd_key = "simple_dispatchable_device"
    devices = Vector()
    for uid in sdd_ids
        device = Dict()
        device["uid"] = uid

        if has_opf_data
            device["p_on"] = [
                opf_data[i][sdd_key][uid]["p_on"] for i in periods
            ]
        else
            # Only include real power if on status is 1
            device["p_on"] = [
                ifelse(
                    round(schedule_data.on_status[uid][i]) == 1.0,
                    schedule_data.real_power[uid][i],
                    0.0,
                )
                for i in periods
            ]
        end

        # Using q from ACOPF for the solution:
        if has_opf_data
            device["q"] = [opf_data[i][sdd_key][uid]["q"] for i in periods]
        else
            device["q"] = Vector(schedule_data.reactive_power[uid])
        end

        device["on_status"] = Int64.(round.(
            [schedule_data.on_status[uid][i] for i in periods]
        ))  
        # Set reserves to zero for now
        device["p_reg_res_up"] = [0.0 for i in periods]
        device["p_reg_res_down"] = [0.0 for i in periods]
        device["p_syn_res"] = [0.0 for i in periods]
        device["p_nsyn_res"] = [0.0 for i in periods]
        device["p_ramp_res_up_online"] = [0.0 for i in periods]
        device["p_ramp_res_down_online"] = [0.0 for i in periods]
        device["p_ramp_res_up_offline"] = [0.0 for i in periods]
        device["p_ramp_res_down_offline"] = [0.0 for i in periods]
        device["q_res_up"] = [0.0 for i in periods]
        device["q_res_down"] = [0.0 for i in periods]
        push!(devices, device)
    end 
    tsout[sdd_key] = devices

    # AC lines
    ac_key = "ac_line"
    ac_lines = Vector()
    for ac_id in input_data.ac_line_ids
        line = Dict{String, Any}()
        line["uid"] = ac_id
        if has_opf_data
            line["on_status"] = [
                opf_data[i][ac_key][ac_id]["on_status"] for i in periods
            ]
        else
            # Using initial_status would be better here as it ensures non-AC
            # solutions don't violate an allow_switching=false condition.
            line["on_status"] = [1 for i in periods]
        end
        push!(ac_lines, line)
    end
    tsout[ac_key] = ac_lines

    # DC lines
    dc_key = "dc_line"
    dc_lines = Vector()
    for dc_id in input_data.dc_line_ids
        line = Dict{String, Any}()
        line["uid"] = dc_id
        if has_opf_data
            line["pdc_fr"] = [opf_data[i][dc_key][dc_id]["pdc_fr"] for i in periods]
            line["qdc_fr"] = [opf_data[i][dc_key][dc_id]["qdc_fr"] for i in periods]
            line["qdc_to"] = [opf_data[i][dc_key][dc_id]["qdc_to"] for i in periods]
        else
            line["pdc_fr"] = [0.0 for i in periods]
            line["qdc_fr"] = [0.0 for i in periods]
            line["qdc_to"] = [0.0 for i in periods]
        end
        push!(dc_lines, line)
    end
    tsout[dc_key] = dc_lines

    # Two winding transformers
    twt_key = "two_winding_transformer"
    transformers = Vector()
    for twt_id in input_data.twt_ids
        twt = Dict{String, Any}()
        twt_params = input_data.twt_lookup[twt_id]
        twt["uid"] = twt_id
        if has_opf_data
            twt["tm"] = [opf_data[i][twt_key][twt_id]["tm"] for i in periods]
            twt["ta"] = [opf_data[i][twt_key][twt_id]["ta"] for i in periods]
            twt["on_status"] = [
                opf_data[i][twt_key][twt_id]["on_status"] for i in periods
            ]
        else
            # Could use the initial status instead here
            twt["tm"] = [0.5*(twt_params["tm_lb"] + twt_params["tm_ub"]) for _ in periods]
            twt["ta"] = [0.5*(twt_params["ta_lb"] + twt_params["ta_ub"]) for _ in periods]
            twt["on_status"] = [1 for i in periods]
        end
        push!(transformers, twt)
    end
    tsout[twt_key] = transformers

    # Shunts
    shunt_key = "shunt"
    shunts = Vector()
    for shunt_id in input_data.shunt_ids
        shunt = Dict{String, Any}()
        shunt["uid"] = shunt_id
        shunt_params = input_data.shunt_lookup[shunt_id]
        if has_opf_data
            steps = [opf_data[i][shunt_key][shunt_id]["step"] for i in periods]
            shunt["step"] = [
                # steps are currently floats (even if they should take integer
                # values). Attempt to make them ints so we don't get errors
                # parsing the solution file.
                # Note that this array will have all float values unless all
                # entries can be converted to integers.
                #
                # steps will potentially have to be rounded (and converted to
                # Int64) in a postprocessing function.
                #ifelse(isinteger(steps[i]), Int64(steps[i]), steps[i])
                Int64(round(steps[i]))
                for i in periods
            ]
        else
            shunt["step"] = [shunt_params["initial_status"]["step"] for _ in periods]
        end
        push!(shunts, shunt)
    end
    tsout[shunt_key] = shunts

    # Buses
    bus_key = "bus"
    buses = Vector()
    for bus_id in input_data.bus_ids
        bus = Dict{String, Any}()
        bus["uid"] = bus_id
        bus_params = input_data.bus_lookup[bus_id]
        if has_opf_data
            bus["vm"] = [opf_data[i][bus_key][bus_id]["vm"] for i in periods]
            bus["va"] = [opf_data[i][bus_key][bus_id]["va"] for i in periods]
        else
            # Note that we could also use initial_status here
            bus["vm"] = [0.5*(bus_params["vm_lb"] + bus_params["vm_ub"]) for _ in periods]
            bus["va"] = [0.0 for _ in periods]
        end
        push!(buses, bus)
    end
    tsout[bus_key] = buses

    # Only write duals if specified. Note that writing duals will prevent
    # evaluating the solution with check_data.py.
    if write_duals
        solfile_dict["time_series_duals"] = Dict{String, Any}()
        solfile_dict["time_series_duals"][bus_key] = Vector{Any}()
        for uid in input_data.bus_ids
            bus_params = input_data.bus_lookup[uid]
            bus_dual_dict = Dict{String, Any}()
            bus_dual_dict["uid"] = uid
            if (
                # Use first interval as proxy to check whether duals have been
                # exported.
                "p_balance_dual" in keys(opf_data[1][bus_key][uid])
                && "q_balance_dual" in keys(opf_data[1][bus_key][uid])
            )
                bus_dual_dict["p_balance_dual"] = [
                    opf_data[i][bus_key][uid]["p_balance_dual"]
                    for i in periods
                ]
                bus_dual_dict["q_balance_dual"] = [
                    opf_data[i][bus_key][uid]["q_balance_dual"]
                    for i in periods
                ]
            else
                # We were requested to write duals, but the OPF data does
                # not include duals. Contact Robby/Carleton if this occurs.
                throw(Exception)
            end
            push!(solfile_dict["time_series_duals"][bus_key], bus_dual_dict)
        end
    end

    #
    # Postprocess solution data
    #
    if postprocess
        solfile_dict, projected_device_intervals = postprocess_solution_data(
            # Note that we are sending this function the "original" processed
            # data rather than the tightened data. For the purpose of assigning
            # reserves, we do not want the tightened bounds.
            solfile_dict, input_data, schedule_data
        )
        if print_projected_devices
            println("Projected devices and intervals:")
            for (uid, i, diff) in projected_device_intervals
                println("$(uid), $(i), $(diff)")
            end
        end
    end

    if include_reserves
        # Solve LP to assign reserves
        if has_opf_data
            # Note that we could perform this solve without OPF data, but the
            # reserve calculation function is not set up for it. All we need to
            # compute reserves is p_on and q.
            #
            # Note that this calculation does not account for any projection
            # that may have been performed.
            #reserves = calculate_reserves_from_generation(
            #    input_data, opf_data, schedule_data.on_status
            #)
            #
            # This calculation, which uses the solution dict we build construct
            # in this function, does account for projection we just performed.
            reserves = calculate_reserves_from_generation(input_data, solfile_dict)
        elseif !(reserve_data === nothing)
            # If we don't have OPF data, we might be able to assign reserves
            # from scheduling model.
            reserves = reserve_data
        else
            # neither opf_data nor reserve_data were provided.
            # TODO: API to compute reserves via LP from only schedule data
            throw(Exception)
        end
        sdds_out = solfile_dict["time_series_output"]["simple_dispatchable_device"]
        for sdd in sdds_out
            uid = sdd["uid"]
            # TODO: Is it possible that, by assigning these reserves, we violate
            # p_max/p_min?
            sdd["p_reg_res_up"] = reserves.p_rgu[uid]
            sdd["p_reg_res_down"] = reserves.p_rgd[uid]
            sdd["p_syn_res"] = reserves.p_scr[uid]
            sdd["p_nsyn_res"] = reserves.p_nsc[uid]
            sdd["p_ramp_res_up_online"] = reserves.p_rru_on[uid]
            sdd["p_ramp_res_down_online"] = reserves.p_rrd_on[uid]
            sdd["p_ramp_res_up_offline"] = reserves.p_rru_off[uid]
            sdd["p_ramp_res_down_offline"] = reserves.p_rrd_off[uid]
            sdd["q_res_up"] = reserves.q_qru[uid]
            sdd["q_res_down"] = reserves.q_qrd[uid]
        end
    end

    return solfile_dict
end


function postprocess_solution_data(solution_data, input_data, schedule_data)
    # TODO: This function currently projects active powers using sequential
    # ramp bounds, rounds shunt steps, and assigns reserves. These should each
    # be moved into their own functions.
    #
    # Should I copy solution data here? May be useful for debugging.
    # solution_data = deepcopy(solution_data)
    tsout = solution_data["time_series_output"]
    sdd_output_array = tsout["simple_dispatchable_device"]
    shunt_output_array = tsout["shunt"]
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup
    #
    # Make sure the computed SDD real powers don't violate ramping constraints
    #
    prev_p_sdd = Dict{String, Float64}(
        uid => input_data.sdd_lookup[uid]["initial_status"]["p"]
        for uid in input_data.sdd_ids
    )
    prev_u_sdd = Dict{String, Int64}(
        uid => input_data.sdd_lookup[uid]["initial_status"]["on_status"]
        for uid in input_data.sdd_ids
    )

    # This will contain tuples of (uid, interval) indicating the devices for
    # which real powers were altered by the projection routine.
    projected_device_intervals = Vector{Tuple}()

    for i in input_data.periods
        # Ramping logic:
        # - if on and not SU, use p_ramp_up_ub
        # - if on and SU, use p_startup_ramp_ub
        # - if on (=> not SD), use p_ramp_down_ub
        # - if not on, use p_shutdown_ramp_ub
        dt = input_data.dt[i]
        for output in sdd_output_array
            uid = output["uid"]
            u_on = output["on_status"][i]

            p_on = output["p_on"][i]
            if Bool(u_on)
                # If device is on, use the power computed by OPF, which is in
                # the solution dict.
                p = output["p_on"][i]
            else
                # If device is off, use the scheduled power, which may include
                # SU/SD power.
                p = schedule_data.real_power[uid][i]
            end
            p_lb = sdd_ts_lookup[uid]["p_lb"][i]
            p_ub = sdd_ts_lookup[uid]["p_ub"][i]

            if Bool(u_on) && !Bool(prev_u_sdd[uid])
                # If we are on, and were not on previously, we are
                # in startup. p_startup_ramp_ub applies
                ramp_ub = prev_p_sdd[uid] + dt*sdd_lookup[uid]["p_startup_ramp_ub"]
                ramp_lb = prev_p_sdd[uid] - dt*sdd_lookup[uid]["p_ramp_down_ub"]
            elseif Bool(u_on)
                ramp_ub = prev_p_sdd[uid] + dt*sdd_lookup[uid]["p_ramp_up_ub"]
                ramp_lb = prev_p_sdd[uid] - dt*sdd_lookup[uid]["p_ramp_down_ub"]
            elseif !Bool(u_on)
                ramp_ub = prev_p_sdd[uid] + dt*sdd_lookup[uid]["p_startup_ramp_ub"]
                ramp_lb = prev_p_sdd[uid] - dt*sdd_lookup[uid]["p_shutdown_ramp_ub"]
            end

            # Look forward one step to decide if I need to further
            # constraint p_on due to an impending shutdown.
            if (
                i < length(input_data.periods)
                && Bool(u_on)
                && !Bool(output["on_status"][i+1])
            )
                # This is the first value of the SD power curve
                p_next = schedule_data.real_power[uid][i+1]
                dt_next = input_data.dt[i+1]
                sd_ramp_ub = p_next + dt_next*sdd_lookup[uid]["p_shutdown_ramp_ub"]

                ramp_ub = min(ramp_ub, sd_ramp_ub)
            end

            if Bool(u_on)
                # Only compute a new power if we have on_status==1. Otherwise,
                # power is fixed, even if we are in SU/SD.

                # Take the more conservative of the two bounds (ramping and
                # pre-specified)
                lb = max(p_lb, ramp_lb)
                ub = min(p_ub, ramp_ub)
                if ub < lb - 1e-8
                    # Print some info for debugging before raising error
                    println("$uid, $i")
                    println(schedule_data.real_power[uid])
                    println(schedule_data.on_status[uid])
                    lbs = [sdd_ts_lookup[uid]["p_lb"][i] for i in input_data.periods]
                    ubs = [sdd_ts_lookup[uid]["p_ub"][i] for i in input_data.periods]
                    println("lbs: ", lbs)
                    println("ubs: ", ubs)
                    println("ramp_lb: ", ramp_lb)
                    println("ramp_ub: ", ramp_ub)
                    println("(lb, ub): ", (lb, ub))
                end
                # skip assert so that solution will run through the evaluator
                #@assert ub >= lb - 1e-8

                if p > ub
                    p = ub
                end
                if p < lb
                    p = lb
                end

                if p != output["p_on"][i]
                    diff = p - output["p_on"][i]
                    output["p_on"][i] = p
                    push!(projected_device_intervals, (uid, i, diff))
                end
            end

            # Update previous status
            prev_p_sdd[uid] = p
            prev_u_sdd[uid] = output["on_status"][i]
        end
    end

    # Round all shunt steps to integer values.
    # Need to make sure the "step" field has type Int64 so that it is displayed
    # as an int in the JSON file.
    shunt_output_array = tsout["shunt"]
    old_shunt_steps = [shunt["step"] for shunt in tsout["shunt"]]
    tsout["shunt"] = [
        Dict("uid" => shunt["uid"], "step" => Int64.(round.(shunt["step"])))
        for shunt in tsout["shunt"]
    ]
    println("Rounded shunt steps:")
    for (old_val, shunt) in zip(old_shunt_steps, tsout["shunt"])
        if old_val != shunt["step"]
            uid = shunt["uid"]
            println("  $uid: $old_val -> $(shunt["step"])")
        end
    end

    #
    # Populate reserve values for SDDs
    #
    compute_reserves = false
    if compute_reserves
        periods = input_data.periods
        producer_set = Set(input_data.sdd_ids_producer)
        consumer_set = Set(input_data.sdd_ids_consumer)

        # This is the tolerance used to determine if a device is in SU/SD when
        # on_status=0. 1e-6 is Gurobi's default primal feasibility tolerance.
        # This is unreliable compared to a using a structural "power curve
        # intervals" that we can compute from on_status time-series.
        # TODO: Don't check a float to determine SU/SD.
        sched_tolerance = 1e-6
        for sdd in sdd_output_array
            uid = sdd["uid"]
            p_on = sdd["p_on"]
            q = sdd["q"]
            p_sched = schedule_data.real_power[uid]

            # Constants necessary to compute reactive power reserves.
            q_bound_cap = input_data.sdd_lookup[uid]["q_bound_cap"]
            q_linear_cap = input_data.sdd_lookup[uid]["q_linear_cap"]

            # Initialize reserves to zero
            p_reg_res_up = [0.0 for i in periods]
            p_reg_res_down = [0.0 for i in periods]
            p_syn_res = [0.0 for i in periods]
            p_nsyn_res = [0.0 for i in periods]
            p_ramp_res_up_online = [0.0 for i in periods]
            p_ramp_res_down_online = [0.0 for i in periods]
            p_ramp_res_up_offline = [0.0 for i in periods]
            p_ramp_res_down_offline = [0.0 for i in periods]
            q_res_up = [0.0 for i in periods]
            q_res_down = [0.0 for i in periods]

            # Compute reserves using generator slacks
            prgu_ub = input_data.sdd_lookup[uid]["p_reg_res_up_ub"]
            pscr_ub = input_data.sdd_lookup[uid]["p_syn_res_ub"]
            prru_on_ub = input_data.sdd_lookup[uid]["p_ramp_res_up_online_ub"]
            prgd_ub = input_data.sdd_lookup[uid]["p_reg_res_down_ub"]
            prrd_on_ub = input_data.sdd_lookup[uid]["p_ramp_res_down_online_ub"]
            pnsc_ub = input_data.sdd_lookup[uid]["p_nsyn_res_ub"]
            prru_off_ub = input_data.sdd_lookup[uid]["p_ramp_res_up_offline_ub"]
            prrd_off_ub = input_data.sdd_lookup[uid]["p_ramp_res_down_offline_ub"]
            for i in periods
                if uid in producer_set
                    #
                    # Branch on on_status to compute real power reserves
                    #
                    if sdd["on_status"][i] == 1
                        #
                        # Attempt to convert pmax_slack into reserve
                        #
                        pmax_slack = input_data.sdd_ts_lookup[uid]["p_ub"][i] - p_on[i]

                        # First try to use ub slack to fill p_rgu
                        p_rgu_ub_slack = max(prgu_ub - p_reg_res_up[i], 0.0)
                        # We need all these slacks, because prgu must respect all three
                        # of these bounds.
                        pscr_ub_slack = max(pscr_ub - p_reg_res_up[i] - p_syn_res[i], 0.0)
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_reg_res_up[i] += min(pmax_slack, p_rgu_ub_slack, pscr_ub_slack, prru_on_ub_slack)
                        pmax_slack -= min(pmax_slack, p_rgu_ub_slack, pscr_ub_slack, prru_on_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        # Make sure slack doesn't become negative
                        pmax_slack = max(pmax_slack, 0.0)

                        # Then try to use ub slack to fill p_scr
                        pscr_ub_slack = max(pscr_ub - p_reg_res_up[i] - p_syn_res[i], 0.0)
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_syn_res[i] += min(pmax_slack, pscr_ub_slack, prru_on_ub_slack)
                        pmax_slack -= min(pmax_slack, pscr_ub_slack, prru_on_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        # Make sure slack doesn't become negative
                        pmax_slack = max(pmax_slack, 0.0)

                        # Then try to use ub slack to fill p_rru_on
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_ramp_res_up_online[i] += min(pmax_slack, prru_on_ub_slack)
                        pmax_slack -= min(pmax_slack, prru_on_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        pmax_slack = max(pmax_slack, 0.0)

                        #
                        # Attempt to convert pmin_slack into reserve
                        #
                        pmin_slack = p_on[i] - input_data.sdd_ts_lookup[uid]["p_lb"][i]

                        # First try to use lb slack to fill p_rgd
                        prgd_slack = max(prgd_ub - p_reg_res_down[i], 0.0)
                        prrd_on_ub_slack = max(
                            (prrd_on_ub - p_reg_res_down[i]
                             - p_ramp_res_down_online[i]),
                            0.0,
                        )
                        p_reg_res_down[i] += min(pmin_slack, prgd_slack, prrd_on_ub_slack)
                        pmin_slack -= p_reg_res_down[i]
                        #@assert pmin_slack >= -1e-8
                        pmin_slack = max(pmin_slack, 0.0)

                        # Then try to use lb slack to fill p_rrd_on
                        prrd_on_ub_slack = max(
                            (prrd_on_ub - p_reg_res_down[i]
                             - p_ramp_res_down_online[i]),
                            0.0,
                        )
                        p_ramp_res_down_online[i] += min(pmin_slack, prrd_on_ub_slack)
                        pmin_slack -= p_ramp_res_down_online[i]
                        #@assert pmin_slack >= -1e-8
                        pmin_slack = max(pmin_slack, 0.0)
                    else
                        # We are not online. Attempt to set offline reserves
                        # These are a little more complicated because I need to know
                        # p_su and p_sd

                        #
                        # Attempt to assign pmax_slack to nsc and rru_off reserve
                        #
                        pmax_slack = input_data.sdd_ts_lookup[uid]["p_ub"][i] - p_sched[i]

                        pnsc_ub_slack = max(pnsc_ub - p_nsyn_res[i], 0.0)
                        # Note: these slacks should always be positive. They were negative
                        # before only because I had a violated assumption.
                        prru_off_ub_slack = max(
                            prru_off_ub - p_nsyn_res[i] - p_ramp_res_up_offline[i],
                            0.0,
                        )
                        p_nsyn_res[i] += min(pmax_slack, pnsc_ub_slack, prru_off_ub_slack)
                        pmax_slack -= min(pmax_slack, pnsc_ub_slack, prru_off_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        pmax_slack = max(pmax_slack, 0.0)

                        prru_off_ub_slack = max(
                            prru_off_ub - p_nsyn_res[i] - p_ramp_res_up_offline[i],
                            0.0,
                        )
                        p_ramp_res_up_offline[i] += min(pmax_slack, prru_off_ub_slack)
                        pmax_slack -= min(pmax_slack, prru_off_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        pmax_slack = max(pmax_slack, 0.0)
                    end

                    #
                    # Branch on p_sched to set reactive power reserves
                    #
                    if p_sched[i] >= sched_tolerance || sdd["on_status"][i] == 1
                        q_max_slack = input_data.sdd_ts_lookup[uid]["q_ub"][i] - q[i]
                        q_min_slack = q[i] - input_data.sdd_ts_lookup[uid]["q_lb"][i]

                        # This is the flag that tells us if the SDD implements
                        # *both* leq and geq inequalities linking p and q.
                        if q_bound_cap == 1
                            q_0_lb = input_data.sdd_lookup[uid]["q_0_lb"]
                            q_0_ub = input_data.sdd_lookup[uid]["q_0_ub"]
                            beta_lb = input_data.sdd_lookup[uid]["beta_lb"]
                            beta_ub = input_data.sdd_lookup[uid]["beta_ub"]

                            pq_max_slack = q_0_ub + beta_ub*p_sched[i] - q[i]
                            pq_min_slack = q[i] - q_0_lb + beta_lb*p_sched[i]

                            q_res_up[i] = min(q_max_slack, pq_max_slack)
                            q_res_down[i] = min(q_min_slack, pq_min_slack)
                        elseif q_linear_cap != 1
                            q_res_up[i] = q_max_slack
                            q_res_down[i] = q_min_slack
                        end
                    end
                elseif uid in consumer_set
                    if sdd["on_status"][i] == 1
                        #
                        # Attempt to convert pmax_slack into rgd and rrd_on
                        #
                        pmax_slack = input_data.sdd_ts_lookup[uid]["p_ub"][i] - p_on[i]

                        # First try to use ub slack to fill p_rgd
                        prgd_ub_slack = max(prgd_ub - p_reg_res_down[i], 0.0)
                        prrd_on_ub_slack = max(
                            (prrd_on_ub - p_reg_res_down[i]
                             - p_ramp_res_down_online[i]),
                            0.0,
                        )
                        p_reg_res_down[i] += min(pmax_slack, prgd_ub_slack, prrd_on_ub_slack)
                        pmax_slack -= min(pmax_slack, prgd_ub_slack, prrd_on_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        # Make sure slack doesn't become negative
                        pmax_slack = max(pmax_slack, 0.0)

                        # Then use remaining slack to fill p_rrd
                        prrd_on_ub_slack = max(
                            (prrd_on_ub - p_reg_res_down[i]
                             - p_ramp_res_down_online[i]),
                            0.0,
                        )
                        p_ramp_res_down_online[i] += min(pmax_slack, prrd_on_ub_slack)
                        pmax_slack -= p_ramp_res_down_online[i]
                        #@assert pmax_slack >= -1e-8
                        pmax_slack = max(pmax_slack, 0.0)

                        #
                        # Attempt to convert pmin_slack into rgu, scr, and rru_on
                        #
                        pmin_slack = p_on[i] - input_data.sdd_ts_lookup[uid]["p_lb"][i]

                        # First try to use ub slack to fill p_rgu
                        p_rgu_ub_slack = max(prgu_ub - p_reg_res_up[i], 0.0)
                        # We need all these slacks, because prgu must respect all three
                        # of these bounds.
                        pscr_ub_slack = max(pscr_ub - p_reg_res_up[i] - p_syn_res[i], 0.0)
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_reg_res_up[i] += min(pmin_slack, p_rgu_ub_slack, pscr_ub_slack, prru_on_ub_slack)
                        pmin_slack -= min(pmin_slack, p_rgu_ub_slack, pscr_ub_slack, prru_on_ub_slack)
                        #@assert pmin_slack >= -1e-8
                        # Make sure slack doesn't become negative
                        pmin_slack = max(pmin_slack, 0.0)

                        # Then try to use ub slack to fill p_scr
                        pscr_ub_slack = max(pscr_ub - p_reg_res_up[i] - p_syn_res[i], 0.0)
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_syn_res[i] += min(pmin_slack, pscr_ub_slack, prru_on_ub_slack)
                        pmin_slack -= min(pmin_slack, pscr_ub_slack, prru_on_ub_slack)
                        #@assert pmin_slack >= -1e-8
                        # Make sure slack doesn't become negative
                        pmin_slack = max(pmin_slack, 0.0)

                        # Then try to use ub slack to fill p_rru_on
                        prru_on_ub_slack = max(
                            (
                                prru_on_ub - p_ramp_res_up_online[i]
                                - p_syn_res[i] - p_reg_res_up[i]
                            ),
                            0.0,
                        )
                        p_ramp_res_up_online[i] += min(pmin_slack, prru_on_ub_slack)
                        pmin_slack -= min(pmin_slack, prru_on_ub_slack)
                        #@assert pmin_slack >= -1e-8
                        pmin_slack = max(pmin_slack, 0.0)

                    else
                        #
                        # Attempt to fill rrd_off with pmax_slack
                        #
                        pmax_slack = input_data.sdd_ts_lookup[uid]["p_ub"][i] - p_on[i]

                        prrd_off_ub_slack = max(prrd_off_ub - p_ramp_res_down_offline[i], 0.0)
                        p_ramp_res_down_offline[i] += min(pmax_slack, prrd_off_ub_slack)
                        pmax_slack -= min(pmax_slack, prrd_off_ub_slack)
                        #@assert pmax_slack >= -1e-8
                        pmax_slack = max(pmax_slack, 0.0)
                    end

                    #
                    # Branch on p_sched to set reactive power reserves
                    #
                    if p_sched[i] >= sched_tolerance || sdd["on_status"][i] == 1
                        q_max_slack = input_data.sdd_ts_lookup[uid]["q_ub"][i] - q[i]
                        q_min_slack = q[i] - input_data.sdd_ts_lookup[uid]["q_lb"][i]

                        # This is the flag that tells us if the SDD implements
                        # *both* leq and geq inequalities linking p and q.
                        if q_bound_cap == 1
                            q_0_lb = input_data.sdd_lookup[uid]["q_0_lb"]
                            q_0_ub = input_data.sdd_lookup[uid]["q_0_ub"]
                            beta_lb = input_data.sdd_lookup[uid]["beta_lb"]
                            beta_ub = input_data.sdd_lookup[uid]["beta_ub"]

                            pq_max_slack = q_0_ub + beta_ub*p_sched[i] - q[i]
                            pq_min_slack = q[i] - q_0_lb + beta_lb*p_sched[i]

                            q_res_down[i] = min(q_max_slack, pq_max_slack)
                            q_res_up[i] = min(q_min_slack, pq_min_slack)
                        elseif q_linear_cap != 1
                            # If q_linear_cap == 1, we cannot provide reactive
                            # reserve.
                            q_res_down[i] = q_max_slack
                            q_res_up[i] = q_min_slack
                        end
                    end
                else
                    throw(Exception)
                end
            end

            reserve_key_arrays = [
                ("p_reg_res_up", p_reg_res_up),
                ("p_reg_res_down", p_reg_res_down),
                ("p_syn_res", p_syn_res),
                ("p_nsyn_res", p_nsyn_res),
                ("p_ramp_res_up_online", p_ramp_res_up_online),
                ("p_ramp_res_down_online", p_ramp_res_down_online),
                ("p_ramp_res_up_offline", p_ramp_res_up_offline),
                ("p_ramp_res_down_offline", p_ramp_res_down_offline),
                ("q_res_up", q_res_up),
                ("q_res_down", q_res_down),
            ]
            # Put computed reserves into output dicts
            for (key, arr) in reserve_key_arrays
                sdd[key] = arr
            end
        end
    end

    return solution_data, projected_device_intervals
end
