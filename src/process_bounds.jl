"""
TODO: This function is intended to tighten bounds before the scheduling problem.
It will need to know the power curve for each device so feasible power regions
consider regions that are feasible wrt SU/SD ramp limits as well.
"""
function tighten_bounds_using_ramp_limits(input_data::NamedTuple)
end

"""
"""
function tighten_bounds_using_ramp_limits(
    input_data::NamedTuple,
    on_status::Dict,
    real_power::Dict;
    tolerance = 1e-8,
)
    sdd_ts_lookup = input_data.sdd_ts_lookup

    # Initialize data structures that will hold tightened bounds
    p_lb = Dict{String, Vector{Float64}}(
        uid => [sdd_ts_lookup[uid]["p_lb"][i] for i in input_data.periods]
        for uid in input_data.sdd_ids
    )
    p_ub = Dict{String, Vector{Float64}}(
        uid => [sdd_ts_lookup[uid]["p_ub"][i] for i in input_data.periods]
        for uid in input_data.sdd_ids
    )

    for uid in input_data.sdd_ids
        uid_lb = p_lb[uid]
        uid_ub = p_ub[uid]
        ru_ub = input_data.sdd_lookup[uid]["p_ramp_up_ub"]
        rd_ub = input_data.sdd_lookup[uid]["p_ramp_down_ub"]
        su_ru_ub = input_data.sdd_lookup[uid]["p_startup_ramp_ub"]
        sd_rd_ub = input_data.sdd_lookup[uid]["p_shutdown_ramp_ub"]

        u_prev = input_data.sdd_lookup[uid]["initial_status"]["on_status"]
        # This will only be used to get the last/first value in a SU/SD
        # power curve.
        p_prev = input_data.sdd_lookup[uid]["initial_status"]["p"]

        prev_lb = p_prev
        prev_ub = p_prev

        #
        # Forward pass
        #
        for i in input_data.periods
            u = on_status[uid][i]
            p = real_power[uid][i]

            if u_prev > tolerance && u > tolerance
                # If we are on, and were on in the previous interval, we must
                # respect ru_ub and rd_ub
                ub = prev_ub + ru_ub*input_data.dt[i]
                lb = prev_lb - rd_ub*input_data.dt[i]
                if ub < p_ub[uid][i] - tolerance
                    p_ub[uid][i] = ub
                end
                if lb > p_lb[uid][i] + tolerance
                    p_lb[uid][i] = lb
                end
                #p_ub[uid][i] = min(ub, p_ub[uid][i])
                #p_lb[uid][i] = max(lb, p_lb[uid][i])
                p_lb[uid][i], p_ub[uid][i] = check_bounds_not_crossing(
                    # I believe using the same tolerance crossing bounds as for
                    # equality with existing bounds (above) is appropriate.
                    p_lb[uid][i], p_ub[uid][i], tolerance = tolerance
                )
            elseif u_prev <= tolerance && u > tolerance
                # If we are in start-up, we must respect su_ru_ub
                ub = p_prev + su_ru_ub*input_data.dt[i]
                # The ramp down bound is still valid, even though we should
                # never ramp down on startup.
                lb = p_prev - rd_ub*input_data.dt[i]
                if ub < p_ub[uid][i] - tolerance
                    p_ub[uid][i] = ub
                end
                if lb > p_lb[uid][i] + tolerance
                    p_lb[uid][i] = lb
                end
                #p_ub[uid][i] = min(ub, p_ub[uid][i])
                #p_lb[uid][i] = max(lb, p_lb[uid][i])
                p_lb[uid][i], p_ub[uid][i] = check_bounds_not_crossing(
                    p_lb[uid][i], p_ub[uid][i], tolerance = tolerance
                )
            end

            # Don't need to attempt to propagate bounds backwards. This will
            # happen in the reverse pass.

            # Set previous u/p/lb/ub
            u_prev = u
            p_prev = p
            prev_ub = p_ub[uid][i]
            prev_lb = p_lb[uid][i]
        end

        #
        # Reverse pass
        #
        last = lastindex(input_data.periods)
        u_next = on_status[uid][last]
        p_next = real_power[uid][last]
        next_ub = p_ub[uid][last]
        next_lb = p_lb[uid][last]
        # Note that I reverse, then apply the slice index.
        for i in reverse(input_data.periods)[2:last]
            u = on_status[uid][i]
            p = real_power[uid][i]

            # If we are on, and will be on in the next interval, we must
            # respect ru_ub and rd_ub
            if u_next > tolerance && u > tolerance
                # Now the UB needs to respect the ramp-down limit, and the
                # LB needs to respect the ramp-up limit
                ub = next_ub + rd_ub*input_data.dt[i+1]
                lb = next_lb - ru_ub*input_data.dt[i+1]
                if ub < p_ub[uid][i] - tolerance
                    p_ub[uid][i] = ub
                end
                if lb > p_lb[uid][i] + tolerance
                    p_lb[uid][i] = lb
                end
                #p_ub[uid][i] = min(ub, p_ub[uid][i])
                #p_lb[uid][i] = max(lb, p_lb[uid][i])
                p_lb[uid][i], p_ub[uid][i] = check_bounds_not_crossing(
                    p_lb[uid][i], p_ub[uid][i], tolerance = tolerance
                )
            elseif u_next <= tolerance && u > tolerance
                # We are about to shut down and need to resepct sd_rd_ub
                ub = p_next + sd_rd_ub*input_data.dt[i+1]
                lb = p_next - ru_ub*input_data.dt[i+1]
                if ub < p_ub[uid][i] - tolerance
                    p_ub[uid][i] = ub
                end
                if lb > p_lb[uid][i] + tolerance
                    p_lb[uid][i] = lb
                end
                #p_ub[uid][i] = min(ub, p_ub[uid][i])
                #p_lb[uid][i] = max(lb, p_lb[uid][i])
                p_lb[uid][i], p_ub[uid][i] = check_bounds_not_crossing(
                    p_lb[uid][i], p_ub[uid][i], tolerance = tolerance
                )
            end

            u_next = u
            p_next = p
            next_ub = p_ub[uid][i]
            next_lb = p_lb[uid][i]
        end
    end

    # TODO: Performance implications of using deepcopy here? We almost
    # definitely don't need to copy everything.
    tightened_data = deepcopy(input_data)
    # Update bounds in sdd_ts_lookup with tightened bounds we just computed.
    for uid in tightened_data.sdd_ids
        tightened_data.sdd_ts_lookup[uid]["p_lb"] = p_lb[uid]
        tightened_data.sdd_ts_lookup[uid]["p_ub"] = p_ub[uid]
    end
    return tightened_data
end


function check_bounds_not_crossing(lb, ub; tolerance = 1e-8)
    if lb > ub && lb - ub <= tolerance
        # lb and ub are equal within tolerance
        # The convention here is to convert both bounds to the average between
        # them.
        bound = 0.5*(lb + ub)
        return bound, bound
    elseif lb > ub
        println("lb: $lb, ub: $ub, tolerance = $tolerance")
        throw(Exception)
    else
        # lb <= ub
        return lb, ub 
    end
end
