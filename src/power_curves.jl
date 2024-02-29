function get_interval_times(input_data::NamedTuple)
    interval_times = Vector{NamedTuple}()
    start = 0.0
    stop = 0.0
    mid = 0.0
    for duration in input_data.dt
        start = stop
        stop += duration
        mid = 0.5*(start + stop)
        push!(interval_times, (start=start, stop=stop, mid=mid))
    end
    return interval_times
end

function get_power_curve_durations(input_data::NamedTuple)
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup
    periods = input_data.periods
    uids = input_data.sdd_ids

    rusu_key = "p_startup_ramp_ub"
    rdsd_key = "p_shutdown_ramp_ub"
    init_key = "initial_status"

    # Duration of a startup power curve *ending at i* is the p_on lower bound at
    # i divided by the maximum startup ramp up rate.
    #
    # Note that SUPCs *end* at the start of interval i, so i *should not* be
    # included.
    supc_duration = Dict(
        uid => [
            sdd_ts_lookup[uid]["p_lb"][i]/sdd_lookup[uid][rusu_key]
            for i in periods
        ] for uid in uids
    )

    # SDPC duration depends on *previous* minimum power, which we need
    # to pre-compute to handle the first interval properly
    prev_p_lb = Dict(
        uid => [
            # Use this instead of ifelse as it does not attempt to evaluate the
            # expressions unless necessary.
            i == 1 ? sdd_lookup[uid][init_key]["p"] : sdd_ts_lookup[uid]["p_lb"][i-1]
            for i in periods
        ] for uid in uids
    )
    # Duration of a shutdown power curve *starting at i* is the p_on lower
    # bound at *i-1* divided by the maximum shutdown ramp down rate.
    #
    # Note that SDPCs *start* at the start of interval i, so i *should* be
    # included.
    sdpc_duration = Dict(
        uid => [prev_p_lb[uid][i]/sdd_lookup[uid][rdsd_key] for i in periods]
        for uid in uids
    )
    return supc_duration, sdpc_duration
end

function get_contiguous_power_curves(input_data::NamedTuple; tolerance = 1e-8)
    interval_times = get_interval_times(input_data)
    supc_duration, sdpc_duration = get_power_curve_durations(input_data)

    uids = input_data.sdd_ids
    periods = input_data.periods

    # Maps UID to a vector containing the *start point* of a power curve
    # for startup at each coordinate.
    # E.g. if supc[uid][4] == 1, then if the device starts up at 4, it
    # must be in a power curve in intervals 1, 2, and 3.
    supc = Dict(uid => [i for i in periods] for uid in uids)

    # Maps UID to a vector containing the *end point* of a power curve
    # for shutdown at each coordinate.
    # E.g. if sdpc[uid][4] == 7, then if the device shuts down at 4, it
    # must be in a power curve in intervals 4, 5, 6, and 7.
    #
    # Note that the interval of SD *is* included in the power curve.
    # We initialize each coordinate to i-1 to signify an empty interval.
    sdpc = Dict(uid => [i-1 for i in periods] for uid in uids)

    # It makes sense to compute p_supc and p_sdpc at this point.
    sdd_lookup = input_data.sdd_lookup
    rusu_key = "p_startup_ramp_ub"
    rdsd_key = "p_shutdown_ramp_ub"
    p_rusu = Dict(uid => sdd_lookup[uid][rusu_key] for uid in uids)
    p_rdsd = Dict(uid => sdd_lookup[uid][rdsd_key] for uid in uids)
    p_supc = Dict{Tuple{String, Int64, Int64}, Float64}()
    p_sdpc = Dict{Tuple{String, Int64, Int64}, Float64}()

    ifirst = firstindex(interval_times)
    ilast = lastindex(interval_times)
    for uid in uids
        for (i, interval) in enumerate(interval_times)

            # Look backward to calculate the startup power curve
            #
            # Note that we *do not* include i, the time of startup, in the
            # power curve.
            for t in reverse(ifirst:i-1)
                t_interval = interval_times[t]
                delta_t = interval.stop - t_interval.stop

                # If we are within the PC duration of the SU time
                ramp_duration_t = supc_duration[uid][i] - delta_t
                if ramp_duration_t > tolerance
                    supc[uid][i] = t
                    p_supc[uid, t, i] = ramp_duration_t * p_rusu[uid]
                else
                    # ramp_duration_t decreases monotonically with t.
                    # Once it is <= 0, we can safely break.
                    break
                end
            end

            # Look forward to calculate shutdown power curve
            #
            # Note that we *do* include i, the time of shutdown, in the
            # power curve.
            for t in i : ilast
                t_interval = interval_times[t]
                delta_t = t_interval.stop - interval.start

                # If we are within the PC duration of the SD time
                left_to_ramp = sdpc_duration[uid][i] - delta_t
                if left_to_ramp > tolerance
                    sdpc[uid][i] = t
                    p_sdpc[uid, t, i] = left_to_ramp * p_rdsd[uid]
                else
                    # left_to_ramp decreases monotonically as t increases.
                    # Once it is <= 0, we can safely break.
                    break
                end
            end

        end
    end

    # Note that supc is "forward-looking" and sdpc is "backward-looking",
    # which are not the T_supc and T_sdpc required by the formulation.
    # These will be computed in a separate function.
    return supc, p_supc, sdpc, p_sdpc
end

"""
Returns u_su and u_sd indicators per the problem formulation. These are 1
on the first interval after a SU curve or the first interval of an SD
curve.
"""
function get_su_sd_from_on_status(input_data::NamedTuple, on_status::Dict)
    uids = input_data.sdd_ids
    periods = input_data.periods
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup
    init_key = "initial_status"
    on_key = "on_status"
    prev_on_status = Dict(
        uid => [
            i == 1 ? sdd_lookup[uid][init_key][on_key] : on_status[uid][i-1]
            for i in periods
        ] for uid in uids
    )
    u_su = Dict(uid => [0 for i in periods] for uid in uids)
    u_sd = Dict(uid => [0 for i in periods] for uid in uids)
    for uid in uids
        for i in periods
            u_on = Bool(round(on_status[uid][i]))
            u_on_prev = Bool(round(prev_on_status[uid][i]))
            if u_on && !u_on_prev
                u_su[uid][i] = 1
            elseif !u_on && u_on_prev
                u_sd[uid][i] = 1
            end
        end
    end
    return u_su, u_sd
end

"""
Returns indicators of whether we are in an SU/SD power curve and lookups
of the corresponding SU/SD power. Note the difference between these indicators
and u_su/u_sd from the problem formulation.
"""
function get_supc_sdpc_lookups(
    input_data::NamedTuple, on_status::Dict; tolerance = 1e-8
)
    supc, p_supc, sdpc, p_sdpc = get_contiguous_power_curves(
        input_data; tolerance = 1e-8
    )
    u_su, u_sd = get_su_sd_from_on_status(input_data, on_status)

    uids = input_data.sdd_ids
    periods = input_data.periods

    # These are the realized u_su, u_sd, p_su, and p_sd for a given
    # commitment schedule.
    supc_status = Dict(uid => [0 for i in periods] for uid in uids)
    sdpc_status = Dict(uid => [0 for i in periods] for uid in uids)
    p_su = Dict(uid => [0.0 for i in periods] for uid in uids)
    p_sd = Dict(uid => [0.0 for i in periods] for uid in uids)
    for uid in uids
        for i in periods

            # If we are starting up
            if Bool(u_su[uid][i])
                # Set the supc flag and the startup power 
                # for all intervals specified by the curve
                #
                # Note that the interval where u_su==1 (i) *is not* included
                # in the power curve.
                for t in range(supc[uid][i], i-1, step=1)
                    supc_status[uid][t] = 1
                    # This assertion makes sure we are not in multiple
                    # "power curves" simultaneously. If so, we have a bug.
                    @assert p_su[uid][t] == 0.0 && p_sd[uid][t] == 0.0
                    p_su[uid][t] = p_supc[uid, t, i]
                end
            end

            # If we are shutting down
            if Bool(u_sd[uid][i])
                # Set the sdpc flag and the shutdown power
                # for all intervals specified by the curve
                #
                # Note that the interval where u_sd==1 (i) *is* included
                # in the power curve.
                for t in range(i, sdpc[uid][i], step=1)
                    sdpc_status[uid][t] = 1
                    # This assertion makes sure we are not in multiple
                    # "power curves" simultaneously. If so, we have a bug.
                    @assert p_su[uid][t] == 0.0 && p_sd[uid][t] == 0.0
                    p_sd[uid][t] = p_sdpc[uid, t, i]
                end
            end

        end
    end
    return supc_status, p_su, sdpc_status, p_sd
end

function get_inverse_power_curve_intervals(
    input_data::NamedTuple; tolerance = 1e-8
)
    # The "power curve intervals" computed above are "standard", contiguous
    # power curve intervals. E.g. all intervals prior to a startup where
    # ramp-up power is needed. For startup at a given point in time, these
    # intervals are guaranteed to be contiguous.
    #
    # The MIP formulation requires "inverse" power curve intervals, e.g. all
    # future intervals where starup implies we need power now. These intervals
    # are *not* guaranteed to be contiguous.
    #
    # Naive computation:
    # for every generator
    #     for every time point
    #         for every following time point
    #             get PC duration at that time point
    #             if delta < duration, add to Tsupc
    #
    # Note that Tsupc needs to be a collection, not just an interval start/end
    # point.
    #
    # for every generator, for every time point
    #     initialize to an empty vector
    # for every generator, for every time point i
    #     for interval specified by supc (t)
    #         add t to Tsupc[uid][i]
    supc, _, sdpc, _ = get_contiguous_power_curves(
        input_data; tolerance = tolerance
    )
    periods = input_data.periods
    sdd_ids = input_data.sdd_ids
    Tsupc = Dict(uid => [Vector{Int64}() for i in periods] for uid in sdd_ids)
    Tsdpc = Dict(uid => [Vector{Int64}() for i in periods] for uid in sdd_ids)
    for uid in sdd_ids
        for i in periods
            # i is the point at which we start up. t is the point where we
            # must have some startup power.
            for t in supc[uid][i] : i-1
                push!(Tsupc[uid][t], i)
            end
            # i is the point where we shut down. t is the point where we must
            # have some shutdown power.
            # Note that t can equal i.
            for t in i : sdpc[uid][i]
                push!(Tsdpc[uid][t], i)
            end
        end
    end
    return Tsupc, Tsdpc
end
