#module SchedulingModel

import JuMP as jmp
using JSON

using HiGHS

# Note that both of these files are being included in multiple
# locations. Is this bad? What is the alternative?
# - Appears to cause a warning "Replacing docs..."
#include("power_curves.jl")
#include("reserves.jl")

function add_power_variables!(model, sdd_ts_lookup, periods, sdd_ids)
    p = @variable(model, 
        min(sdd_ts_lookup[uid]["p_lb"][i], 0.0) <=
        p[uid in sdd_ids, i in periods]
        <= max(sdd_ts_lookup[uid]["p_ub"][i], 0.0)
    )
    q = @variable(
        model,
        # Do these bounds make sense for reactive power?
        min(sdd_ts_lookup[uid]["q_lb"][i], 0.0) <=
        q[uid in sdd_ids, i in periods]
        <= max(sdd_ts_lookup[uid]["q_ub"][i], 0.0)
    )
    return p, q
end

function add_p_on_status_variables!(model, sdd_ts_lookup, periods, sdd_ids)
    p_on_status = @variable(model, 
        sdd_ts_lookup[uid]["on_status_lb"][i] <=
        p_on_status[uid in sdd_ids, i in periods]
        <= sdd_ts_lookup[uid]["on_status_ub"][i],
        binary = true
    )
    return p_on_status
end

"""
The parts of (44-45) associated with initial uptime/downtime
"""
function add_init_mustrun_outage_constraints!(
    model, sdd_lookup, periods, sdd_ids, interval_times;
    tolerance = 1e-8
)
    # WARNING: This constraint compares computed floating point values
    # without a tolerance. This is unsafe, and could cause bugs.
    # TODO: Update comparisons to respect equality-within-tolerance
    Tout0 = Dict{String, Any}(uid => Vector() for uid in sdd_ids)
    Tmr0 = Dict{String, Any}(uid => Vector() for uid in sdd_ids)
    for uid in sdd_ids
        Tjout0 = Tout0[uid]
        Tjmr0 = Tmr0[uid]
        init_dtime = sdd_lookup[uid]["initial_status"]["accu_down_time"]
        init_utime = sdd_lookup[uid]["initial_status"]["accu_up_time"]
        min_dtime = sdd_lookup[uid]["down_time_lb"]
        min_utime = sdd_lookup[uid]["in_service_time_lb"]
        remaining_dtime = min_dtime - init_dtime
        remaining_utime = min_utime - init_utime
        @assert init_dtime > 0 || init_utime > 0
        @assert init_dtime == 0 || init_utime == 0
        if init_dtime > 0
            # We start with some downtime. We need to not violate minimum
            # downtime.
            for i in periods
                a = interval_times[i].start
                if a < remaining_dtime - tolerance
                    push!(Tjout0, i)
                else
                    break
                end
            end
            # We started with some down time. No need to enforce initial
            # uptime
            delete!(Tmr0, uid)
            # If we have already satisfied the minimum downtime, delete this
            # key so we don't have an empty constraint.
            if length(Tjout0) == 0
                delete!(Tout0, uid)
            end
        elseif init_utime > 0
            # We start with some uptime. Need to not violate minimum.
            for i in periods
                a = interval_times[i].start
                if a < remaining_utime - tolerance
                    push!(Tjmr0, i)
                else
                    break
                end
            end
            delete!(Tout0, uid)
            if length(Tjmr0) == 0
                delete!(Tmr0, uid)
            end
        end
    end

    outage_init_uids = keys(Tout0)
    outage_init_uid_t = [
        (uid, t) for uid in outage_init_uids for t in Tout0[uid]
    ]
    mustrun_init_uids = keys(Tmr0)
    mustrun_init_uid_t = [
        (uid, t) for uid in mustrun_init_uids for t in Tmr0[uid]
    ]
    if length(outage_init_uid_t) > 0
        outage_init = @constraint(
            model,
            outage_init[(uid, t) in outage_init_uid_t],
            model[:p_on_status][uid, t] == 0,
        )
    else
        outage_init = nothing
    end

    if length(mustrun_init_uid_t) > 0
        mustrun_init = @constraint(
            model,
            mustrun_init[(uid, t) in mustrun_init_uid_t],
            model[:p_on_status][uid, t] == 1,
        )
    else
        mustrun_init = nothing
    end

    return mustrun_init, outage_init
end

"""
Constraints (46-48)
"""
function add_su_sd_variables_constraints!(model, sdd_lookup, periods, sdd_ids)
    # TODO: This constraint will need to be applied to ac lines when we relax
    # the copper plate assumption

    # TODO: bounds for SU/SD?
    u_su = @variable(model, u_su[u in sdd_ids, i in periods], binary=true)
    u_sd = @variable(model, u_sd[u in sdd_ids, i in periods], binary=true)

    u = model[:p_on_status]
    prev_u_lookup = Dict(
        (uid, i) => (i == 1)
            ? sdd_lookup[uid]["initial_status"]["on_status"]
            : u[uid, i-1]
        for uid in sdd_ids for i in periods
    )
    on_status_evolution = @constraint(
        model,
        on_status_evolution[uid in sdd_ids, i in periods],
        u[uid, i] - prev_u_lookup[uid, i] == u_su[uid, i] - u_sd[uid, i],
    )
    prohibit_su_sd = @constraint(
        model,
        prohibit_su_sd[uid in sdd_ids, i in periods],
        u_su[uid, i] + u_sd[uid, i] <= 1,
    )
    return u_su, u_sd, on_status_evolution, prohibit_su_sd
end

"""
Constraints (56-57)
"""
function add_min_up_down_constraints!(model, sdd_lookup, periods, sdd_ids, interval_times)
    # WARNING: This constraint compares computed floating point values
    # without a tolerance. This is unsafe, and could cause bugs.
    # TODO: Update comparisons to respect equality-within-tolerance

    # These map UIDs to an array of lower bounds of the represented
    # intervals. The upper bound of the Tdnmin interval at position i is
    # the start of period i.
    Tdnmin_lb = Dict{String, Vector{Int64}}(
        uid => Vector(periods) for uid in sdd_ids
    )
    Tupmin_lb = Dict{String, Vector{Int64}}(
        uid => Vector(periods) for uid in sdd_ids
    )

    for uid in sdd_ids
        min_dn = sdd_lookup[uid]["down_time_lb"]
        min_up = sdd_lookup[uid]["in_service_time_lb"]
        for i in periods
            # We set the lower bound equal to the upper bound (i)
            # to represent an empty interval. Need to be careful when
            # summing over Tdnmin_lb_i:i.
            Tdnmin_lb_i = i
            Tupmin_lb_i = i
            start = interval_times[i].start
            # Walk backwards to find all other intervals in the
            # min up/downtime windows
            #
            # Note that binary search may be more efficient here.
            # Also, since the same min up/down times may repeat,
            # it may improve performance to implement a cache.
            for t in reverse(1:i-1)
                diff = start - interval_times[t].start
                if diff < min_dn
                    Tdnmin_lb_i = t
                end
                if diff < min_up
                    Tupmin_lb_i = t
                end
                if diff >= min_dn && diff >= min_up
                    break
                end
            end
            Tdnmin_lb[uid][i] = Tdnmin_lb_i
            Tupmin_lb[uid][i] = Tupmin_lb_i
        end
    end
    u_su = model[:u_su]
    u_sd = model[:u_sd]
    min_up = @constraint(
        model,
        min_up[uid in sdd_ids, i in periods],
        u_sd[uid, i] <= 1 - sum(u_su[uid, t] for t in Tupmin_lb[uid][i] : i-1, init=0)
    )
    min_dn = @constraint(
        model,
        min_dn[uid in sdd_ids, i in periods],
        u_su[uid, i] <= 1 - sum(u_sd[uid, t] for t in Tdnmin_lb[uid][i] : i-1, init=0)
    )
    return min_up, min_dn
end

"""
Constraint (58)
"""
function add_maximum_starts_over_intervals!(
    model, sdd_lookup, periods, sdd_ids, interval_times;
    tolerance=1e-8,
)
    su_ub = "startups_ub"
    i_init = first(periods)
    i_final = last(periods)
    # Need to construct the max_su_info data structure from interval times
    # Need to locate the intervals where w_start and w_end occur
    # Binary search seems to be the best option here.
    interval_starts = [a.start for a in interval_times]
    max_su_info = Dict{Tuple{String, Int64}, NamedTuple}()
    for uid in sdd_ids
        for w in 1:length(sdd_lookup[uid][su_ub])
            w_a_start = sdd_lookup[uid][su_ub][w][1]
            w_a_stop = sdd_lookup[uid][su_ub][w][2]

            # From the Julia docs: Returns first index *greater than or equal to*.
            # This is what we want
            w_i_start = searchsortedfirst(interval_starts, w_a_start)

            # Postprocess to correct for potential roundoff error
            if w_i_start > 1 && w_a_start - interval_starts[w_i_start-1] <= tolerance
                @assert w_a_start - interval_starts[w_i_start-1] > 0
                # We are close enough to the previous interval start that we should
                # return use it as the start of our "restricted interval". This can
                # happen due to floating point error.
                w_i_start = w_i_start - 1
            end

            # This can be violated if w_a_start > all interval starts.
            # We should return an empty interval. This happens in (0, 1, 3, 2).
            if w_i_start > lastindex(interval_starts)
                max_su_info[uid, w] = (
                    interval=[],
                    max_su=sdd_lookup[uid][su_ub][w][3],
                )
                break
            end

            # Returns the last index *less than or equal to* w_a_stop. We want
            # a_start < w_a_stop
            w_i_stop = searchsortedlast(interval_starts, w_a_stop)

            # This will be violated if w_a_stop < all interval starts.
            # In this case the restricted interval is empty. This should
            # be checked *before and after* processing with tolerance.
            # See above for how to return an empty interval. I don't think
            # this ever happens, however, because it would require w_a_stop
            # to be negative.
            @assert w_i_stop >= firstindex(interval_starts)

            # Postprocess to correct for potentially equal (within tolerance)
            # values.
            @assert interval_starts[w_i_stop] <= w_a_stop
            if w_a_stop - interval_starts[w_i_stop] <= tolerance
                w_i_stop -= 1
                @assert w_i_stop >= firstindex(interval_starts)
            end

            # This is not necessarily the case. if w_a_start==w_a_stop (or they
            # are both strictly within the same interval, these indices will cross)
            #@assert w_i_start <= w_i_stop
            w_interval = w_i_start:w_i_stop

            max_su_info[uid, w] = (
                interval=w_interval,
                max_su=sdd_lookup[uid][su_ub][w][3],
            )
        end
    end

    device_intervals = sort([k for k in keys(max_su_info)])
    u_su = model[:u_su]
    max_starts_over_intervals = @constraint(
        model,
        max_starts_over_intervals[(uid, w) in device_intervals],
        # Use init=0 in case we have an empty interval. This can happen
        # if w_a_start == w_a_stop, or the entire restricted interval is
        # outside our horizon.
        sum(u_su[uid, t] for t in max_su_info[uid, w].interval, init=0)
            <= max_su_info[uid, w].max_su,
    )
    return max_starts_over_intervals
end

function _get_max_min_energy_interval_info(
    sdd_lookup, periods, sdd_ids, interval_times; tolerance=1e-8
)
    max_en_info = Dict{Tuple{String, Int64}, NamedTuple}()
    min_en_info = Dict{Tuple{String, Int64}, NamedTuple}()

    # Populate dicts with intervals and max/min energy
    i_init = firstindex(interval_times)
    i_final = lastindex(interval_times)

    max_en_key = "energy_req_ub"
    min_en_key = "energy_req_lb"

    interval_mids = [a.mid for a in interval_times]

    # Finds the first coordinate strictly greater than t, within tolerance
    function _binary_search_gt(array, t; tolerance=1e-8)
        i = searchsortedfirst(array, t)
        @assert i > i_final || array[i] >= t
        if i <= i_final && array[i] - t <= tolerance
            # If we are equal-within-tolerance to some array element, use
            # the next element. If this is i_final+1, the interval will be
            # empty.
            i += 1
        end
        return i
    end

    # Finds the last coordinate less than or equal to t, within tolerance
    function _binary_search_leq(array, t; tolerance=1e-8)
        i = searchsortedlast(array, t)
        @assert i < i_init || array[i] <= t
        if i < i_final && array[i+1] - t <= tolerance
            @assert t < array[i+1]
            i += 1
        end
        return i
    end

    for uid in sdd_ids
        for (w, (start, stop, max_en)) in enumerate(sdd_lookup[uid][max_en_key])
            # Note that these indices could be outside the range of
            # [i_init, i_final]. In these cases, however, we will get empty
            # intervals and never access variables at these invalid indices.
            w_i_start = _binary_search_gt(
                interval_mids, start; tolerance=tolerance
            )
            w_i_stop = _binary_search_leq(
                interval_mids, stop; tolerance=tolerance
            )
            # Note that this interval could be empty
            max_en_info[uid, w] = (interval=(w_i_start:w_i_stop), max_en=max_en)
        end
        for (w, (start, stop, min_en)) in enumerate(sdd_lookup[uid][min_en_key])
            w_i_start = _binary_search_gt(
                interval_mids, start; tolerance=tolerance
            )
            w_i_stop = _binary_search_leq(
                interval_mids, stop; tolerance=tolerance
            )
            min_en_info[uid, w] = (interval=(w_i_start:w_i_stop), min_en=min_en)
        end
    end
    return max_en_info, min_en_info
end

"""
Constraints (68-69)
"""
function add_max_min_energy_over_intervals!(
    model, sdd_lookup, periods, sdd_ids, interval_times, dt;
    tolerance=1e-8,
)
    max_en_info, min_en_info = _get_max_min_energy_interval_info(
        sdd_lookup, periods, sdd_ids, interval_times; tolerance=tolerance
    )
    # Construct constraints from processed data
    p = model[:p]
    uid_intervals = sort([k for k in keys(max_en_info)])
    max_en_over_intervals = @constraint(
        model,
        max_energy_over_intervals[(uid, w) in uid_intervals],
        sum(dt[t]*p[uid, t] for t in max_en_info[uid, w].interval, init=0)
            <= max_en_info[uid, w].max_en,
    )
    uid_intervals = sort([k for k in keys(min_en_info)])
    min_en_over_intervals = @constraint(
        model,
        min_energy_over_intervals[(uid, w) in uid_intervals],
        sum(dt[t]*p[uid, t] for t in min_en_info[uid, w].interval, init=0)
            >= min_en_info[uid, w].min_en,
    )
    return max_en_over_intervals, min_en_over_intervals
end

"""
Constraints (62-63)
"""
function add_startup_shutdown_power!(
    model, sdd_ts_lookup, periods, sdd_ids, p_su_pc, T_su_pc, p_sd_pc, T_sd_pc
)
    p_su = @variable(model,
        min(sdd_ts_lookup[uid]["p_lb"][i], 0.0) <=
        p_su[uid in sdd_ids, i in periods]
        <= max(sdd_ts_lookup[uid]["p_ub"][i], 0.0)
    )
    p_sd = @variable(model,
        min(sdd_ts_lookup[uid]["p_lb"][i], 0.0) <=
        p_sd[uid in sdd_ids, i in periods]
        <= max(sdd_ts_lookup[uid]["p_ub"][i], 0.0)
    )
    u_su = model[:u_su]
    u_sd = model[:u_sd]
    su_power_curve = @constraint(
        model,
        su_power_curve[uid in sdd_ids, i in periods],
        p_su[uid, i] == sum(
            p_su_pc[uid, i, t]*u_su[uid, t] for t in T_su_pc[uid][i],
            init=0
        ),
    )
    sd_power_curve = @constraint(
        model,
        sd_power_curve[uid in sdd_ids, i in periods],
        p_sd[uid, i] == sum(
            p_sd_pc[uid, i, t]*u_sd[uid, t] for t in T_sd_pc[uid][i],
            init=0
        ),
    )
    return p_su, p_sd, su_power_curve, sd_power_curve
end

"""
Constraint (61). This also adds the p_on variable.
"""
function add_p_on_su_sd_aggregation!(model, sdd_ts_lookup, sdd_ids, periods)
    p_on = @variable(model,
        min(sdd_ts_lookup[uid]["p_lb"][i], 0.0) <=
        p_on[uid in sdd_ids, i in periods]
        <= max(sdd_ts_lookup[uid]["p_ub"][i], 0.0)
    )
    p = model[:p]
    p_su = model[:p_su]
    p_sd = model[:p_sd]
    p_on_su_sd_aggregation = @constraint(
        model,
        p_disaggregation[uid in sdd_ids, i in periods],
        p[uid, i] == p_on[uid, i] + p_su[uid, i] + p_sd[uid, i],
    )
    return p_on, p_on_su_sd_aggregation
end

"""
Constraints (64-67)
"""
function add_ramp_constraints!(
    model, sdd_lookup, periods, sdd_ids, dt;
    p=nothing,
    u_on=nothing,
    u_su=nothing,
    u_sd=nothing,
    tightening_fraction = 1.0
)
    if p === nothing
        p = model[:p]
    end
    # Note that this constraint doesn't need sdd_ts_lookup
    prev_p_lookup = Dict(
        (uid, i) => (i == 1) ? sdd_lookup[uid]["initial_status"]["p"] : p[uid, i-1]
        for uid in sdd_ids for i in periods
    )
    pjru = Dict(uid => sdd_lookup[uid]["p_ramp_up_ub"] for uid in sdd_ids)
    pjrd = Dict(uid => sdd_lookup[uid]["p_ramp_down_ub"] for uid in sdd_ids)
    pjrusu = Dict(uid => sdd_lookup[uid]["p_startup_ramp_ub"] for uid in sdd_ids)
    pjrdsd = Dict(uid => sdd_lookup[uid]["p_shutdown_ramp_ub"] for uid in sdd_ids)
    if u_on === nothing
        u_on = model[:p_on_status]
    end
    if u_su === nothing
        u_su = model[:u_su]
    end
    if u_sd === nothing
        u_sd = model[:u_sd]
    end
    ramp_ub = @constraint(
        model,
        ramp_ub[uid in sdd_ids, i in periods],
        p[uid, i] - prev_p_lookup[uid, i] <= dt[i]*(
            tightening_fraction*pjru[uid]*(u_on[uid, i] - u_su[uid, i])
            + pjrusu[uid]*(u_su[uid, i] + 1 - u_on[uid, i])
        ),
    )
    ramp_lb = @constraint(
        model,
        ramp_lb[uid in sdd_ids, i in periods],
        p[uid, i] - prev_p_lookup[uid, i] >= -dt[i]*(
            tightening_fraction*pjrd[uid]*u_on[uid, i] + pjrdsd[uid]*(1 - u_on[uid, i])
        ),
    )
    return ramp_ub, ramp_lb
end

"""
Constraints (99-101) and (109-111)
"""
function add_on_su_sd_implication_constraints!(
    model, input_data;
    include_reserves = false,
    p_on = nothing,
    u_on = nothing,
    p_su = nothing,
    p_sd = nothing,
)
    sdd_ts_lookup = input_data.sdd_ts_lookup
    periods = input_data.periods
    sdd_ids = input_data.sdd_ids
    sdd_ids_producer = input_data.sdd_ids_producer
    sdd_ids_consumer = input_data.sdd_ids_consumer
    if u_on === nothing
        u_on = model[:p_on_status]
    end
    if p_on === nothing
        p_on = model[:p_on]
    end
    if p_su === nothing
        p_su = model[:p_su]
    end
    if p_sd === nothing
        p_sd = model[:p_sd]
    end
    if include_reserves
        p_rgu = model[:p_rgu]
        p_rgd = model[:p_rgd]
        p_scr = model[:p_scr]
        p_nsc = model[:p_nsc]
        p_rru_on = model[:p_rru_on]
        p_rrd_on = model[:p_rrd_on]
        p_rru_off = model[:p_rru_off]
        p_rrd_off = model[:p_rrd_off]
        q_qru = model[:q_qru]
        q_qrd = model[:q_qrd]
        reserve_ubs = _get_max_reserves(input_data)
        (p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
         p_rru_off_max, p_rrd_off_max) = reserve_ubs
        # "Relative reserve limits"
        @constraint(model,
            p_on_ub_implication_prod[uid in sdd_ids_producer, i in periods],
            (
                p_on[uid, i] + p_rgu[uid, i] + p_scr[uid, i] + p_rru_on[uid, i]
            ) <= sdd_ts_lookup[uid]["p_ub"][i] * u_on[uid, i],
        )
        @constraint(model,
            p_on_ub_implication_cons[uid in sdd_ids_consumer, i in periods],
            (
                p_on[uid, i] + p_rgd[uid, i] + p_rrd_on[uid, i]
            ) <= sdd_ts_lookup[uid]["p_ub"][i] * u_on[uid, i],
        )
        @constraint(model,
            p_on_lb_implication_prod[uid in sdd_ids_producer, i in periods],
            sdd_ts_lookup[uid]["p_lb"][i] * u_on[uid, i] <= (
                p_on[uid, i] - p_rgd[uid, i] - p_rrd_on[uid, i]
            )
        )
        @constraint(model,
            p_on_lb_implication_cons[uid in sdd_ids_consumer, i in periods],
            sdd_ts_lookup[uid]["p_lb"][i] * u_on[uid, i] <= (
                p_on[uid, i] - p_rgu[uid, i] - p_scr[uid, i] - p_rru_on[uid, i]
            )
        )
        @constraint(model,
            p_su_sd_ub_implication_prod[uid in sdd_ids_producer, i in periods],
            p_su[uid, i] + p_sd[uid, i] + p_nsc[uid, i] + p_rru_off[uid, i] <= (
                sdd_ts_lookup[uid]["p_ub"][i] * (1 - u_on[uid, i])
            )
        )
        @constraint(model,
            p_su_sd_ub_implication_cons[uid in sdd_ids_consumer, i in periods],
            p_su[uid, i] + p_sd[uid, i] + p_rrd_off[uid, i] <= (
                sdd_ts_lookup[uid]["p_ub"][i] * (1 - u_on[uid, i])
            )
        )
        # "Absolute reserve limits"
        @constraint(model,
            p_rgu_implication[uid in sdd_ids, i in periods],
            p_rgu[uid, i] <= p_rgu_max[uid]*u_on[uid, i]
        )
        @constraint(model,
            p_rgd_implication[uid in sdd_ids, i in periods],
            p_rgd[uid, i] <= p_rgd_max[uid]*u_on[uid, i]
        )
        @constraint(model,
            p_scr_implication[uid in sdd_ids, i in periods],
            p_rgu[uid, i] + p_scr[uid, i] <= p_scr_max[uid]*u_on[uid, i]
        )
        @constraint(model,
            p_nsc_implication[uid in sdd_ids, i in periods],
            p_nsc[uid, i] <= p_nsc_max[uid]*(1 - u_on[uid, i])
        )
        @constraint(model,
            p_rru_on_implication[uid in sdd_ids, i in periods],
            p_rgu[uid, i] + p_scr[uid, i] + p_rru_on[uid, i] <= (
                p_rru_on_max[uid]*u_on[uid, i]
            )
        )
        @constraint(model,
            p_rru_off_implication[uid in sdd_ids, i in periods],
            p_nsc[uid, i] + p_rru_off[uid, i] <= (
                p_rru_off_max[uid]*(1 - u_on[uid, i])
            )
        )
        @constraint(model,
            p_rrd_on_implication[uid in sdd_ids, i in periods],
            p_rgd[uid, i] + p_rrd_on[uid, i] <= p_rrd_on_max[uid]*u_on[uid, i]
        )
        @constraint(model,
            p_rrd_off_implication[uid in sdd_ids, i in periods],
            p_rrd_off[uid, i] <= p_rrd_off_max[uid]*(1 - u_on[uid, i])
        )
        for uid in sdd_ids_producer
            for i in periods
                jmp.fix(p_rrd_off[uid, i], 0.0; force = true)
            end
        end
        for uid in sdd_ids_consumer
            for i in periods
                jmp.fix(p_nsc[uid, i], 0.0; force = true)
                jmp.fix(p_rru_off[uid, i], 0.0; force = true)
            end
        end
        return nothing
    else
        p_on_ub_implication = @constraint(
            model,
            p_on_ub_implication[uid in sdd_ids, i in periods],
            p_on[uid, i] <= sdd_ts_lookup[uid]["p_ub"][i] * u_on[uid, i],
        )
        p_on_lb_implication = @constraint(
            model,
            p_on_lb_implication[uid in sdd_ids, i in periods],
            sdd_ts_lookup[uid]["p_lb"][i] * u_on[uid, i] <= p_on[uid, i],
        )
        p_su_sd_ub_implication = @constraint(
            model,
            p_su_sd_ub_implication[uid in sdd_ids, i in periods],
            p_su[uid, i] + p_sd[uid, i] <= (
                sdd_ts_lookup[uid]["p_ub"][i] * (1 - u_on[uid, i])
            )
        )
        return p_on_ub_implication, p_on_lb_implication, p_su_sd_ub_implication
    end
end

"""
Constraints (102-103) and (112-113)
"""
function add_reactive_power_implication_constraints!(
    model, input_data, T_su_pc, T_sd_pc;
    include_reserves = false,
    q=nothing,
    u_on=nothing,
    u_su=nothing,
    u_sd=nothing,
)
    sdd_ts_lookup = input_data.sdd_ts_lookup
    periods = input_data.periods
    sdd_ids = input_data.sdd_ids
    sdd_ids_consumer = input_data.sdd_ids_consumer
    sdd_ids_producer = input_data.sdd_ids_producer
    q_max = Dict(uid => sdd_ts_lookup[uid]["q_ub"] for uid in sdd_ids)
    q_min = Dict(uid => sdd_ts_lookup[uid]["q_lb"] for uid in sdd_ids)
    if q === nothing
        q = model[:q]
    end
    if u_on === nothing
        u_on = model[:p_on_status]
    end
    if u_su === nothing
        u_su = model[:u_su]
    end
    if u_sd === nothing
        u_sd = model[:u_sd]
    end
    if include_reserves
        q_qru = model[:q_qru]
        q_qrd = model[:q_qrd]
        @constraint(model,
            q_implication_max_prod[uid in sdd_ids_producer, i in periods],
            q[uid, i] + q_qru[uid, i] <= q_max[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        @constraint(model,
            q_implication_max_cons[uid in sdd_ids_consumer, i in periods],
            q[uid, i] + q_qrd[uid, i] <= q_max[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        @constraint(model,
            q_implication_min_prod[uid in sdd_ids_producer, i in periods],
            q[uid, i] - q_qrd[uid, i] >= q_min[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        @constraint(model,
            q_implication_min_cons[uid in sdd_ids_consumer, i in periods],
            q[uid, i] - q_qru[uid, i] >= q_min[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        return nothing
    else
        q_implication_max = @constraint(
            model,
            q_implication_max[uid in sdd_ids, i in periods],
            q[uid, i] <= q_max[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        q_implication_min = @constraint(
            model,
            q_implication_min[uid in sdd_ids, i in periods],
            q[uid, i] >= q_min[uid][i] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            )
        )
        return q_implication_max, q_implication_min
    end
end

"""
Constraints (104-106) and (114-116)
"""
# TODO: These constraints will need to be updated (differentiate between
# producers and consumers) when qrd/qru start being used.
function add_real_reactive_linking_constraints!(
    model,
    input_data,
    T_su_pc,
    T_sd_pc;
    pq = nothing,
    u_on = nothing,
    u_susd = nothing,
    include_reserves = false
)
    sdd_lookup = input_data.sdd_lookup
    periods = input_data.periods
    sdd_ids = input_data.sdd_ids
    sdd_ids_producer = input_data.sdd_ids_producer
    sdd_ids_consumer = input_data.sdd_ids_consumer
    producer_set = Set(sdd_ids_producer)
    consumer_set = Set(sdd_ids_consumer)
    q_bound = Dict(uid => sdd_lookup[uid]["q_bound_cap"] for uid in sdd_ids)
    q_linear = Dict(uid => sdd_lookup[uid]["q_linear_cap"] for uid in sdd_ids)
    pq_bound_sdds = [uid for uid in sdd_ids if q_bound[uid] == 1]
    pq_bound_sdds_prod = filter(uid -> uid in producer_set, pq_bound_sdds)
    pq_bound_sdds_cons = filter(uid -> uid in consumer_set, pq_bound_sdds)
    pq_linear_sdds = [uid for uid in sdd_ids if q_linear[uid] == 1]
    q_0_ub = Dict(uid => sdd_lookup[uid]["q_0_ub"] for uid in pq_bound_sdds)
    q_0_lb = Dict(uid => sdd_lookup[uid]["q_0_lb"] for uid in pq_bound_sdds)
    beta_ub = Dict(uid => sdd_lookup[uid]["beta_ub"] for uid in pq_bound_sdds)
    beta_lb = Dict(uid => sdd_lookup[uid]["beta_lb"] for uid in pq_bound_sdds)
    q_0 = Dict(uid => sdd_lookup[uid]["q_0"] for uid in pq_linear_sdds)
    beta = Dict(uid => sdd_lookup[uid]["beta"] for uid in pq_linear_sdds)

    # Get default variables from model if not explicitly provided
    if pq === nothing
        p = model[:p]
        q = model[:q]
    else
        p, q = pq
    end
    if u_on === nothing
        u_on = model[:p_on_status]
    end
    if u_susd === nothing
        u_su = model[:u_su]
        u_sd = model[:u_sd]
    else
        u_su, u_sd = u_susd
    end

    if include_reserves
        q_qru = model[:q_qru]
        q_qrd = model[:q_qrd]
        @constraint(model,
            pq_ub_prod[uid in pq_bound_sdds_prod, i in periods],
            q[uid, i] + q_qru[uid, i] <= q_0_ub[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_ub[uid]*p[uid, i]
        )
        @constraint(model,
            pq_ub_cons[uid in pq_bound_sdds_cons, i in periods],
            q[uid, i] + q_qrd[uid, i] <= q_0_ub[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_ub[uid]*p[uid, i]
        )
        @constraint(model,
            pq_lb_prod[uid in pq_bound_sdds_prod, i in periods],
            q[uid, i] - q_qrd[uid, i] >= q_0_lb[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_lb[uid]*p[uid, i]
        )
        @constraint(model,
            pq_lb_cons[uid in pq_bound_sdds_cons, i in periods],
            q[uid, i] - q_qru[uid, i] >= q_0_lb[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_lb[uid]*p[uid, i]
        )
        pq_eq = @constraint(
            model,
            pq_eq[uid in pq_linear_sdds, i in periods],
            q[uid, i] == q_0[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta[uid]*p[uid, i]
        )
        for uid in pq_linear_sdds
            for i in periods
                jmp.fix(q_qru[uid, i], 0.0; force = true)
                jmp.fix(q_qrd[uid, i], 0.0; force = true)
            end
        end
    else
        pq_ub = @constraint(
            model,
            pq_ub[uid in pq_bound_sdds, i in periods],
            q[uid, i] <= q_0_ub[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_ub[uid]*p[uid, i]
        )
        pq_lb = @constraint(
            model,
            pq_lb[uid in pq_bound_sdds, i in periods],
            q[uid, i] >= q_0_lb[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta_lb[uid]*p[uid, i]
        )
        pq_eq = @constraint(
            model,
            pq_eq[uid in pq_linear_sdds, i in periods],
            q[uid, i] == q_0[uid] * (
                u_on[uid, i]
                + sum(u_su[uid, t] for t in T_su_pc[uid][i], init=0)
                + sum(u_sd[uid, t] for t in T_sd_pc[uid][i], init=0)
            ) + beta[uid]*p[uid, i]
        )
        return pq_ub, pq_lb, pq_eq
    end
end

function add_cost_block_aggregation_constraint!(model, sdd_ts_lookup, periods, sdd_ids, cost_blocks)
    p_disagg = @constraint(
        model,
        p_disagg[uid=sdd_ids, i=periods],
        model[:p][uid, i]
        == sum(
            model[:cost_block_p][uid, i, j] for j in 1:length(cost_blocks[i, uid])
        ),
    )
    return p_disagg
end

function add_p_on_status_implication!(m, sdd_ts_lookup, periods, sdd_ids)
    # NOTE: I believe these constraints will be rendered invalid by the
    # p_on, p_su, and p_sd implications above.
    p_on_implies_ub = @constraint(
        m,
        p_on_implies_ub[uid in sdd_ids, i in periods],
        m[:p][uid, i] <= sdd_ts_lookup[uid]["p_ub"][i]*m[:p_on_status][uid, i],
    )
    p_on_implies_lb = @constraint(
        m,
        p_on_implies_lb[uid in sdd_ids, i in periods],
        sdd_ts_lookup[uid]["p_lb"][i]*m[:p_on_status][uid, i] <= m[:p][uid, i],
    )
    return p_on_implies_lb, p_on_implies_ub
end

function add_copperplate_power_balance!(
    m,
    periods,
    sdd_ids_producer,
    sdd_ids_consumer;
    relax = false,
    overcommitment_factor = 1.0,
)
    if relax
        @variable(m, 0 <= p_balance_slack_pos[i in periods])
        @variable(m, 0 <= p_balance_slack_neg[i in periods])
        @variable(m, 0 <= q_balance_slack_pos[i in periods])
        @variable(m, 0 <= q_balance_slack_neg[i in periods])
        copperplate_p_balance = @constraint(
            m,
            copperplate_p_balance[i in periods],
            (
                sum(m[:p][uid, i] for uid in sdd_ids_producer)
                == (
                    overcommitment_factor * sum(m[:p][uid, i] for uid in sdd_ids_consumer)
                    + p_balance_slack_pos[i] - p_balance_slack_neg[i]
                )
            )
        )
        copperplate_q_balance = @constraint(
            m,
            copperplate_q_balance[i in periods],
            (
                sum(m[:q][uid, i] for uid in sdd_ids_producer)
                == (
                    overcommitment_factor * sum(m[:q][uid, i] for uid in sdd_ids_consumer)
                    + q_balance_slack_pos[i] - q_balance_slack_neg[i]
                )
            )
        )
    else
        copperplate_p_balance = @constraint(
            m,
            copperplate_p_balance[i in periods],
            (sum(m[:p][uid, i] for uid in sdd_ids_producer)
                == overcommitment_factor * sum(m[:p][uid, i] for uid in sdd_ids_consumer)),
        )
        copperplate_q_balance = @constraint(
            m,
            copperplate_q_balance[i in periods],
            (sum(m[:q][uid, i] for uid in sdd_ids_producer)
                == overcommitment_factor * sum(m[:q][uid, i] for uid in sdd_ids_consumer)),
        )
    end
    return copperplate_p_balance, copperplate_q_balance
end

function add_device_power_cost!(m, sdd_ts_lookup, periods, sdd_ids)
    cost_blocks = Dict(
        (i, uid) => sdd_ts_lookup[uid]["cost"][i]
        for i in periods
        for uid in sdd_ids
    )
    cost_block_p = @variable(
        m,
        cost_block_p[
            uid in sdd_ids,
            i in periods,
            j in 1:length(cost_blocks[i, uid]),
        ],
        lower_bound = 0.0,
        upper_bound = cost_blocks[i, uid][j][2],
    )
    p_disagg = add_cost_block_aggregation_constraint!(
        m, sdd_ts_lookup, periods, sdd_ids, cost_blocks
    )
    device_cost = @expression(
        m,
        device_cost[uid=sdd_ids, i=periods],
        sum(cb[1] * cost_block_p[uid, i, b] for (b, cb) in enumerate(cost_blocks[i, uid]))
    )
    return cost_blocks, cost_block_p, p_disagg, device_cost
end

function add_on_su_sd_cost!(model, data)
    sdd_ids = data.sdd_ids
    periods = data.periods
    dt = data.dt
    c_on = Dict(uid => data.sdd_lookup[uid]["on_cost"] for uid in sdd_ids)
    c_su = Dict(uid => data.sdd_lookup[uid]["startup_cost"] for uid in sdd_ids)
    c_sd = Dict(uid => data.sdd_lookup[uid]["shutdown_cost"] for uid in sdd_ids)
    u_on = model[:p_on_status]
    u_su = model[:u_su]
    u_sd = model[:u_sd]

    on_cost_affexpr = JuMP.AffExpr(0.0)
    su_cost_affexpr = JuMP.AffExpr(0.0)
    sd_cost_affexpr = JuMP.AffExpr(0.0)
    for uid in sdd_ids
        for i in periods
            add_to_expression!(on_cost_affexpr, dt[i]*c_on[uid], u_on[uid, i])
            add_to_expression!(su_cost_affexpr, c_su[uid], u_su[uid, i])
            add_to_expression!(sd_cost_affexpr, c_sd[uid], u_sd[uid, i])
        end
    end
    on_cost = @expression(
        model,
        on_cost,
        on_cost_affexpr
    )
    su_cost = @expression(
        model,
        su_cost,
        su_cost_affexpr
    )
    sd_cost = @expression(
        model,
        sd_cost,
        sd_cost_affexpr
    )
    return on_cost, su_cost, sd_cost
end


"""
Return a JuMP model that implements the constraints and objective necessary to
schedule power generation/consumption with a copper-plate assumption.

`input_data` is the NamedTuple returned by `process_input_data`.

"""
function get_copperplate_scheduling_model(
    input_data::NamedTuple;
    include_balances = true,
    relax_balances = false,
    include_reserves = false,
    relax_reserves = true,
    args = nothing,
    penalize_only_reserve_shortfalls = false,
    include_objective = true,
    overcommitment_factor = 1.0,
)
    if args === nothing
        args = Dict()
    end
    print_program_info = get(args, "print_program_info", false)

    ### Data Preparation ###
    if print_program_info
        println("data preparation")
    end

    #
    # Extract useful data from processed data
    #
    (
     dt, periods, bus_lookup, bus_ids, shunt_lookup, shunt_ids, ac_line_lookup,
     ac_line_ids, twt_lookup, twt_ids, dc_line_lookup, dc_line_ids, sdd_lookup,
     sdd_ts_lookup, sdd_ids, sdd_ids_producer, sdd_ids_consumer, violation_cost,
    ) = input_data
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

    # NOTE: instead of using the UIDs as keys in the math model,
    # we may prefer to give them integer ids, some UIDs in the datasets
    # are long and confusing

    ### Mathematical Model ###
    if print_program_info
        println("model building")
    end

    model = JuMP.Model()

    p, q = add_power_variables!(model, sdd_ts_lookup, periods, sdd_ids)
    p_on_status = add_p_on_status_variables!(model, sdd_ts_lookup, periods, sdd_ids)

    #
    # Add constraints from section 2.4
    #
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

    #
    # Add constraints from section 2.6
    #
    _, p_su_ru, _, p_sd_rd = get_contiguous_power_curves(input_data)
    T_su_pc, T_sd_pc = get_inverse_power_curve_intervals(input_data)
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

    #
    # Add reserves, which are necessary for implication constraints if present
    #
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    if include_reserves
        reserve_vars = add_reserve_variables!(model, input_data)
        (p_rgu, p_rgd, p_scr, p_nsc, p_rru_on, p_rrd_on, p_rru_off,
         p_rrd_off, q_qru, q_qrd) = reserve_vars

        res_costs = _get_reserve_costs(input_data)
        (c_rgu, c_rgd, c_scr, c_nsc, c_rru_on, c_rrd_on, c_rru_off,
         c_rrd_off, c_qru, c_qrd) = res_costs

        res_shortfalls = add_reserve_shortfall_variables!(model, input_data)
        (p_rgu_slack, p_rgd_slack, p_scr_slack, p_nsc_slack, p_rru_slack,
         p_rrd_slack, q_qru_slack, q_qrd_slack) = res_shortfalls

        res_penalties = _get_reserve_shortfall_penalties(input_data)
        z_rgu, z_rgd, z_scr, z_nsc, z_rru, z_rrd, z_qru, z_qrd = res_penalties

        add_reserve_balances!(
            model, input_data, model[:p]; relax = relax_reserves
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
    ###

    #
    # These constraints need to be updated if we are including reserves
    #
    add_on_su_sd_implication_constraints!(
        model, input_data; include_reserves = include_reserves
    )
    add_reactive_power_implication_constraints!(
        model, input_data, T_su_pc, T_sd_pc; include_reserves = include_reserves
    )
    add_real_reactive_linking_constraints!(
        model, input_data, T_su_pc, T_sd_pc; include_reserves = include_reserves
    )
    ###
    if include_balances
        copperplate_p_balance, copperplate_q_balance = add_copperplate_power_balance!(
            model, periods, sdd_ids_producer, sdd_ids_consumer;
            relax = relax_balances,
            overcommitment_factor = overcommitment_factor,
        )
    end

    ### Mathematical Model - Objective ###

    # TODO: Change name of p_disagg (cost_block_aggregation?)
    cost_blocks, cost_block_p, p_disagg, device_cost = add_device_power_cost!(
        model, sdd_ts_lookup, periods, sdd_ids 
    )
    cost_by_period = @expression(
        model,
        cost_by_period[i=periods],
        sum(device_cost[uid, i] for uid in sdd_ids_consumer) -
        sum(device_cost[uid, i] for uid in sdd_ids_producer)
    )

    on_cost, su_cost, sd_cost = add_on_su_sd_cost!(model, input_data)

    if penalize_only_reserve_shortfalls
        if !include_reserves
            throw(Exception)
        end
        # It is imprecise to call this a "cost expression", as it also
        # includes some penalties on the copper-plate balance violations
        cost_expr = @expression(model, cost_expr, 0.0)
        reserve_cost_expr = @expression(model, reserve_cost_expr, 0.0)
    elseif relax_balances && include_balances
        e_vio_cost = violation_cost["e_vio_cost"]
        p_balance_slack_pos = model[:p_balance_slack_pos]
        p_balance_slack_neg = model[:p_balance_slack_neg]
        q_balance_slack_pos = model[:q_balance_slack_pos]
        q_balance_slack_neg = model[:q_balance_slack_neg]
        cost_expr = @expression(model,
            cost_expr,
            sum(dt[i]*cost_by_period[i] for i in periods)
            - on_cost - su_cost - sd_cost
            - e_vio_cost*sum(
                p_balance_slack_pos[i] + p_balance_slack_neg[i]
                + q_balance_slack_pos[i] + q_balance_slack_neg[i]
                for i in periods
            )
        )
        reserve_cost_expr = @expression(model,
            reserve_cost_expr,
            # Reserve cost expressions
            # If reserves are not included, these are zero
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
    else
        cost_expr = @expression(model,
            cost_expr,
            sum(dt[i]*cost_by_period[i] for i in periods)
            - on_cost - su_cost - sd_cost
        )
        reserve_cost_expr = @expression(model,
            reserve_cost_expr,
            # Reserve cost expressions
            # If reserves are not included, these are zero
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
    end

    #
    # Market surplus objective
    #
    if include_objective
        @objective(model,
            Max,
            # These are zero if we are only penalizing reserve shortfall
            cost_expr + reserve_cost_expr
            # Reserve shortfall expressions
            # If reserves are not included, these are zero
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
        )
    end

    return model
end


function solve_scheduling_model(
    model::jmp.Model;
    solver = HiGHS.Optimizer,
    relax_integrality::Bool = false,

    # This parameter is deprecated. This should be specified as an attribute
    # on the provided optimizer.
    #accept_first_feasible::Bool = false,
)
    if relax_integrality
        jmp.relax_integrality(model)
    end
    jmp.set_optimizer(model, solver)
    return optimize!(model)
end


function extract_data_from_scheduling_model(
    input_data::NamedTuple,
    model::JuMP.Model;
    include_reserves = false,
)
    periods = input_data.periods
    ids = input_data.sdd_ids
    u_on = model[:p_on_status]
    p = model[:p]
    q = model[:q]
    on_status = Dict(uid => [value(u_on[uid, i]) for i in periods] for uid in ids)
    real_power = Dict(uid => [value(p[uid, i]) for i in periods] for uid in ids)
    reactive_power = Dict(uid => [value(q[uid, i]) for i in periods] for uid in ids)

    if include_reserves
        p_rgu = model[:p_rgu]
        p_rgd = model[:p_rgd]
        p_scr = model[:p_scr]
        p_nsc = model[:p_nsc]
        p_rru_on = model[:p_rru_on]
        p_rrd_on = model[:p_rrd_on]
        p_rru_off = model[:p_rru_off]
        p_rrd_off = model[:p_rrd_off]
        q_qru = model[:q_qru]
        q_qrd = model[:q_qrd]
        p_rgu = Dict(uid => [value(p_rgu[uid, i]) for i in periods] for uid in ids)
        p_rgd = Dict(uid => [value(p_rgd[uid, i]) for i in periods] for uid in ids)
        p_scr = Dict(uid => [value(p_scr[uid, i]) for i in periods] for uid in ids)
        p_nsc = Dict(uid => [value(p_nsc[uid, i]) for i in periods] for uid in ids)
        p_rru_on = Dict(uid => [value(p_rru_on[uid, i]) for i in periods] for uid in ids)
        p_rru_off = Dict(uid => [value(p_rru_off[uid, i]) for i in periods] for uid in ids)
        p_rrd_on = Dict(uid => [value(p_rrd_on[uid, i]) for i in periods] for uid in ids)
        p_rrd_off = Dict(uid => [value(p_rrd_off[uid, i]) for i in periods] for uid in ids)
        q_qru = Dict(uid => [value(q_qru[uid, i]) for i in periods] for uid in ids)
        q_qrd = Dict(uid => [value(q_qrd[uid, i]) for i in periods] for uid in ids)
        return (
            on_status = on_status,
            real_power = real_power,
            reactive_power = reactive_power,
            p_rgu = p_rgu,
            p_rgd = p_rgd,
            p_scr = p_scr,
            p_nsc = p_nsc,
            p_rru_on = p_rru_on,
            p_rru_off = p_rru_off,
            p_rrd_on = p_rrd_on,
            p_rrd_off = p_rrd_off,
            q_qru = q_qru,
            q_qrd = q_qrd,
        )
    else
        return (
            on_status = on_status,
            real_power = real_power,
            reactive_power = reactive_power,
        )
    end
end


function extract_reserves_from_scheduling_model(
    input_data::NamedTuple,
    model::JuMP.Model,
)
    ids = input_data.sdd_ids
    periods = input_data.periods
    p_rgu = model[:p_rgu]
    p_rgd = model[:p_rgd]
    p_scr = model[:p_scr]
    p_nsc = model[:p_nsc]
    p_rru_on = model[:p_rru_on]
    p_rrd_on = model[:p_rrd_on]
    p_rru_off = model[:p_rru_off]
    p_rrd_off = model[:p_rrd_off]
    q_qru = model[:q_qru]
    q_qrd = model[:q_qrd]
    p_rgu = Dict(uid => [value(p_rgu[uid, i]) for i in periods] for uid in ids)
    p_rgd = Dict(uid => [value(p_rgd[uid, i]) for i in periods] for uid in ids)
    p_scr = Dict(uid => [value(p_scr[uid, i]) for i in periods] for uid in ids)
    p_nsc = Dict(uid => [value(p_nsc[uid, i]) for i in periods] for uid in ids)
    p_rru_on = Dict(uid => [value(p_rru_on[uid, i]) for i in periods] for uid in ids)
    p_rru_off = Dict(uid => [value(p_rru_off[uid, i]) for i in periods] for uid in ids)
    p_rrd_on = Dict(uid => [value(p_rrd_on[uid, i]) for i in periods] for uid in ids)
    p_rrd_off = Dict(uid => [value(p_rrd_off[uid, i]) for i in periods] for uid in ids)
    q_qru = Dict(uid => [value(q_qru[uid, i]) for i in periods] for uid in ids)
    q_qrd = Dict(uid => [value(q_qrd[uid, i]) for i in periods] for uid in ids)
    return (
        p_rgu = p_rgu,
        p_rgd = p_rgd,
        p_scr = p_scr,
        p_nsc = p_nsc,
        p_rru_on = p_rru_on,
        p_rru_off = p_rru_off,
        p_rrd_on = p_rrd_on,
        p_rrd_off = p_rrd_off,
        q_qru = q_qru,
        q_qrd = q_qrd,
    )    
end

function compute_device_reserve_values(
    input_data::NamedTuple,
    schedule_data::NamedTuple,
)
    # "Reserve value" is defined as:
    # (shortfall penalty - reserve cost) * amount of reserve
    # summed over time. Each device can have a reserve value in the positive
    # or negative direction.
    #
    # More specifically, this is the penalty associated with changing this
    # device's power level. We do not compute this for offline (or ramping)
    # devices, as the power level does not change during OPF solves. Should
    # we even include these devices in the dict?
    pos_reserve_value = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    neg_reserve_value = Dict(uid => 0.0 for uid in input_data.sdd_ids)

    costs = _get_reserve_costs(input_data)
    penalties = _get_reserve_shortfall_penalties(input_data)
    zone_data = _get_sdd_in_zones(input_data)

    sdd_rgu_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    sdd_rgd_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    sdd_scr_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    sdd_nsc_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    sdd_rrd_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)
    sdd_rru_penalty = Dict(uid => 0.0 for uid in input_data.sdd_ids)

    for azr_id in input_data.azr_ids
        for sdd_id in zone_data.sdd_in_azone[azr_id]
            # What if the sdd is in more than one active reserve zone?
            # Pretty sure this shouldn't happen
            sdd_rgu_penalty[sdd_id] = penalties.z_rgu[azr_id]
            sdd_rgd_penalty[sdd_id] = penalties.z_rgd[azr_id]
            sdd_scr_penalty[sdd_id] = penalties.z_scr[azr_id]
            sdd_nsc_penalty[sdd_id] = penalties.z_nsc[azr_id]
            sdd_rru_penalty[sdd_id] = penalties.z_rru[azr_id]
            sdd_rrd_penalty[sdd_id] = penalties.z_rrd[azr_id]
        end
    end

    for uid in input_data.sdd_ids
        up_reserve_value = sum(
            input_data.dt[i] * (
                (
                    sdd_rgu_penalty[uid] + sdd_scr_penalty[uid]
                    + sdd_nsc_penalty[uid] - costs.c_rgu[uid][i]
                ) * schedule_data.p_rgu[uid][i]
                + (
                    sdd_scr_penalty[uid] + sdd_nsc_penalty[uid]
                    - costs.c_scr[uid][i]
                ) * schedule_data.p_scr[uid][i]
                + (
                    sdd_rru_penalty[uid] - costs.c_rru_on[uid][i]
                ) * schedule_data.p_rru_on[uid][i]
            )
            for i in input_data.periods
        )
        down_reserve_value = sum(
            input_data.dt[i] * (
                (
                    sdd_rgd_penalty[uid] - costs.c_rgd[uid][i]
                ) * schedule_data.p_rgd[uid][i]
                + (
                    sdd_rrd_penalty[uid] - costs.c_rrd_on[uid][i]
                ) * schedule_data.p_rrd_on[uid][i]
            )
            for i in input_data.periods
        )
        if input_data.sdd_lookup[uid]["device_type"] == "consumer"
            # Positive change violates rgd and rrd_on
            # Negative change violates rgu, scr, and rru_on
            pos_reserve_value[uid] = down_reserve_value
            neg_reserve_value[uid] = up_reserve_value
        elseif input_data.sdd_lookup[uid]["device_type"] == "producer"
            # Positive change violates rgu, scr, and rru_on
            # Negative change violates rgd, rrd_on
            pos_reserve_value[uid] = up_reserve_value
            neg_reserve_value[uid] = down_reserve_value
        else
            throw(Exception)
        end
    end
    reserve_value = Dict(
        uid => max(pos_reserve_value[uid], neg_reserve_value[uid], 0.0)
        for uid in input_data.sdd_ids
    )
    #return reserve_value
    # Reserve violations associated with positive and negative changes
    return pos_reserve_value, neg_reserve_value
end

function get_point_from_solution_data(data::Dict, model::JuMP.Model, input_data::NamedTuple)
    ts_key = "time_series_output"
    sdd_key = "simple_dispatchable_device"
    sdd_ids = input_data.sdd_ids
    periods = input_data.periods
    tsout = data[ts_key]
    sdds = tsout[sdd_key]
    point = Dict{JuMP.VariableRef, Float64}()

    p_on = model[:p_on]
    p = model[:p]
    q = model[:q]
    u = model[:p_on_status]
    # TODO: Compute u_su, u_sd, p_su, and p_sd
    for sdd in sdds
        uid = sdd["uid"]
        for i in periods
            point[p_on[uid, i]] = sdd["p_on"][i]
            point[q[uid, i]] = sdd["q"][i]
            point[u[uid, i]] = sdd["on_status"][i]

            # Hack to get around computing p, which depends on p_su and p_sd,
            # which depend on power curves.
            if Bool(sdd["on_status"][i])
                point[p[uid, i]] = sdd["p_on"][i]
            end
        end
    end
    return point
end

function warmstart!(
    model::JuMP.Model,
    input_data::NamedTuple,
    schedule_data::NamedTuple,
)
    for i in input_data.periods
        for uid in input_data.sdd_ids
            p = model[:p]
            u = model[:p_on_status]
            set_start_value(p[uid, i], schedule_data.real_power[uid][i])
            set_start_value(u[uid, i], schedule_data.on_status[uid][i])
        end
    end
end

#end # SchedulingModel
