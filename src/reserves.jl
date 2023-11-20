using JuMP
using HiGHS

include("power_curves.jl")

"""
Usage:

reserve_ubs = _get_max_reserves(input_data)
(p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
 p_rru_off_max, p_rrd_off_max) = reserve_ubs

"""
function _get_max_reserves(input_data::NamedTuple)
    sdd_ids = input_data.sdd_ids
    sdd_lookup = input_data.sdd_lookup
    p_rgu_key = "p_reg_res_up_ub"
    p_rgd_key = "p_reg_res_down_ub"
    p_scr_key = "p_syn_res_ub"
    p_nsc_key = "p_nsyn_res_ub"
    p_rru_on_key = "p_ramp_res_up_online_ub"
    p_rrd_on_key = "p_ramp_res_down_online_ub"
    p_rru_off_key = "p_ramp_res_up_offline_ub"
    p_rrd_off_key = "p_ramp_res_down_offline_ub"
    p_rgu_max = Dict(uid => sdd_lookup[uid][p_rgu_key] for uid in sdd_ids)
    p_rgd_max = Dict(uid => sdd_lookup[uid][p_rgd_key] for uid in sdd_ids)
    p_scr_max = Dict(uid => sdd_lookup[uid][p_scr_key] for uid in sdd_ids)
    p_nsc_max = Dict(uid => sdd_lookup[uid][p_nsc_key] for uid in sdd_ids)
    p_rru_on_max = Dict(uid => sdd_lookup[uid][p_rru_on_key] for uid in sdd_ids)
    p_rrd_on_max = Dict(uid => sdd_lookup[uid][p_rrd_on_key] for uid in sdd_ids)
    p_rru_off_max = Dict(uid => sdd_lookup[uid][p_rru_off_key] for uid in sdd_ids)
    p_rrd_off_max = Dict(uid => sdd_lookup[uid][p_rrd_off_key] for uid in sdd_ids)
    return (
        p_rgu_max = p_rgu_max,
        p_rgd_max = p_rgd_max,
        p_scr_max = p_scr_max,
        p_nsc_max = p_nsc_max,
        p_rru_on_max = p_rru_on_max,
        p_rrd_on_max = p_rrd_on_max,
        p_rru_off_max = p_rru_off_max,
        p_rrd_off_max = p_rrd_off_max,
    )
end

function _get_reserve_costs(input_data::NamedTuple)
    sdd_ids = input_data.sdd_ids
    sdd_ts_lookup = input_data.sdd_ts_lookup
    c_rgu_key = "p_reg_res_up_cost"
    c_rgd_key = "p_reg_res_down_cost"
    c_scr_key = "p_syn_res_cost"
    c_nsc_key = "p_nsyn_res_cost"
    c_rru_on_key = "p_ramp_res_up_online_cost"
    c_rrd_on_key = "p_ramp_res_down_online_cost"
    c_rru_off_key = "p_ramp_res_up_offline_cost"
    c_rrd_off_key = "p_ramp_res_down_offline_cost"
    c_qru_key = "q_res_up_cost"
    c_qrd_key = "q_res_down_cost"
    c_rgu = Dict(uid => sdd_ts_lookup[uid][c_rgu_key] for uid in sdd_ids)
    c_rgd = Dict(uid => sdd_ts_lookup[uid][c_rgd_key] for uid in sdd_ids)
    c_scr = Dict(uid => sdd_ts_lookup[uid][c_scr_key] for uid in sdd_ids)
    c_nsc = Dict(uid => sdd_ts_lookup[uid][c_nsc_key] for uid in sdd_ids)
    c_rru_on = Dict(uid => sdd_ts_lookup[uid][c_rru_on_key] for uid in sdd_ids)
    c_rrd_on = Dict(uid => sdd_ts_lookup[uid][c_rrd_on_key] for uid in sdd_ids)
    c_rru_off = Dict(uid => sdd_ts_lookup[uid][c_rru_off_key] for uid in sdd_ids)
    c_rrd_off = Dict(uid => sdd_ts_lookup[uid][c_rrd_off_key] for uid in sdd_ids)
    c_qru = Dict(uid => sdd_ts_lookup[uid][c_qru_key] for uid in sdd_ids)
    c_qrd = Dict(uid => sdd_ts_lookup[uid][c_qrd_key] for uid in sdd_ids)
    return (
        c_rgu = c_rgu,
        c_rgd = c_rgd,
        c_scr = c_scr,
        c_nsc = c_nsc,
        c_rru_on = c_rru_on,
        c_rrd_on = c_rrd_on,
        c_rru_off = c_rru_off,
        c_rrd_off = c_rrd_off,
        c_qru = c_qru,
        c_qrd = c_qrd,
    )
end

function _get_sdd_in_zones(input_data::NamedTuple)
    cons_ids = input_data.sdd_ids_consumer
    prod_ids = input_data.sdd_ids_producer
    cons_set = Set(cons_ids)
    prod_set = Set(prod_ids)
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    sdd_in_azone = Dict(uid => [] for uid in azr_ids)
    sdd_in_rzone = Dict(uid => [] for uid in rzr_ids)
    bus_sdds = Dict(uid => [] for uid in input_data.bus_ids)
    for uid in input_data.sdd_ids
        bus_id = input_data.sdd_lookup[uid]["bus"]
        push!(bus_sdds[bus_id], uid)
    end
    for bus_id in input_data.bus_ids
        for zone_id in input_data.bus_lookup[bus_id]["active_reserve_uids"]
            for sdd_id in bus_sdds[bus_id]
                push!(sdd_in_azone[zone_id], sdd_id)
            end
        end
        for zone_id in input_data.bus_lookup[bus_id]["reactive_reserve_uids"]
            for sdd_id in bus_sdds[bus_id]
                push!(sdd_in_rzone[zone_id], sdd_id)
            end
        end
    end
    c_sdd_in_azone = Dict(
        azr_id => [sdd_id for sdd_id in sdd_in_azone[azr_id] if sdd_id in cons_set]
        for azr_id in azr_ids
    )
    p_sdd_in_azone = Dict(
        azr_id => [sdd_id for sdd_id in sdd_in_azone[azr_id] if sdd_id in prod_set]
        for azr_id in azr_ids
    )
    return (
        sdd_in_azone = sdd_in_azone,
        sdd_in_rzone = sdd_in_rzone,
        c_sdd_in_azone = c_sdd_in_azone,
        p_sdd_in_azone = p_sdd_in_azone,
    )
    # Usage:
    # zone_data = _get_sdd_in_zones(input_data)
    # sdd_in_azone, sdd_in_rzone, c_sdd_in_azone, p_sdd_in_azone = zone_data
end

function _get_reserve_shortfall_penalties(input_data::NamedTuple)
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    azr_lookup = input_data.azr_lookup
    rzr_lookup = input_data.rzr_lookup
    # Keys for the zonal violation penalties
    z_rgu_key = "REG_UP_vio_cost"
    z_rgd_key = "REG_DOWN_vio_cost"
    z_scr_key = "SYN_vio_cost"
    z_nsc_key = "NSYN_vio_cost"
    z_rru_key = "RAMPING_RESERVE_UP_vio_cost"
    z_rrd_key = "RAMPING_RESERVE_DOWN_vio_cost"
    z_qru_key = "REACT_UP_vio_cost"
    z_qrd_key = "REACT_DOWN_vio_cost"
    z_rgu = Dict(uid => azr_lookup[uid][z_rgu_key] for uid in azr_ids)
    z_rgd = Dict(uid => azr_lookup[uid][z_rgd_key] for uid in azr_ids)
    z_scr = Dict(uid => azr_lookup[uid][z_scr_key] for uid in azr_ids)
    z_nsc = Dict(uid => azr_lookup[uid][z_nsc_key] for uid in azr_ids)
    z_rru = Dict(uid => azr_lookup[uid][z_rru_key] for uid in azr_ids)
    z_rrd = Dict(uid => azr_lookup[uid][z_rrd_key] for uid in azr_ids)
    z_qru = Dict(uid => rzr_lookup[uid][z_qru_key] for uid in rzr_ids)
    z_qrd = Dict(uid => rzr_lookup[uid][z_qrd_key] for uid in rzr_ids)
    return (
        z_rgu = z_rgu,
        z_rgd = z_rgd,
        z_scr = z_scr,
        z_nsc = z_nsc,
        z_rru = z_rru,
        z_rrd = z_rrd,
        z_qru = z_qru,
        z_qrd = z_qrd,
    )
    # Usage:
    # res_penalties = _get_reserve_shortfall_penalties(input_data)
    # z_rgu, z_rgd, z_scr, z_nsc, z_rru, z_rrd, z_qru, z_qrd = res_penalties
end

function _get_reserve_requirements(
    input_data::NamedTuple, zone_data::NamedTuple, interval::Int64, p_sdd::Dict
)
    i = interval
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    azr_ts_lookup = input_data.azr_ts_lookup
    rzr_ts_lookup = input_data.rzr_ts_lookup
    azr_lookup = input_data.azr_lookup
    c_sdd_in_azone = zone_data.c_sdd_in_azone
    p_sdd_in_azone = zone_data.p_sdd_in_azone
    rru_key = "RAMPING_RESERVE_UP"
    rrd_key = "RAMPING_RESERVE_DOWN"
    qru_key = "REACT_UP"
    qrd_key = "REACT_DOWN"
    p_rru_req = Dict(uid => azr_ts_lookup[uid][rru_key][i] for uid in azr_ids)
    p_rrd_req = Dict(uid => azr_ts_lookup[uid][rrd_key][i] for uid in azr_ids)
    q_qru_req = Dict(uid => rzr_ts_lookup[uid][qru_key][i] for uid in rzr_ids)
    q_qrd_req = Dict(uid => rzr_ts_lookup[uid][qrd_key][i] for uid in rzr_ids)

    rgu_factor = Dict(uid => azr_lookup[uid]["REG_UP"] for uid in azr_ids)
    rgd_factor = Dict(uid => azr_lookup[uid]["REG_DOWN"] for uid in azr_ids)
    scr_factor = Dict(uid => azr_lookup[uid]["SYN"] for uid in azr_ids)
    nsc_factor = Dict(uid => azr_lookup[uid]["NSYN"] for uid in azr_ids)

    # An interval was provided. We assume p_sdd maps uid -> float.

    # Note that these endogenous reserve requirements are calculated
    # using p_sdd = p_on + p_su + p_sd, not just p_on
    p_rgu_req = Dict(
        azr_id => rgu_factor[azr_id]*sum(
            p_sdd[sdd_id] for sdd_id in c_sdd_in_azone[azr_id]; init = 0
        ) for azr_id in azr_ids
    )
    p_rgd_req = Dict(
        azr_id => rgd_factor[azr_id]*sum(
            p_sdd[sdd_id] for sdd_id in c_sdd_in_azone[azr_id]; init = 0
        ) for azr_id in azr_ids
    )
    p_scr_req = Dict(
        azr_id => scr_factor[azr_id]*maximum(
            p_sdd[sdd_id] for sdd_id in p_sdd_in_azone[azr_id]; init = 0
        ) for azr_id in azr_ids
    )
    p_nsc_req = Dict(
        azr_id => nsc_factor[azr_id]*maximum(
            p_sdd[sdd_id] for sdd_id in p_sdd_in_azone[azr_id]; init = 0
        ) for azr_id in azr_ids
    )
    return (
        # Exogenous requirements
        p_rru_req = p_rru_req,
        p_rrd_req = p_rrd_req,
        q_qru_req = q_qru_req,
        q_qrd_req = q_qrd_req,
        # Endogenous requirements
        p_rgu_req = p_rgu_req,
        p_rgd_req = p_rgd_req,
        p_scr_req = p_scr_req,
        p_nsc_req = p_nsc_req,
    )
    # Usage:
    # res_req = _get_reserve_requirements(input_data, zone_data, i, p_sdd)
    # (p_rru_req, p_rrd_req, q_qru_req, q_qrd_req,
    # p_rgu_req, p_rgd_req, p_scr_req, p_nsc_req) = res_req
end

# This is the "time-indexed" implementation. No interval is provided, and p_sdd
# is an array of JuMP VariableRef.
function _get_reserve_requirements(
    input_data::NamedTuple, zone_data::NamedTuple, p_sdd::JuMP.Containers.DenseAxisArray
)
    periods = input_data.periods
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    azr_ts_lookup = input_data.azr_ts_lookup
    rzr_ts_lookup = input_data.rzr_ts_lookup
    azr_lookup = input_data.azr_lookup
    c_sdd_in_azone = zone_data.c_sdd_in_azone
    p_sdd_in_azone = zone_data.p_sdd_in_azone
    rru_key = "RAMPING_RESERVE_UP"
    rrd_key = "RAMPING_RESERVE_DOWN"
    qru_key = "REACT_UP"
    qrd_key = "REACT_DOWN"
    p_rru_req = Dict(uid => azr_ts_lookup[uid][rru_key] for uid in azr_ids)
    p_rrd_req = Dict(uid => azr_ts_lookup[uid][rrd_key] for uid in azr_ids)
    q_qru_req = Dict(uid => rzr_ts_lookup[uid][qru_key] for uid in rzr_ids)
    q_qrd_req = Dict(uid => rzr_ts_lookup[uid][qrd_key] for uid in rzr_ids)

    rgu_factor = Dict(uid => azr_lookup[uid]["REG_UP"] for uid in azr_ids)
    rgd_factor = Dict(uid => azr_lookup[uid]["REG_DOWN"] for uid in azr_ids)
    scr_factor = Dict(uid => azr_lookup[uid]["SYN"] for uid in azr_ids)
    nsc_factor = Dict(uid => azr_lookup[uid]["NSYN"] for uid in azr_ids)

    # Note that these endogenous reserve requirements are calculated
    # using p_sdd = p_on + p_su + p_sd, not just p_on
    #
    # NOTE: These dicts contain vectors of JuMP expressions. Should the
    # expressions be created with add_to_expression! ? Probably.
    sum_cons_in_zone = Dict(
        azr_id => [JuMP.AffExpr(0.0) for i in periods] for azr_id in azr_ids
    )
    #p_rgu_req = Dict(
    #    azr_id => [JuMP.AffExpr(0.0) for i in periods] for azr_id in azr_ids
    #)
    #p_rgd_req = Dict(
    #    azr_id => [JuMP.AffExpr(0.0) for i in periods] for azr_id in azr_ids
    #)
    for azr_id in azr_ids
        for sdd_id in c_sdd_in_azone[azr_id]
            for i in periods
                JuMP.add_to_expression!(
                    sum_cons_in_zone[azr_id][i],
                    p_sdd[sdd_id, i],
                )
            end
        end
    end
    p_rgu_req = Dict(
        azr_id => [
            # Note that this constructs a new expression
            rgu_factor[azr_id]*sum_cons_in_zone[azr_id][i] for i in periods
        ]
        for azr_id in azr_ids
    )
    p_rgd_req = Dict(
        azr_id => [
            # Note that this constructs a new expression
            rgd_factor[azr_id]*sum_cons_in_zone[azr_id][i] for i in periods
        ]
        for azr_id in azr_ids
    )
    #p_rgu_req = Dict(
    #    azr_id => [
    #        rgu_factor[azr_id]*sum(
    #            p_sdd[sdd_id, i] for sdd_id in c_sdd_in_azone[azr_id]; init = 0
    #        )
    #        for i in periods
    #    ]
    #    for azr_id in azr_ids
    #)
    #p_rgd_req = Dict(
    #    azr_id => [
    #        rgd_factor[azr_id]*sum(
    #            p_sdd[sdd_id, i] for sdd_id in c_sdd_in_azone[azr_id]; init = 0
    #        )
    #        for i in periods
    #    ]
    #    for azr_id in azr_ids
    #)
    #
    # With p_sdd a JuMP variable, these maxes will have to be reformulated.
    # These will need to be dicts: azr_id -> sdd_id -> [scr_factor*p_sdd[1],...]
    #
    p_scr_req = Dict(
        azr_id => Dict(
            sdd_id => [scr_factor[azr_id]*p_sdd[sdd_id, i] for i in periods]
            for sdd_id in p_sdd_in_azone[azr_id]
        ) for azr_id in azr_ids
    )
    p_nsc_req = Dict(
        azr_id => Dict(
            sdd_id => [nsc_factor[azr_id]*p_sdd[sdd_id, i] for i in periods]
            for sdd_id in p_sdd_in_azone[azr_id]
        ) for azr_id in azr_ids
    )
    return (
        # Exogenous requirements
        p_rru_req = p_rru_req,
        p_rrd_req = p_rrd_req,
        q_qru_req = q_qru_req,
        q_qrd_req = q_qrd_req,
        # Endogenous requirements
        p_rgu_req = p_rgu_req,
        p_rgd_req = p_rgd_req,
        p_scr_req = p_scr_req,
        p_nsc_req = p_nsc_req,
    )
    # Usage:
    # res_req = _get_reserve_requirements(input_data, zone_data, i, p_sdd)
    # (p_rru_req, p_rrd_req, q_qru_req, q_qrd_req,
    # p_rgu_req, p_rgd_req, p_scr_req, p_nsc_req) = res_req
end

function add_reserve_variables!(model::JuMP.Model, input_data::NamedTuple)
    sdd_ids = input_data.sdd_ids
    sdd_ts_lookup = input_data.sdd_ts_lookup
    periods = input_data.periods

    reserve_ubs = _get_max_reserves(input_data)
    (p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
     p_rru_off_max, p_rrd_off_max) = reserve_ubs
    q_res_max = Dict(
        # q reserves do not have an "absolute limit", so we use ub-lb, which
        # is also a valid bound.
        uid => [
            sdd_ts_lookup[uid]["q_ub"][i] - sdd_ts_lookup[uid]["q_lb"][i]
            for i in periods
        ]
        for uid in sdd_ids
    )

    p_rgu = @variable(model, 0 <= p_rgu[uid in sdd_ids, periods] <= p_rgu_max[uid])
    p_rgd = @variable(model, 0 <= p_rgd[uid in sdd_ids, periods] <= p_rgd_max[uid])
    p_scr = @variable(model, 0 <= p_scr[uid in sdd_ids, periods] <= p_scr_max[uid])
    p_nsc = @variable(model, 0 <= p_nsc[uid in sdd_ids, periods] <= p_nsc_max[uid])
    p_rru_on = @variable(model, 0 <= p_rru_on[uid in sdd_ids, periods] <= p_rru_on_max[uid])
    p_rrd_on = @variable(model, 0 <= p_rrd_on[uid in sdd_ids, periods] <= p_rrd_on_max[uid])
    p_rru_off = @variable(model, 0 <= p_rru_off[uid in sdd_ids, periods] <= p_rru_off_max[uid])
    p_rrd_off = @variable(model, 0 <= p_rrd_off[uid in sdd_ids, periods] <= p_rrd_off_max[uid])
    q_qru = @variable(model, 0 <= q_qru[uid in sdd_ids, i in periods] <= q_res_max[uid][i])
    q_qrd = @variable(model, 0 <= q_qrd[uid in sdd_ids, i in periods] <= q_res_max[uid][i])
    return (
        p_rgu=p_rgu,
        p_rgd=p_rgd,
        p_scr=p_scr,
        p_nsc=p_nsc,
        p_rru_on=p_rru_on,
        p_rrd_on=p_rrd_on,
        p_rru_off=p_rru_off,
        p_rrd_off=p_rrd_off,
        q_qru=q_qru,
        q_qrd=q_qrd,
    )
end

function add_reserve_variables!(
    model::JuMP.Model, input_data::NamedTuple, interval::Int64
)
    sdd_ids = input_data.sdd_ids
    sdd_ts_lookup = input_data.sdd_ts_lookup
    i = interval

    reserve_ubs = _get_max_reserves(input_data)
    (p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
     p_rru_off_max, p_rrd_off_max) = reserve_ubs
    q_res_max = Dict(
        # q reserves do not have an "absolute limit", so we use ub-lb, which
        # is also a valid bound.
        uid => (sdd_ts_lookup[uid]["q_ub"][i] - sdd_ts_lookup[uid]["q_lb"][i])
        for uid in sdd_ids
    )

    # Note e.g. that p_rgu <= min(p_rgu_max, p_scr_max, p_rru_on_max)
    # TODO: Tighten these bounds using the min of all applicable bounds.
    p_rgu = @variable(model, 0 <= p_rgu[uid in sdd_ids] <= p_rgu_max[uid])
    p_rgd = @variable(model, 0 <= p_rgd[uid in sdd_ids] <= p_rgd_max[uid])
    p_scr = @variable(model, 0 <= p_scr[uid in sdd_ids] <= p_scr_max[uid])
    p_nsc = @variable(model, 0 <= p_nsc[uid in sdd_ids] <= p_nsc_max[uid])
    p_rru_on = @variable(model, 0 <= p_rru_on[uid in sdd_ids] <= p_rru_on_max[uid])
    p_rrd_on = @variable(model, 0 <= p_rrd_on[uid in sdd_ids] <= p_rrd_on_max[uid])
    p_rru_off = @variable(model, 0 <= p_rru_off[uid in sdd_ids] <= p_rru_off_max[uid])
    p_rrd_off = @variable(model, 0 <= p_rrd_off[uid in sdd_ids] <= p_rrd_off_max[uid])
    q_qru = @variable(model, 0 <= q_qru[uid in sdd_ids] <= q_res_max[uid])
    q_qrd = @variable(model, 0 <= q_qrd[uid in sdd_ids] <= q_res_max[uid])
    return (
        p_rgu=p_rgu,
        p_rgd=p_rgd,
        p_scr=p_scr,
        p_nsc=p_nsc,
        p_rru_on=p_rru_on,
        p_rrd_on=p_rrd_on,
        p_rru_off=p_rru_off,
        p_rrd_off=p_rrd_off,
        q_qru=q_qru,
        q_qrd=q_qrd,
    )
end

function add_reserve_shortfall_variables!(model::Model, input_data::NamedTuple)
    periods = input_data.periods
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    # No interval was provided. Index by UIDs and periods
    p_rgu_slack = @variable(model, 0 <= p_rgu_slack[uid in azr_ids, periods])
    p_rgd_slack = @variable(model, 0 <= p_rgd_slack[uid in azr_ids, periods])
    p_scr_slack = @variable(model, 0 <= p_scr_slack[uid in azr_ids, periods])
    p_nsc_slack = @variable(model, 0 <= p_nsc_slack[uid in azr_ids, periods])
    p_rru_slack = @variable(model, 0 <= p_rru_slack[uid in azr_ids, periods])
    p_rrd_slack = @variable(model, 0 <= p_rrd_slack[uid in azr_ids, periods])
    q_qru_slack = @variable(model, 0 <= q_qru_slack[uid in rzr_ids, periods])
    q_qrd_slack = @variable(model, 0 <= q_qrd_slack[uid in rzr_ids, periods])
    return (
        p_rgu_slack,
        p_rgd_slack,
        p_scr_slack,
        p_nsc_slack,
        p_rru_slack,
        p_rrd_slack,
        q_qru_slack,
        q_qrd_slack,
    )
end

function add_reserve_shortfall_variables!(
    model::Model, input_data::NamedTuple, interval::Int64
)
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    # An interval was provided. Index only by UIDs
    p_rgu_slack = @variable(model, 0 <= p_rgu_slack[uid in azr_ids])
    p_rgd_slack = @variable(model, 0 <= p_rgd_slack[uid in azr_ids])
    p_scr_slack = @variable(model, 0 <= p_scr_slack[uid in azr_ids])
    p_nsc_slack = @variable(model, 0 <= p_nsc_slack[uid in azr_ids])
    p_rru_slack = @variable(model, 0 <= p_rru_slack[uid in azr_ids])
    p_rrd_slack = @variable(model, 0 <= p_rrd_slack[uid in azr_ids])
    q_qru_slack = @variable(model, 0 <= q_qru_slack[uid in rzr_ids])
    q_qrd_slack = @variable(model, 0 <= q_qrd_slack[uid in rzr_ids])
    return (
        p_rgu_slack,
        p_rgd_slack,
        p_scr_slack,
        p_nsc_slack,
        p_rru_slack,
        p_rrd_slack,
        q_qru_slack,
        q_qrd_slack,
    )
    # Usage:
    # res_shortfalls = add_reserve_shortfall_variables!(...)
    # (p_rgu_slack, p_rgd_slack, p_scr_slack, p_nsc_slack, p_rru_slack,
    #  p_rrd_slack, q_qru_slack, q_qrd_slack) = res_shortfalls
end

function add_reserve_balances!(
    model::Model, input_data::NamedTuple, interval::Int64, p_sdd::Dict
)
    i = interval
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    zone_data = _get_sdd_in_zones(input_data)
    sdd_in_azone, sdd_in_rzone, _, _ = zone_data
    res_req = _get_reserve_requirements(input_data, zone_data, i, p_sdd)
    (p_rru_req, p_rrd_req, q_qru_req, q_qrd_req,
     p_rgu_req, p_rgd_req, p_scr_req, p_nsc_req) = res_req
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
    p_rgu_slack = model[:p_rgu_slack]
    p_rgd_slack = model[:p_rgd_slack]
    p_scr_slack = model[:p_scr_slack]
    p_nsc_slack = model[:p_nsc_slack]
    p_rru_slack = model[:p_rru_slack]
    p_rrd_slack = model[:p_rrd_slack]
    q_qru_slack = model[:q_qru_slack]
    q_qrd_slack = model[:q_qrd_slack]
    #
    # Reserve balance constraints
    #
    # Endogenous reserve requirements
    p_rgu_balance = @constraint(model,
        p_rgu_balance[azr_id in azr_ids],
        (
            sum(p_rgu[sdd_id] for sdd_id in sdd_in_azone[azr_id], init = 0)
            + p_rgu_slack[azr_id]
        ) >= p_rgu_req[azr_id]
    )
    p_rgd_balance = @constraint(model,
        p_rgd_balance[azr_id in azr_ids],
        (
            sum(p_rgd[sdd_id] for sdd_id in sdd_in_azone[azr_id], init = 0)
            + p_rgd_slack[azr_id]
        ) >= p_rgd_req[azr_id]
    )
    p_scr_balance = @constraint(model,
        p_scr_balance[azr_id in azr_ids],
        (
            sum(
                (p_rgu[sdd_id] + p_scr[sdd_id])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_scr_slack[azr_id]
        ) >= p_rgu_req[azr_id] + p_scr_req[azr_id]
    )
    p_nsc_balance = @constraint(model,
        p_nsc_balance[azr_id in azr_ids],
        (
            sum(
                (p_rgu[sdd_id] + p_scr[sdd_id] + p_nsc[sdd_id])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_nsc_slack[azr_id]
        ) >= p_rgu_req[azr_id] + p_scr_req[azr_id] + p_nsc_req[azr_id]
    )
    # Exogenous reserve requirements
    p_rru_balance = @constraint(model,
        p_rru_balance[azr_id in azr_ids],
        (
            sum(
                (p_rru_on[sdd_id] + p_rru_off[sdd_id])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_rru_slack[azr_id]
        ) >= p_rru_req[azr_id]
    )
    p_rrd_balance = @constraint(model,
        p_rrd_balance[azr_id in azr_ids],
        (
            sum(
                (p_rrd_on[sdd_id] + p_rrd_off[sdd_id])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_rrd_slack[azr_id]
        ) >= p_rrd_req[azr_id]
    )
    q_qru_balance = @constraint(model,
        q_qru_balance[rzr_id in rzr_ids],
        (
            sum(q_qru[sdd_id] for sdd_id in sdd_in_rzone[rzr_id], init = 0)
            + q_qru_slack[rzr_id]
        ) >= q_qru_req[rzr_id]
    )
    q_qrd_balance = @constraint(model,
        q_qrd_balance[rzr_id in rzr_ids],
        (
            sum(q_qrd[sdd_id] for sdd_id in sdd_in_rzone[rzr_id], init = 0)
            + q_qrd_slack[rzr_id]
        ) >= q_qrd_req[rzr_id]
    )
    # Shouldn't need to return these balance equations
    return nothing
end

# This is the "time-indexed" implementation. No interval is provided, and p_sdd
# is an array of JuMP VariableRef (model[:p] in the scheduling model).
#
# By default, reserve requirements are relaxed and penalized in the objective.
# if relax==false, we fix the slack variables to zero to require reserve
# requirements to be strictly satisfied.
function add_reserve_balances!(
    model::Model,
    input_data::NamedTuple,
    p_sdd::JuMP.Containers.DenseAxisArray;
    relax::Bool = true,
)
    periods = input_data.periods
    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    zone_data = _get_sdd_in_zones(input_data)
    sdd_in_azone, sdd_in_rzone, c_sdd_in_azone, p_sdd_in_azone = zone_data
    # In this time-indexed implementation, we don't pass an interval to
    # _get_reserve_requirements, and p_sdd is indexed by time.
    # p_sdd shouldd not need to be used for the rest of this function.
    res_req = _get_reserve_requirements(input_data, zone_data, p_sdd)
    (p_rru_req, p_rrd_req, q_qru_req, q_qrd_req,
     p_rgu_req, p_rgd_req, p_scr_req, p_nsc_req) = res_req
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
    p_rgu_slack = model[:p_rgu_slack]
    p_rgd_slack = model[:p_rgd_slack]
    p_scr_slack = model[:p_scr_slack]
    p_nsc_slack = model[:p_nsc_slack]
    p_rru_slack = model[:p_rru_slack]
    p_rrd_slack = model[:p_rrd_slack]
    q_qru_slack = model[:q_qru_slack]
    q_qrd_slack = model[:q_qrd_slack]
    if !relax
        # Fix slack variables to zero to enforce requirements rigorously
        for uid in azr_ids
            for i in periods
                JuMP.fix(p_rgu_slack[uid, i], 0.0; force = true)
                JuMP.fix(p_rgd_slack[uid, i], 0.0; force = true)
                JuMP.fix(p_scr_slack[uid, i], 0.0; force = true)
                JuMP.fix(p_nsc_slack[uid, i], 0.0; force = true)
                JuMP.fix(p_rru_slack[uid, i], 0.0; force = true)
                JuMP.fix(p_rrd_slack[uid, i], 0.0; force = true)
            end
        end
        for uid in rzr_ids
            for i in periods
                JuMP.fix(q_qru_slack[uid, i], 0.0; force = true)
                JuMP.fix(q_qrd_slack[uid, i], 0.0; force = true)
            end
        end
    end
    #
    # Reserve balance constraints
    # Now, all reserve balance constraints need to be indexed by time.
    #
    # Endogenous reserve requirements
    p_rgu_balance = @constraint(model,
        p_rgu_balance[azr_id in azr_ids, i in periods],
        (
            sum(p_rgu[sdd_id, i] for sdd_id in sdd_in_azone[azr_id], init = 0)
            + p_rgu_slack[azr_id, i]
        ) >= p_rgu_req[azr_id][i]
    )
    p_rgd_balance = @constraint(model,
        p_rgd_balance[azr_id in azr_ids, i in periods],
        (
            sum(p_rgd[sdd_id, i] for sdd_id in sdd_in_azone[azr_id], init = 0)
            + p_rgd_slack[azr_id, i]
        ) >= p_rgd_req[azr_id][i]
    )

    #
    # These balance constraints are reformulated with an intermediate variable
    # representing the max producer active power in a given reserve zone.
    # The reformulation is done to avoid repeating the balance constraint
    # a potentially quadratic number of times. This would be a heavy performance
    # hit as the balance constraint already contains a sum over all consumers in
    # the zone (p_rgu_req), which is expensive to rewrite for a large number
    # of constraints.
    #
    max_prod_in_azone_ub = Dict(
        azr_id => [
            maximum(
                input_data.sdd_ts_lookup[sdd_id]["p_ub"][i]
                for sdd_id in p_sdd_in_azone[azr_id];
                init = 0.0
            )
            for i in periods
        ]
        for azr_id in azr_ids
    )
    @variable(model,
        0 <=
        max_prod_in_azone[azr_id in azr_ids, i in periods]
        <= max_prod_in_azone_ub[azr_id][i]
    )
    @constraint(model,
        [azr_id in azr_ids, i in periods, p_sdd_id in p_sdd_in_azone[azr_id]],
        max_prod_in_azone[azr_id, i] >= p_sdd[p_sdd_id, i]
    )
    scr_factor = Dict(uid => input_data.azr_lookup[uid]["SYN"] for uid in azr_ids)
    nsc_factor = Dict(uid => input_data.azr_lookup[uid]["NSYN"] for uid in azr_ids)
    p_scr_balance = @constraint(model,
        p_scr_balance[azr_id in azr_ids, i in periods],
        (
            sum(
                (p_rgu[sdd_id, i] + p_scr[sdd_id, i])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_scr_slack[azr_id, i]
        ) >= p_rgu_req[azr_id][i] + scr_factor[azr_id]*max_prod_in_azone[azr_id, i]
    )
    p_nsc_balance = @constraint(model,
        p_nsc_balance[azr_id in azr_ids, i in periods],
        (
            sum(
                (p_rgu[sdd_id, i] + p_scr[sdd_id, i] + p_nsc[sdd_id, i])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_nsc_slack[azr_id, i]
        ) >= (
            p_rgu_req[azr_id][i]
            + scr_factor[azr_id]*max_prod_in_azone[azr_id, i]
            + nsc_factor[azr_id]*max_prod_in_azone[azr_id, i]
        )
    )
    # Exogenous reserve requirements
    p_rru_balance = @constraint(model,
        p_rru_balance[azr_id in azr_ids, i in periods],
        (
            sum(
                (p_rru_on[sdd_id, i] + p_rru_off[sdd_id, i])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_rru_slack[azr_id, i]
        ) >= p_rru_req[azr_id][i]
    )
    p_rrd_balance = @constraint(model,
        p_rrd_balance[azr_id in azr_ids, i in periods],
        (
            sum(
                (p_rrd_on[sdd_id, i] + p_rrd_off[sdd_id, i])
                for sdd_id in sdd_in_azone[azr_id],
                init = 0
            )
            + p_rrd_slack[azr_id, i]
        ) >= p_rrd_req[azr_id][i]
    )
    q_qru_balance = @constraint(model,
        q_qru_balance[rzr_id in rzr_ids, i in periods],
        (
            sum(q_qru[sdd_id, i] for sdd_id in sdd_in_rzone[rzr_id], init = 0)
            + q_qru_slack[rzr_id, i]
        ) >= q_qru_req[rzr_id][i]
    )
    q_qrd_balance = @constraint(model,
        q_qrd_balance[rzr_id in rzr_ids, i in periods],
        (
            sum(q_qrd[sdd_id, i] for sdd_id in sdd_in_rzone[rzr_id], init = 0)
            + q_qrd_slack[rzr_id, i]
        ) >= q_qrd_req[rzr_id][i]
    )
    # Shouldn't need to return these balance equations
    return nothing
end

function make_reserve_model(
    input_data::NamedTuple,
    interval::Int64,
    opf_data::Dict,
    on_status::Dict,
    supc_status::Dict,
    sdpc_status::Dict,
    p_su::Dict,
    p_sd::Dict,
)
    #sched_sdds = opf_data["simple_dispatchable_device"]
    model = JuMP.Model()

    sdd_ids = input_data.sdd_ids
    cons_ids = input_data.sdd_ids_consumer
    prod_ids = input_data.sdd_ids_producer
    cons_set = Set(cons_ids)
    prod_set = Set(prod_ids)
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup
    # Round on_status as sometimes solvers can return inexact integer results
    on_status = Dict(sdd_id => round(on_status[sdd_id]) for sdd_id in sdd_ids)
    i = interval
    dt = input_data.dt[i]

    azr_ids = input_data.azr_ids
    rzr_ids = input_data.rzr_ids
    azr_lookup = input_data.azr_lookup
    azr_ts_lookup = input_data.azr_ts_lookup
    rzr_lookup = input_data.rzr_lookup
    rzr_ts_lookup = input_data.rzr_ts_lookup

    # These should be the original bounds, not bounds that are tightened e.g.
    # due to ramping limits.
    p_max = Dict(sdd_id => sdd_ts_lookup[sdd_id]["p_ub"][i] for sdd_id in sdd_ids)
    p_min = Dict(sdd_id => sdd_ts_lookup[sdd_id]["p_lb"][i] for sdd_id in sdd_ids)
    q_max = Dict(sdd_id => sdd_ts_lookup[sdd_id]["q_ub"][i] for sdd_id in sdd_ids)
    q_min = Dict(sdd_id => sdd_ts_lookup[sdd_id]["q_lb"][i] for sdd_id in sdd_ids)

    # Get the SDDs that are in each reserve zone
    zone_data = _get_sdd_in_zones(input_data)
    sdd_in_azone, sdd_in_rzone, c_sdd_in_azone, p_sdd_in_azone = zone_data

    sdd_key = "simple_dispatchable_device"

    #
    # Create useful "parameter" lookups
    #
    reserve_ubs = _get_max_reserves(input_data)
    (p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
     p_rru_off_max, p_rrd_off_max) = reserve_ubs
    q_res_max = Dict(
        # q reserves do not have an "absolute limit", so we use ub-lb, which
        # is also a valid bound.
        uid => (sdd_ts_lookup[uid]["q_ub"][i] - sdd_ts_lookup[uid]["q_lb"][i])
        for uid in sdd_ids
    )

    res_costs = _get_reserve_costs(input_data)
    c_rgu = Dict(uid => res_costs.c_rgu[uid][i] for uid in sdd_ids)
    c_rgd = Dict(uid => res_costs.c_rgd[uid][i] for uid in sdd_ids)
    c_scr = Dict(uid => res_costs.c_scr[uid][i] for uid in sdd_ids)
    c_nsc = Dict(uid => res_costs.c_nsc[uid][i] for uid in sdd_ids)
    c_rru_on = Dict(uid => res_costs.c_rru_on[uid][i] for uid in sdd_ids)
    c_rrd_on = Dict(uid => res_costs.c_rrd_on[uid][i] for uid in sdd_ids)
    c_rru_off = Dict(uid => res_costs.c_rru_off[uid][i] for uid in sdd_ids)
    c_rrd_off = Dict(uid => res_costs.c_rrd_off[uid][i] for uid in sdd_ids)
    c_qru = Dict(uid => res_costs.c_qru[uid][i] for uid in sdd_ids)
    c_qrd = Dict(uid => res_costs.c_qrd[uid][i] for uid in sdd_ids)

    # Zonal penalty lookups. Use z instead of c here to distinguish.
    res_penalties = _get_reserve_shortfall_penalties(input_data)
    z_rgu, z_rgd, z_scr, z_nsc, z_rru, z_rrd, z_qru, z_qrd = res_penalties

    p_on = Dict(uid => opf_data[sdd_key][uid]["p_on"] for uid in sdd_ids)
    q_sdd = Dict(uid => opf_data[sdd_key][uid]["q"] for uid in sdd_ids)
    # Note that p_su and p_sd were computed from on_status and have been
    # provided directly to this function
    p_sdd = Dict(uid => (p_on[uid] + p_su[uid] + p_sd[uid]) for uid in sdd_ids)

    #
    # Create reserve variables
    #
    res_vars = add_reserve_variables!(model, input_data, i)
    (p_rgu, p_rgd, p_scr, p_nsc, p_rru_on, p_rrd_on, p_rru_off,
     p_rrd_off, q_qru, q_qrd) = res_vars
    #
    # Create zonal reserve shortfall variables
    #
    sf_vars = add_reserve_shortfall_variables!(model, input_data, i)
    (p_rgu_slack, p_rgd_slack, p_scr_slack, p_nsc_slack, p_rru_slack,
     p_rrd_slack, q_qru_slack, q_qrd_slack) = sf_vars

    #
    # Add reserve balance (requirement) constraints
    #
    add_reserve_balances!(model, input_data, i, p_sdd)

    #
    # "Absolute" reserve limits
    #
    # Note that on_status has been rounded at this point.
    p_rgu_limit = @constraint(model,
        p_rgu_limit[sdd_id in sdd_ids],
        p_rgu[sdd_id] <= p_rgu_max[sdd_id]*on_status[sdd_id]
    )
    p_rgd_limit = @constraint(model,
        p_rgd_limit[sdd_id in sdd_ids],
        p_rgd[sdd_id] <= p_rgd_max[sdd_id]*on_status[sdd_id]
    )
    p_scr_limit = @constraint(model,
        p_scr_limit[sdd_id in sdd_ids],
        p_rgu[sdd_id] + p_scr[sdd_id] <= p_scr_max[sdd_id]*on_status[sdd_id]
    )
    p_nsc_limit = @constraint(model,
        p_nsc_limit[sdd_id in sdd_ids],
        p_nsc[sdd_id] <= p_nsc_max[sdd_id]*(1 - on_status[sdd_id])
    )
    p_rru_on_limit = @constraint(model,
        p_rru_on_limit[sdd_id in sdd_ids],
        (
            p_rgu[sdd_id] + p_scr[sdd_id] + p_rru_on[sdd_id]
        ) <= p_rru_on_max[sdd_id]*on_status[sdd_id]
    )
    p_rru_off_limit = @constraint(model,
        p_rru_off_limit[sdd_id in sdd_ids],
        (
            p_nsc[sdd_id] + p_rru_off[sdd_id]
        ) <= p_rru_off_max[sdd_id]*(1 - on_status[sdd_id])
    )
    p_rrd_on_limit = @constraint(model,
        p_rrd_on_limit[sdd_id in sdd_ids],
        (
            p_rgd[sdd_id] + p_rrd_on[sdd_id]
        ) <= p_rrd_on_max[sdd_id]*on_status[sdd_id]
    )
    p_rrd_off_limit = @constraint(model,
        p_rrd_off_limit[sdd_id in sdd_ids],
        p_rrd_off[sdd_id] <= p_rrd_off_max[sdd_id]*(1 - on_status[sdd_id])
    )
    zero_p_rrd_off_producer = @constraint(model,
        zero_p_rrd_off_producer[sdd_id in prod_ids],
        p_rrd_off[sdd_id] == 0
    )
    zero_p_nsc_consumer = @constraint(model,
        zero_p_nsc_consumer[sdd_id in cons_ids],
        p_nsc[sdd_id] == 0
    )
    zero_p_rru_off_consumer = @constraint(model,
        zero_p_rru_off_consumer[sdd_id in cons_ids],
        p_rru_off[sdd_id] == 0
    )

    #
    # "Relative" reserve limits
    #
    # Producers
    on_p_max_producer_headroom = @constraint(model,
        on_p_max_producer_headroom[sdd_id in prod_ids],
        (
            p_on[sdd_id] + p_rgu[sdd_id] + p_scr[sdd_id] + p_rru_on[sdd_id]
        ) <= p_max[sdd_id]*on_status[sdd_id]
    )
    on_p_min_producer_headroom = @constraint(model,
        p_min_on_producer_headroom[sdd_id in prod_ids],
        (
            p_on[sdd_id] - p_rgd[sdd_id] - p_rrd_on[sdd_id]
        ) >= p_min[sdd_id]*on_status[sdd_id]
    )
    off_producer_headroom = @constraint(model,
        off_producer_headroom[sdd_id in prod_ids],
        (
            p_su[sdd_id] + p_sd[sdd_id] + p_nsc[sdd_id] + p_rru_off[sdd_id]
        ) <= p_max[sdd_id]*(1 - on_status[sdd_id])
    )
    # Note that supc/sdpc_status are computed from on_status and
    # provided to this function.
    q_max_producer_headroom = @constraint(model,
        q_max_producer_headroom[sdd_id in prod_ids],
        q_sdd[sdd_id] + q_qru[sdd_id] <= q_max[sdd_id]*(
            on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
        )
    )
    q_min_producer_headroom = @constraint(model,
        q_min_producer_headroom[sdd_id in prod_ids],
        q_sdd[sdd_id] - q_qrd[sdd_id] >= q_min[sdd_id]*(
            on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
        )
    )
    # Branch on q_bound_cap/q_linear_cap on a per-generator basis to determine
    # which additional real-reactive linking constraints must be satisfied by
    # reserves
    for sdd_id in prod_ids
        q_bound_cap = sdd_lookup[sdd_id]["q_bound_cap"]
        q_linear_cap = sdd_lookup[sdd_id]["q_linear_cap"]
        if q_bound_cap == 1
            q_0_lb = sdd_lookup[sdd_id]["q_0_lb"]
            q_0_ub = sdd_lookup[sdd_id]["q_0_ub"]
            beta_lb = sdd_lookup[sdd_id]["beta_lb"]
            beta_ub = sdd_lookup[sdd_id]["beta_ub"]
            @constraint(model,
                q_sdd[sdd_id] + q_qru[sdd_id] <= q_0_ub*(
                    on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
                ) + beta_ub*p_sdd[sdd_id]
            )
            @constraint(model,
                q_sdd[sdd_id] - q_qrd[sdd_id] >= q_0_lb*(
                    on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
                ) + beta_lb*p_sdd[sdd_id]
            )
        elseif q_linear_cap == 1
            JuMP.fix(q_qru[sdd_id], 0.0; force = true)
            JuMP.fix(q_qrd[sdd_id], 0.0; force = true)
        end
    end

    # Consumers
    on_p_max_consumer_headroom = @constraint(model,
        on_p_max_consumer_headroom[sdd_id in cons_ids],
        (
            p_on[sdd_id] + p_rgd[sdd_id] + p_rrd_on[sdd_id]
        ) <= p_max[sdd_id]*on_status[sdd_id]
    )
    on_p_min_consumer_headroom = @constraint(model,
        p_min_on_consumer_headroom[sdd_id in cons_ids],
        (
            p_on[sdd_id] - p_rgu[sdd_id] - p_scr[sdd_id] - p_rru_on[sdd_id]
        ) >= p_min[sdd_id]*on_status[sdd_id]
    )
    off_consumer_headroom = @constraint(model,
        off_consumer_headroom[sdd_id in cons_ids],
        (
            p_su[sdd_id] + p_sd[sdd_id] + p_rrd_off[sdd_id]
        ) <= p_max[sdd_id]*(1 - on_status[sdd_id])
    )
    # Note that supc/sdpc_status are computed from on_status and
    # provided to this function.
    q_max_consumer_headroom = @constraint(model,
        q_max_consumer_headroom[sdd_id in cons_ids],
        q_sdd[sdd_id] + q_qrd[sdd_id] <= q_max[sdd_id]*(
            on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
        )
    )
    q_min_consumer_headroom = @constraint(model,
        q_min_consumer_headroom[sdd_id in cons_ids],
        q_sdd[sdd_id] - q_qru[sdd_id] >= q_min[sdd_id]*(
            on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
        )
    )
    # Branch on q_bound_cap/q_linear_cap on a per-generator basis to determine
    # which additional real-reactive linking constraints must be satisfied by
    # reserves
    for sdd_id in cons_ids
        q_bound_cap = sdd_lookup[sdd_id]["q_bound_cap"]
        q_linear_cap = sdd_lookup[sdd_id]["q_linear_cap"]
        if q_bound_cap == 1
            q_0_lb = sdd_lookup[sdd_id]["q_0_lb"]
            q_0_ub = sdd_lookup[sdd_id]["q_0_ub"]
            beta_lb = sdd_lookup[sdd_id]["beta_lb"]
            beta_ub = sdd_lookup[sdd_id]["beta_ub"]
            @constraint(model,
                q_sdd[sdd_id] + q_qrd[sdd_id] <= q_0_ub*(
                    on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
                ) + beta_ub*p_sdd[sdd_id]
            )
            @constraint(model,
                q_sdd[sdd_id] - q_qru[sdd_id] >= q_0_lb*(
                    on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
                ) + beta_lb*p_sdd[sdd_id]
            )
        elseif q_linear_cap == 1
            JuMP.fix(q_qru[sdd_id], 0.0; force = true)
            JuMP.fix(q_qrd[sdd_id], 0.0; force = true)
        end
    end

    #
    # Create objective
    #
    @objective(
        model,
        Max,
        # Reserve costs
        - dt*sum(
            (
                c_rgu[uid]*p_rgu[uid]
                + c_rgd[uid]*p_rgd[uid]
                + c_scr[uid]*p_scr[uid]
                + c_nsc[uid]*p_nsc[uid]
                + c_rru_on[uid]*p_rru_on[uid]
                + c_rru_off[uid]*p_rru_off[uid]
                + c_rrd_on[uid]*p_rrd_on[uid]
                + c_rrd_off[uid]*p_rrd_off[uid]
                + c_qru[uid]*q_qru[uid]
                + c_qrd[uid]*q_qrd[uid]
            ) for uid in sdd_ids;
            init = 0
        )
        # Reserve shortfall penalties
        - dt*sum(
            (
                z_rgu[uid]*p_rgu_slack[uid]
                + z_rgd[uid] * p_rgd_slack[uid]
                + z_scr[uid] * p_scr_slack[uid]
                + z_nsc[uid] * p_nsc_slack[uid]
                + z_rru[uid] * p_rru_slack[uid]
                + z_rrd[uid] * p_rrd_slack[uid]
            ) for uid in azr_ids
        )
        - dt*sum(
            (z_qru[uid]*q_qru_slack[uid] + z_qrd[uid]*q_qrd_slack[uid])
            for uid in rzr_ids
        )
    )

    return model
end

function extract_data_from_reserve_model(input_data, model::JuMP.Model)
    # input_data is only necessary to get the SDD UID keys
    sdd_ids = input_data.sdd_ids
    p_rgu = Dict(uid => value(model[:p_rgu][uid]) for uid in sdd_ids)
    p_rgd = Dict(uid => value(model[:p_rgd][uid]) for uid in sdd_ids)
    p_scr = Dict(uid => value(model[:p_scr][uid]) for uid in sdd_ids)
    p_nsc = Dict(uid => value(model[:p_nsc][uid]) for uid in sdd_ids)
    p_rru_on = Dict(uid => value(model[:p_rru_on][uid]) for uid in sdd_ids)
    p_rru_off = Dict(uid => value(model[:p_rru_off][uid]) for uid in sdd_ids)
    p_rrd_on = Dict(uid => value(model[:p_rrd_on][uid]) for uid in sdd_ids)
    p_rrd_off = Dict(uid => value(model[:p_rrd_off][uid]) for uid in sdd_ids)
    q_qru = Dict(uid => value(model[:q_qru][uid]) for uid in sdd_ids)
    q_qrd = Dict(uid => value(model[:q_qrd][uid]) for uid in sdd_ids)
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

function _curtail_reserves(
    reserves, ub;
    primal_feas_tol = 1e-6,
    acceptable_violation = 1e-10,
)
    @assert ub >= -primal_feas_tol
    slack = ub - sum(reserves; init = 0.0)
    # If constraint is violated by < our acceptable violation
    # (1e-8 in formulation), we don't need to curtail any slacks.
    # This additional tolerance was added because of the case
    # where all reserves are zero, but we have a 1e-16 violation.
    # Subtracting from the max reserve yields a negative reserve,
    # which we would like to avoid (even though we will "round up"
    # negative reserves to 0) later.
    if -slack > acceptable_violation
        # We are violating the constraint. If we are not within
        # primal_feas_tol, this is a bug.
        @assert -slack < primal_feas_tol
        n_res = length(reserves)
        # Order reserves from greatest to least
        reserve_order = reverse(sort(1:n_res, by=idx->reserves[idx]))
        # Here we assume at least one reserve was provided
        idx_of_max = reserve_order[1]
        # Slack is negative. We subtract from max reserve, assuming that
        # the result is nonnegative (reserve >= headroom violation)
        #
        # We never want this procedure to result in a negative reserve.
        @assert reserves[idx_of_max] + slack >= 0
        reserves[idx_of_max] += slack
    end
    return reserves
end

function project_reserves_onto_bounds(
    input_data::NamedTuple,
    interval::Int64,
    reserves::NamedTuple,
    opf_data::Dict,
    on_status::Dict,
    supc_status::Dict,
    sdpc_status::Dict,
    p_su::Dict,
    p_sd::Dict;
    primal_feas_tol = 1e-6,
)
    # - For each headroom constraint, compute slack
    # - If slack is negative, distribute among reserves
    # - As these constraints were imposed in the LP, any violation
    #   greater than primal_feas_tol should be considered a bug

    sdd_ids = input_data.sdd_ids
    prod_ids = input_data.sdd_ids_producer
    cons_ids = input_data.sdd_ids_consumer
    sdd_lookup = input_data.sdd_lookup
    sdd_ts_lookup = input_data.sdd_ts_lookup
    i = interval
    sdd_key = "simple_dispatchable_device"
    q_bound_ids = [uid for uid in sdd_ids if sdd_lookup[uid]["q_bound_cap"] == 1]
    q_linear_ids = [uid for uid in sdd_ids if sdd_lookup[uid]["q_linear_cap"] == 1]
    q_bound_set = Set(q_bound_ids)
    q_bound_prod_ids = [uid for uid in prod_ids if uid in q_bound_set]
    q_bound_cons_ids = [uid for uid in cons_ids if uid in q_bound_set]

    p_max = Dict(sdd_id => sdd_ts_lookup[sdd_id]["p_ub"][i] for sdd_id in sdd_ids)
    p_min = Dict(sdd_id => sdd_ts_lookup[sdd_id]["p_lb"][i] for sdd_id in sdd_ids)
    q_max = Dict(sdd_id => sdd_ts_lookup[sdd_id]["q_ub"][i] for sdd_id in sdd_ids)
    q_min = Dict(sdd_id => sdd_ts_lookup[sdd_id]["q_lb"][i] for sdd_id in sdd_ids)

    q_0_ub = Dict(uid => sdd_lookup[uid]["q_0_ub"] for uid in q_bound_ids)
    q_0_lb = Dict(uid => sdd_lookup[uid]["q_0_lb"] for uid in q_bound_ids)
    beta_lb = Dict(uid => sdd_lookup[uid]["beta_lb"] for uid in q_bound_ids)
    beta_ub = Dict(uid => sdd_lookup[uid]["beta_ub"] for uid in q_bound_ids)

    p_on = Dict(uid => opf_data[sdd_key][uid]["p_on"] for uid in sdd_ids)
    q_sdd = Dict(uid => opf_data[sdd_key][uid]["q"] for uid in sdd_ids)
    # Note that p_su and p_sd were computed from on_status and have been
    # provided directly to this function
    p_sdd = Dict(uid => (p_on[uid] + p_su[uid] + p_sd[uid]) for uid in sdd_ids)
    p_su_sd = Dict(uid => (p_su[uid] + p_sd[uid]) for uid in sdd_ids)

    pq_max = Dict(uid => q_0_ub[uid] + beta_ub[uid]*p_sdd[uid] for uid in q_bound_ids)
    pq_min = Dict(uid => q_0_lb[uid] + beta_lb[uid]*p_sdd[uid] for uid in q_bound_ids)

    p_rgu = reserves.p_rgu
    p_rgd = reserves.p_rgd
    p_scr = reserves.p_scr
    p_nsc = reserves.p_nsc
    p_rru_on = reserves.p_rru_on
    p_rru_off = reserves.p_rru_off
    p_rrd_on = reserves.p_rrd_on
    p_rrd_off = reserves.p_rrd_off
    q_qru = reserves.q_qru
    q_qrd = reserves.q_qrd

    reserve_ubs = _get_max_reserves(input_data)
    (p_rgu_max, p_rgd_max, p_scr_max, p_nsc_max, p_rru_on_max, p_rrd_on_max,
     p_rru_off_max, p_rrd_off_max) = reserve_ubs

    res_dicts = [
        p_rgu, p_rgd, p_scr, p_nsc, p_rru_on, p_rru_off, p_rrd_on, p_rrd_off,
        q_qru, q_qrd,
    ]
    # "Project up" reserve values so they are nonnegative.
    # This needs to happen *before* "projecting down" reserves, as otherwise
    # this projection can re-introduce an infeasibility.
    for res in res_dicts
        for uid in sdd_ids
            # TODO: Should this be res[uid] < some tolerance?
            if res[uid] < 0.0
                res[uid] = 0.0
            end
        end
    end

    # (sdd type, reserves, LHS, bound, bound type, required status)
    reserves_ubs = [
        # Absolute reserve limits
        ("any", [p_rgu],                  nothing, p_rgu_max,     "ub", "on"),
        ("any", [p_rgd],                  nothing, p_rgd_max,     "ub", "on"),
        ("any", [p_rgu, p_scr],           nothing, p_scr_max,     "ub", "on"),
        ("any", [p_rgu, p_scr, p_rru_on], nothing, p_rru_on_max,  "ub", "on"),
        ("any", [p_nsc],                  nothing, p_nsc_max,     "ub", "off"),
        ("any", [p_nsc, p_rru_off],       nothing, p_rru_off_max, "ub", "off"),
        ("any", [p_rgd, p_rrd_on],        nothing, p_rrd_on_max,  "ub", "on"),
        ("any", [p_rrd_off],              nothing, p_rrd_off_max, "ub", "off"),
        # This enforces an upper bound of zero
        ("producer", [p_rrd_off],         nothing, nothing,       "ub", "any"),
        ("consumer", [p_nsc],             nothing, nothing,       "ub", "any"),
        ("consumer", [p_rru_off],         nothing, nothing,       "ub", "any"),

        # Producers: active reserves
        ("producer", [p_rgu, p_scr, p_rru_on], p_on, p_max, "ub", "on"),
        ("producer", [p_rgd, p_rrd_on], p_on, p_min, "lb", "on"),
        ("producer", [p_nsc, p_rru_off], p_su_sd, p_max, "ub", "off"),

        # Producers: reactive reserves
        ("producer", [q_qru], q_sdd, q_max, "ub", "on_su_sd"),
        ("producer", [q_qrd], q_sdd, q_min, "lb", "on_su_sd"),
        ("qbd-producer", [q_qru], q_sdd, pq_max, "ub", "on_su_sd"),
        ("qbd-producer", [q_qrd], q_sdd, pq_min, "lb", "on_su_sd"),

        # Consumers: active reserves
        ("consumer", [p_rgd, p_rrd_on], p_on, p_max, "ub", "on"),
        ("consumer", [p_rgu, p_scr, p_rru_on], p_on, p_min, "lb", "on"),
        ("consumer", [p_rrd_off], p_su_sd, p_max, "ub", "off"),

        # Consumer: reactive reserves
        ("consumer", [q_qrd], q_sdd, q_max, "ub", "on_su_sd"),
        ("consumer", [q_qru], q_sdd, q_min, "lb", "on_su_sd"),
        ("qbd-consumer", [q_qrd], q_sdd, pq_max, "ub", "on_su_sd"),
        ("qbd-consumer", [q_qru], q_sdd, pq_min, "lb", "on_su_sd"),

        # Do not need to distinguish between producer and consumer here
        ("qlin", [q_qru], nothing, nothing, "ub", "any"),
        ("qlin", [q_qrd], nothing, nothing, "ub", "any"),
    ]

    for (idx, res_data) in enumerate(reserves_ubs)
        (sdd_type, reserve_dicts, baseline, bound, bd_type, status) = res_data
        if sdd_type == "consumer"
            ids = cons_ids
        elseif sdd_type == "producer"
            ids = prod_ids
        elseif sdd_type == "qbd-producer"
            ids = q_bound_prod_ids
        elseif sdd_type == "qbd-consumer"
            ids = q_bound_cons_ids
        elseif sdd_type == "any"
            ids = sdd_ids
        elseif sdd_type == "qlin"
            ids = q_linear_ids
        end
        for sdd_id in ids
            if bd_type == "ub"
                res_lhs = (baseline === nothing ? 0.0 : baseline[sdd_id])
                res_bound = (bound === nothing ? 0.0 : bound[sdd_id])
                res_ub = res_bound - res_lhs
            elseif bd_type == "lb"
                res_lhs = (baseline === nothing ? 0.0 : baseline[sdd_id])
                res_bound = (bound === nothing ? 0.0 : bound[sdd_id])
                res_ub = res_lhs - res_bound
            end

            # Status here is the "required status" in order to provide these
            # reserves.
            if status == "on"
                applicable = Bool(round(on_status[sdd_id]))
            elseif status == "off"
                applicable = !Bool(round(on_status[sdd_id]))
            elseif status == "on_su_sd"
                applicable = Bool(round(
                    on_status[sdd_id] + supc_status[sdd_id] + sdpc_status[sdd_id]
                ))
            elseif status == "any"
                applicable = true
            end

            res_values = [res[sdd_id] for res in reserve_dicts]
            if applicable
                proj_reserves = _curtail_reserves(
                    res_values, res_ub, primal_feas_tol = primal_feas_tol
                )
                for (res, val) in zip(reserve_dicts, proj_reserves)
                    res[sdd_id] = val
                end
            else
                # Force reserves to zero if this SDD cannot provide these
                # reserves at this interval (e.g. it is offline)
                for res in reserve_dicts
                    res[sdd_id] = 0.0
                end
            end
        end
    end

    return reserves
end

function calculate_reserves_from_generation(
    input_data::NamedTuple,
    interval::Int64,
    # These arguments could be grouped together into a "realized power
    # curve data" NamedTuple, as they are used together frequently.
    opf_data::Dict,
    on_status::Dict,
    supc_status::Dict,
    sdpc_status::Dict,
    p_su::Dict,
    p_sd::Dict;
    optimizer = HiGHS.Optimizer,
    set_silent = false,
    return_model = true,
)
    model = make_reserve_model(
        input_data,
        interval,
        opf_data,
        on_status,
        supc_status,
        sdpc_status,
        p_su,
        p_sd,
    )
    JuMP.set_optimizer(model, optimizer)
    if set_silent
        JuMP.set_silent(model)
    end
    optimize!(model)
    data = extract_data_from_reserve_model(input_data, model)
    data = project_reserves_onto_bounds(
        input_data,
        interval,
        data,
        opf_data,
        on_status,
        supc_status,
        sdpc_status,
        p_su,
        p_sd,
    )
    return (model, data)
end

"""
Gets p_sdd in a convenient format for calling _get_reserve_requirements.
The input is in a similar format to that required by construct_solution_dict.
"""
function _get_p_sdd(
    input_data::NamedTuple,
    on_status::Dict,
    opf_data::Vector,
)::Vector{Dict{String, Float64}}
    supc_status, p_su, sdpc_status, p_sd = get_supc_sdpc_lookups(
        input_data, on_status
    )
    uids = input_data.sdd_ids
    periods = input_data.periods
    p_su = [Dict(uid => p_su[uid][i] for uid in uids) for i in periods]
    p_sd = [Dict(uid => p_sd[uid][i] for uid in uids) for i in periods]
    sdd_key = "simple_dispatchable_device"
    p_on = [
        Dict(uid => opf_data[i][sdd_key][uid]["p_on"] for uid in uids)
        for i in periods
    ]
    p_sdd = [
        Dict(uid => p_su[i][uid] + p_sd[i][uid] + p_on[i][uid] for uid in uids)
        for i in periods
    ]
end

function _solve_reserve_models(
    input_data::NamedTuple,
    opf_data::Vector,
    on_status::Dict;
    optimizer = HiGHS.Optimizer,
    return_models = true,
)
    uids = input_data.sdd_ids
    periods = input_data.periods
    # Note that this function assumes on_status is Dict{String, Vector}
    supc_status, p_su, sdpc_status, p_sd = get_supc_sdpc_lookups(
        input_data, on_status
    )
    # "Reshape" dicts of vectors in to vectors of dicts so data is easier
    # for the reserve model to deal with.
    on_status = [Dict(uid => on_status[uid][i] for uid in uids) for i in periods]
    supc_status = [Dict(uid => supc_status[uid][i] for uid in uids) for i in periods]
    sdpc_status = [Dict(uid => sdpc_status[uid][i] for uid in uids) for i in periods]
    p_su = [Dict(uid => p_su[uid][i] for uid in uids) for i in periods]
    p_sd = [Dict(uid => p_sd[uid][i] for uid in uids) for i in periods]

    model_data_array = Vector{Any}([nothing for _ in periods])
    println("BEGINNING PARALLEL RESERVE CALCULATIONS")
    Threads.@threads for i in periods
        model_data_array[i] = calculate_reserves_from_generation(
            input_data,
            i,
            opf_data[i],
            on_status[i],
            supc_status[i],
            sdpc_status[i],
            p_su[i],
            p_sd[i];
            optimizer = optimizer,
            # If we are computing reserves in parallel, don't display solver log
            set_silent = true,
            return_model = return_models,
        )
    end
    println("DONE WITH PARALLEL RESERVE CALCULATIONS")
    return model_data_array
end

function calculate_reserves_from_generation(
    input_data::NamedTuple,
    opf_data::Vector,
    on_status::Dict;
    optimizer = HiGHS.Optimizer,
)
    model_data_array = _solve_reserve_models(
        input_data,
        opf_data,
        on_status;
        optimizer = optimizer,
        return_models = false,
    )
    data = [data for (_, data) in model_data_array]

    uids = input_data.sdd_ids
    periods = input_data.periods
    p_rgu = Dict(uid => [data[i].p_rgu[uid] for i in periods] for uid in uids)
    p_rgd = Dict(uid => [data[i].p_rgd[uid] for i in periods] for uid in uids)
    p_scr = Dict(uid => [data[i].p_scr[uid] for i in periods] for uid in uids)
    p_nsc = Dict(uid => [data[i].p_nsc[uid] for i in periods] for uid in uids)
    p_rru_on = Dict(uid => [data[i].p_rru_on[uid] for i in periods] for uid in uids)
    p_rru_off = Dict(uid => [data[i].p_rru_off[uid] for i in periods] for uid in uids)
    p_rrd_on = Dict(uid => [data[i].p_rrd_on[uid] for i in periods] for uid in uids)
    p_rrd_off = Dict(uid => [data[i].p_rrd_off[uid] for i in periods] for uid in uids)
    q_qru = Dict(uid => [data[i].q_qru[uid] for i in periods] for uid in uids)
    q_qrd = Dict(uid => [data[i].q_qrd[uid] for i in periods] for uid in uids)
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

function _get_mock_opf_data(
    input_data::NamedTuple,
    solution_dict::Dict{String, Any},
)::Vector{Dict{String, Any}}
    tsout = solution_dict["time_series_output"]
    sdd_key = "simple_dispatchable_device"
    sdd_array = tsout["simple_dispatchable_device"]
    pq_uid_lookups = [
        Dict(
            sdd["uid"] => Dict("p_on" => sdd["p_on"][i], "q" => sdd["q"][i])
            for sdd in sdd_array
        )
        for i in input_data.periods
    ]
    # Mock up the "opf data" array returned by the OPF functions with the minimum
    # data required for reserve calculations
    opf_data = [Dict(sdd_key => pq_uid_lookups[i]) for i in input_data.periods]
    return opf_data
end

function calculate_reserves_from_generation(
    input_data::NamedTuple,
    solution_dict::Dict{String, Any};
    optimizer = HiGHS.Optimizer,
)::NamedTuple
    tsout = solution_dict["time_series_output"]
    sdd_key = "simple_dispatchable_device"
    sdd_array = tsout["simple_dispatchable_device"]
    on_status = Dict(sdd["uid"] => sdd["on_status"] for sdd in sdd_array)
    opf_data = _get_mock_opf_data(input_data, solution_dict)
    return calculate_reserves_from_generation(
        input_data,
        opf_data,
        on_status;
        optimizer = optimizer,
    )
end
