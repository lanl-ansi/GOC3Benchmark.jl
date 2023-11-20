using JSON

"""
Process data from input (problem) file into a more convenient form.

`data` should be the dict obtained by parsing an input JSON file.
The return type is a NamedTuple with fields that are useful for constructing
the scheduling and ACOPF models.

"""
function process_input_data(data::Dict)::NamedTuple
    dt = data["time_series_input"]["general"]["interval_duration"]
    periods = 1:data["time_series_input"]["general"]["time_periods"]
    @assert(length(dt) == length(periods))

    bus_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["bus"])
    # NOTE: We never check inclusion in these IDs, so no reason for them
    # to be arbitrary-order KeySets as opposed to sorted Vectors.
    bus_ids = sort([uid for uid in keys(bus_lookup)])

    shunt_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["shunt"])
    shunt_ids = sort([uid for uid in keys(shunt_lookup)])

    ac_line_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["ac_line"])
    ac_line_ids = sort([uid for uid in keys(ac_line_lookup)])

    twt_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["two_winding_transformer"])
    twt_ids = sort([uid for uid in keys(twt_lookup)])

    dc_line_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["dc_line"])
    dc_line_ids = sort([uid for uid in keys(dc_line_lookup)])

    sdd_lookup = Dict(sdd["uid"] => sdd for sdd in data["network"]["simple_dispatchable_device"])
    sdd_ts_lookup = Dict(sdd["uid"] => sdd for sdd in data["time_series_input"]["simple_dispatchable_device"])
    sdd_ids = sort([uid for uid in keys(sdd_lookup)])

    sdd_ids_producer = sort([uid for (uid,sdd) in sdd_lookup if sdd["device_type"]=="producer"])
    sdd_ids_consumer = sort([uid for (uid,sdd) in sdd_lookup if sdd["device_type"]=="consumer"])

    azr_key = "active_zonal_reserve"
    azr_lookup = Dict(zone["uid"] => zone for zone in data["network"][azr_key])
    azr_ts_lookup = Dict(zone["uid"] => zone for zone in data["time_series_input"][azr_key])
    azr_ids = sort([zone for zone in keys(azr_lookup)])

    rzr_key = "reactive_zonal_reserve"
    rzr_lookup = Dict(zone["uid"] => zone for zone in data["network"][rzr_key])
    rzr_ts_lookup = Dict(zone["uid"] => zone for zone in data["time_series_input"][rzr_key])
    rzr_ids = sort([zone for zone in keys(rzr_lookup)])

    data_tuple = (
        dt = dt,
        periods = periods,
        bus_lookup = bus_lookup,
        bus_ids = bus_ids,
        shunt_lookup = shunt_lookup,
        shunt_ids = shunt_ids,
        ac_line_lookup = ac_line_lookup,
        ac_line_ids = ac_line_ids,
        twt_lookup = twt_lookup,
        twt_ids = twt_ids,
        dc_line_lookup = dc_line_lookup,
        dc_line_ids = dc_line_ids,
        sdd_lookup = sdd_lookup,
        sdd_ts_lookup = sdd_ts_lookup,
        sdd_ids = sdd_ids,
        sdd_ids_producer = sdd_ids_producer,
        sdd_ids_consumer = sdd_ids_consumer,
        violation_cost = data["network"]["violation_cost"],
        azr_lookup = azr_lookup,
        azr_ts_lookup = azr_ts_lookup,
        azr_ids = azr_ids,
        rzr_lookup = rzr_lookup,
        rzr_ts_lookup = rzr_ts_lookup,
        rzr_ids = rzr_ids,
    )
    return data_tuple
end


function get_data_from_file(input_file::String; process_data = true)
    data = Dict{String, Any}()
    open(input_file, "r") do io
        data = JSON.parse(io)
    end
    if process_data
        # NOTE: Type-unstable
        data = process_input_data(data)
    end
    return data
end
