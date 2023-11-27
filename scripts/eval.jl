# First time setup
#
# NOTE, not able to use Conda solution becouse data utils does not exist in Conda package manager
#
# ] add PyCall
# ENV["PYTHON"] = "python3"
# Pkg.build("PyCall")
#

using PyCall

# With this line, we would rely on C3DataUtilities being located in a specific
# location. We will instead rely on it having been installed
#pushfirst!(PyVector(pyimport("sys")."path"), joinpath(@__DIR__, "..", "C3DataUtilities"))

goc_checkdata = pyimport("check_data")
goc_validation = pyimport("datautilities.validation")

check_data_configuration = PyObject(nothing)
check_data_pop_solution = PyObject(nothing)
check_data_parameters = "{}"

# Again, we avoid relying on C3DataUtilities's location
#config_file_path = joinpath(@__DIR__, "..", "C3DataUtilities", goc_checkdata.default_config_file)
config_file_path = joinpath(dirname(goc_checkdata.__file__), goc_checkdata.default_config_file)

function evaluation_summary(case_file::AbstractString, solution_file::AbstractString, workers::Int, solver_runtime::Real)

    summary_csv_file_path = replace(case_file, ".json" => "_$(goc_checkdata.summary_csv_file)")
    summary_json_file_path = replace(case_file, ".json" => "_$(goc_checkdata.summary_json_file)")
    data_errors_file_path = replace(case_file, ".json" => "_$(goc_checkdata.data_errors_file)")
    ignored_errors_file_path = replace(case_file, ".json" => "_$(goc_checkdata.ignored_errors_file)")
    solution_errors_file_path = replace(case_file, ".json" => "_$(goc_checkdata.solution_errors_file)")

    result = goc_validation.check_data(
        case_file,
        solution_file,
        config_file_path,
        check_data_configuration,
        check_data_parameters,
        summary_csv_file_path,
        summary_json_file_path,
        data_errors_file_path,
        ignored_errors_file_path,
        solution_errors_file_path,
        check_data_pop_solution,
    )

    case_id = replace(basename(case_file), ".json" => "")

    evaluation = result["evaluation"]

    feasible = evaluation["feas"]
    objective = evaluation["z"]
    #objective_base = evaluation["z_base"]
    objective_value = evaluation["z_value"]
    objective_cost = -evaluation["z_cost"]
    objective_z_penalty = -evaluation["z_penalty"]
    objective_p_penalty = -evaluation["sum_bus_t_z_p"]
    objective_q_penalty = -evaluation["sum_bus_t_z_q"]
    objective_reserve_penalty = -(
        evaluation["sum_prz_t_z_nsc"] +
        evaluation["sum_prz_t_z_rgd"] +
        evaluation["sum_prz_t_z_rgu"] +
        evaluation["sum_prz_t_z_rrd"] +
        evaluation["sum_prz_t_z_rru"] +
        evaluation["sum_prz_t_z_scr"] +
        evaluation["sum_qrz_t_z_qrd"] +
        evaluation["sum_qrz_t_z_qru"])
    objective_energy_penalty = -(evaluation["z_max_energy"] + evaluation["z_min_energy"])
    objective_cont_avg = evaluation["z_k_average_case"]
    objective_cont_max = evaluation["z_k_worst_case"]

    data = [
        "--------",
        "case id",
        #"num bus",
        #"num branch",
        #"num sdd",
        #"num contingencies",
        "feasible",
        "objective",
        "objective_value",
        "objective_cost",
        "objective_penalty",
        "p_balance_penalty",
        "q_balance_penalty",
        "reserve_penalty",
        "energy_penalty",
        "cont_penalty_avg",
        "cont_penalty_max",
        "workers",
        "solver runtime (sec.)",
    ]
    println(join(data, ", "))

    data = [
        "DATA_EVAL",
        case_id,
        #length(pm_data["bus"]),
        #length(pm_data["branch"]),
        #length(pm_data["sdd"]),
        #length(pm_data["contingencies"]),

        feasible,
        round(Int,objective),
        round(Int,objective_value),
        round(Int,objective_cost),
        round(Int,objective_z_penalty),
        round(Int,objective_p_penalty),
        round(Int,objective_q_penalty),
        round(Int,objective_reserve_penalty),
        round(Int,objective_energy_penalty),
        round(Int,objective_cont_avg),
        round(Int,objective_cont_max),

        workers,
        round(Int,solver_runtime),
    ]
    println(join(data, ", "))

    return result
end
