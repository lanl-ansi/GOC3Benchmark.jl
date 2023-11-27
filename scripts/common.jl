#using Distributed
using ArgParse
using JSON

using HiGHS
using Gurobi

"""common.jl
This file contains common option parsing functionality to help drive command
line-executable Julia scripts that execute the benchmark algorithm.
These scripts may be used as potentially more user-friendly alternatives
to running the MyJulia1.jl script directly.
"""

function get_output_file_path(input_file_path::String)
    return replace(input_file_path, ".json" => "_solution.json")
end

function get_optimizer(solver_name::String)
    name_to_optimizer = Dict(
        "gurobi" => Gurobi.Optimizer,
        "highs" => HiGHS.Optimizer,
    )
    if !(solver_name in keys(name_to_optimizer))
        throw(ArgumentError(
            "Unsupported solver specified. Supported options are: gurobi, highs"
        ))
    end
    return name_to_optimizer[solver_name]
end

function main(args)
    case = args["case"]

    division = 2
    model = "MODEL"
    switching_allowed = 1

    output_file_path = get_output_file_path(case)

    println("case size: $(filesize(case))")

    println("output file: $(output_file_path)")
    println("no hsl: $(args["no-hsl"])")

    println()
    println("running code1 with:")
    println("  $(case)")
    println("  $(args["time-limit"])")
    println("  $(division)")
    println("  $(model)")
    println("  $(switching_allowed)")

    println("MIP solver: $(args["mip-solver"])")

    mip_optimizer = get_optimizer(args["mip-solver"])

    time_start = time()
    code1(
        case,
        args["time-limit"],
        division,
        model,
        switching_allowed,
        output_file_path,
        use_hsl=!args["no-hsl"],
        mip_optimizer = mip_optimizer,
    )
    code1_time = time() - time_start
end


function parse_goc_c3_args()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--case", "-c"
            help = "the case file"
            required = true
        "--time-limit", "-t"
            help = "solution 1 runtime limit (seconds)"
            arg_type = Int
            default = 600000 #2700
        "--hard-time-limit"
            help = "solution 1 and 2 runtime limit (seconds)"
            arg_type = Int
            default = 600000 #34200
        "--no-hsl"
            help = "removes the use of hsl with Ipopt"
            action = :store_true
        "--remove-solution"
            help = "delete the solution files after competition"
            action = :store_true
        "--evaluate-solution"
            help = "Evaluate the solution using the C3DataUtilities program.
                    For this to work, C3DataUtilities must be installed in
                    the Python installation that PyCall knows about."
            action = :store_true
        "--mip-solver"
            help = "Solver to use for the MILP unit commitment subproblem.
                    (default: gurobi)"
            arg_type = String
            default = "gurobi"
    end

    return parse_args(s)
end
