#using Distributed
using ArgParse
using JSON

"""common.jl
This file contains common option parsing functionality to help drive command
line-executable Julia scripts that execute the benchmark algorithm.
These scripts may be used as potentially more user-friendly alternatives
to running the MyJulia1.jl script directly.
"""


function main(args)
    case = args["case"]

    division = 2
    model = "MODEL"
    switching_allowed = 1

    output_file_path = replace(case, ".json" => "_solution.json")

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

    time_start = time()
    code1(case, args["time-limit"], division, model, switching_allowed, output_file_path, use_hsl=!args["no-hsl"])
    code1_time = time() - time_start

    #evaluation_summary(args["case"], length(Distributed.workers()), code1_time)

    if args["remove-solution"]
        @warn("removing solution file $(output_file_path)")
        rm(output_file_path)
    end
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
    end

    return parse_args(s)
end
