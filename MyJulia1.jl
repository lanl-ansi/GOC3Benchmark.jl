start_init = time()
start_pkg = time()

using Pkg
Pkg.activate(DEPOT_PATH[1])

include("scripts/common.jl")
#include("ac-solver-011.jl")
include("scripts/ac-uc-solver.jl")

println("script startup time: $(time() - start_init)")
println("")

function MyJulia1(case::String, time_limit::Int, division::Int, model::String, switching_allowed::Int)
    println("running MyJulia1")
    println("  $(case)")
    println("  $(time_limit)")
    println("  $(division)")
    println("  $(model)")
    println("  $(switching_allowed)")

    startup_time = time() - start_init
    remaining_time = trunc(Int, time_limit-startup_time-10)
    println("")
    println("time limit: $(time_limit)")
    println("startup time: $(startup_time)")
    println("remaining time: $(remaining_time)")

    code1(case, remaining_time, division, model, switching_allowed, "solution.json")
end
