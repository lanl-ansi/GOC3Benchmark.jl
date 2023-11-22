using Test

Test.@testset "GOC-3-Benchmark" begin

    include("ac_uc_solve_io.jl")

    include("solve_models.jl")

end
