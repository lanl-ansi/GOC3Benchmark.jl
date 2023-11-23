# GOC Challenge 3 Benchmark Algorithm
The ARPA-E Benchmark Algorithm for Challenge 3 of the Grid Optimization
Competition. This repository contains Julia code and command line-executable
Julia scripts for solving the GOC-3 AC Unit Commitment problem.

The GOC-3 AC Unit Commitment problem formulation may be found on the GO Competition
[website](https://gocompetition.energy.gov/challenges/challenge-3/formulation).
A solver for the competition problem accepts input data as a JSON file and produces
its output as another JSON file.
Solutions may be evaluated using the
[C3DataUtilities](https://github.com/GOCompetition/C3DataUtilities)
Python package. The JSON files accepted and produced by the solver are described
[here](https://gocompetition.energy.gov/challenges/challenge-3/data_format).
Problem data files used for testing are included in the `test/data` directory of
this repository. All problem files used in the competition are available
[here](https://gocompetition.energy.gov/challenges/600650/datasets).

## Installing this package
This package has the following dependencies:
- JuMP.jl
- HiGHS.jl
- Ipopt.jl
- MathOptSymbolicAD.jl
- ArgParse.jl
- JSON.jl
- Printf.jl

This package is not registered and must be installed using a copy of this repository.
For example:
```
$ git clone https://github.com/lanl-ansi/GOC3Benchmark.jl.git
$ cd GOC3Benchmark.jl
```
```julia-repl
julia> ]
(@v1.X) pkg> add .
```
To make sure that the package is installed correctly, please run the tests with
```
$ julia test/runtests.jl
```

## Using the solver
The primary method of calling this solver is to use one the the
command line-executable scripts in the `scripts` directory.
Alternatively, to reproduce the calling method used by the GOC evaluation
platform, you may interact with the `scripts/MyJulia1.jl` script directly.
```
$ julia --compiled-modules=no -e 'include("MyJulia1.jl"); MyJulia1(ProblemFile, TimeLimitInSeconds, Division, NetworkModel, AllowSwitching)' &>MyJulia1.log
```
where arguments to `MyJulia1` are as follows:
- `ProblemFile`
- `TimeLimitInSeconds`
- `Division`
- `NetworkModel`
- `AllowSwitching`

When running via this API, solutions are written to `solution.json` in the
working directory, as in the competition.

For a more realistic representation of the benchmark code used by the
competition, set the `JULIA_NUM_THREADS` environment variable before
running the solver via either of the above methods. For the competition,
`JULIA_NUM_THREADS=50` was used to allow all 48 ACOPF subproblems to run in
parallel for Division 2, although the optimal number of threads will of course
depend on your machine.

## Test datasets

## Structure of this repository
