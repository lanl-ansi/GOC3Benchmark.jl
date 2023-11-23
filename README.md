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
```julia
julia> ]
(@v1.X) pkg> add .
```
To make sure that the package is installed correctly, please run the tests with
```
$ julia test/runtests.jl
```

## Using the solver
Challenge 3 of the GO Competition expects Julia submissions in the form of a
single `MyJulia1.jl` file, which contains a function `MyJulia1`. This file
is provided in the top-level of this repository. The APIs expected by the
competition evaluation platform are documented
[here](https://gocompetition.energy.gov/languages).
The solver may be run with:
```
$ julia --compiled-modules=no -e 'include("MyJulia1.jl"); MyJulia1(ProblemFile, TimeLimitInSeconds, Division, NetworkModel, AllowSwitching)' &>MyJulia1.log
```
where arguments to `MyJulia1` are as follows:
| Argument | Julia type | Example |
| -------- | ---------- | ------- |
| `ProblemFile` | `String` | `"test/data/C3E4N00073D1_scenario_303.json"` |
| `TimeLimitInSeconds` | `Int` | 600 |
| `Division` | `Int` | 1 |
| `NetworkModel` | `String` | `"C3E4N00073"` |
| `AllowSwitching` | `Int` | 1 |

When running via this API, solutions are written to `solution.json` in the
working directory, as in the competition.
Note that, in this solver, the "basename" of the `ProblemFile` string **must contain**
the network name as a 6-character substring `N#####`. In this solver,
the number of buses is extracted from this string and used to set some parameters.

For a more realistic representation of the benchmark code used by the
competition, set the `JULIA_NUM_THREADS` environment variable before
running the solver. For the competition,
`JULIA_NUM_THREADS=50` was used to allow all 48 ACOPF subproblems to run in
parallel for Division 2, although the optimal number of threads will of course
depend on your machine.

Alternatively to using `MyJulia1.jl`, the solver may be run via the
`scripts/ac-uc-solver.jl` script.
The only required argument to this script is `-c CASE`, which must specify
the path to the input data file.
For example:
```
$ julia scripts/ac-uc-solver.jl -c test/data/C3E4N00073D1_scenario_303.json
```
This script writes the solution file to the input file's directory, with the
same name as the input file other than `.json` replaced with `_solution.json`.
The script will remove solution files if the `--remove-solution` argument is set.

## Structure of this repository
