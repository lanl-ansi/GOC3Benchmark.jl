# GOC Challenge 3 Benchmark Algorithm
The ARPA-E Benchmark Algorithm for Challenge 3 of the Grid Optimization
Competition. This repository contains Julia code and a command line-executable
Julia script for solving the GOC-3 AC Unit Commitment problem.

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
- Gurobi.jl (only necessary for using the command line-executable script)
- Ipopt.jl
- MathOptSymbolicAD.jl
- ArgParse.jl
- JSON.jl
- Printf.jl

This package is not registered and must be installed from this repository.
For example:
```julia
julia> ]
(@v1.X) pkg> add https://github.com/lanl-ansi/GOC3Benchmark.jl.git
```
To make sure that the package is installed correctly, please run the tests with
```julia
julia> ]
(@v1.X) pkg> test GOC3Benchmark
```

## Using the solver
Challenge 3 of the GO Competition expects Julia submissions in the form of a
single `MyJulia1.jl` file, which contains a function `MyJulia1`. This file
is provided in the top-level of this repository. The APIs expected by the
competition evaluation platform are documented
[here](https://gocompetition.energy.gov/languages).
The solver may be run with:
```shell
julia --compiled-modules=no -e 'include("MyJulia1.jl"); MyJulia1(ProblemFile, TimeLimitInSeconds, Division, NetworkModel, AllowSwitching)' &>MyJulia1.log
```
where arguments to `MyJulia1` are as follows:
| Argument | Julia type | Example |
| -------- | ---------- | ------- |
| `ProblemFile` | `String` | `"test/data/C3E4N00073D1_scenario_303.json"` |
| `TimeLimitInSeconds` | `Int` | 600 |
| `Division` | `Int` | 1 |
| `NetworkModel` | `String` | `"C3E4N00073"` |
| `AllowSwitching` | `Int` | 1 |

For example:
```shell
julia --compiled-modules=no -e 'include("MyJulia1.jl"); MyJulia1("test/data/C3E4N00073D1_scenario_303.json", 600, 1, "C3E4N00073", 1)' &>MyJulia1.log
```

When running via this API, solutions are written to `solution.json` in the
working directory, as in the competition.
Note that, in this solver, the "basename" of the `ProblemFile` string **must contain**
the network name as a 6-character substring `N#####`. In this solver,
the number of buses is extracted from this string and used to set some parameters.

For a more realistic representation of the benchmark code used by the
competition, set the `JULIA_NUM_THREADS` environment variable before
running the solver. For the competition,
`JULIA_NUM_THREADS=50` was used to allow all 48 ACOPF subproblems to run in
parallel for Division 2.

Alternatively to using `MyJulia1.jl`, the solver may be run via the
`scripts/ac-uc-solver.jl` script.
The only required argument to this script is `-c CASE`, which must specify
the path to the input data file.
For example:
```
$ julia scripts/ac-uc-solver.jl -c test/data/C3E4N00073D1_scenario_303.json
```
For a list of all options, run
```
$ julia scripts/ac-uc-solver.jl --help
```

This script writes the solution file to the input file's directory, with the
same name as the input file other than `.json` replaced with `_solution.json`.
The script will remove solution files if the `--remove-solution` argument is set.

If the `--evaluate-solution` argument is set, this script will
use the `C3DataUtilities` Python package to evaluate the solution and
display results. This relies on `PyCall.jl`, and on the ability to find
a Python installation with access to the `C3DataUtilities` package.

The default solver for the unit commitment subproblem (when running via the
command line) is Gurobi. If Gurobi is not available, and you would like to run
with an open-source MIP solver, use the `--mip-solver=highs` option.
Note that the tests do not use Gurobi.

## Structure of this repository

Of the Julia "library" code contained in this package, the `run_ac_uc_solver`
function in the `src/ac_uc_solver.jl` file is the primary driver.
This driver operates primarily by calling "subroutines,"  which accept
input data and return output data for the three subproblems:
copper-plate unit commitment, ACOPF, and reserve allocation.
The "subroutines" may be found in the `scheduler.jl`, `opf.jl`, and `reserves.jl`
files. Each of these subproblems build and solve one or more JuMP models.
The code to construct the models themselves may be found in the
`scheduling_model.jl`, `opf_model.jl`, and `reserves.jl` files.

## Citing this repository

If you use this software in your research, we would appreciate you citing the following
publication:
```bibtex
@inproceedings{parker2024goc,
author = {Parker, Robert and Coffrin, Carleton},
title = {Managing Power Balance and Reserve Feasibility in the {AC} Unit Commitment Problem},
booktitle = {2024 Power Systems Computation Conference (PSCC)},
year = {2024},
month = {June},
doi = {https://doi.org/10.1016/j.epsr.2024.110670}
}
```

## License
This software is provided under a BSD license as part of the Grid Optimization
Competition Solvers project, C19076. See LICENSE.md.
