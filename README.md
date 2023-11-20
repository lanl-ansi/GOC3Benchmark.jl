# GOC Challenge 3 Benchmark Algorithm
The ARPA-E Benchmark Algorithm for Challenge 3 of the Grid Optimization
Competition. This repository contains Julia code and command line-executable
Julia scripts for solving the GOC-3 AC Unit Commitment problem.

## Installing this package
- Dependencies
- Installation
- Run tests

## Using the solver
The primary method of calling this solver is to use one the the
command line-executable scripts in the `scripts` directory.
Alternatively, to reproduce the calling method used by the GOC evaluation
platform, you may interact with the `scripts/MyJulia1.jl` script directly.
```
$ julia --compiled-modules=no -e `include("MyJulia1.jl"); MyJulia1(ProblemFile, TimeLimitInSeconds, Division, NetworkModel, AllowSwitching)' &>MyJulia1.log
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
