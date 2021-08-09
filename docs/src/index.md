# HarmonicPowerFlow.jl
[![codecov](https://codecov.io/gh/pweigmann/HarmonicPowerFlow.jl/branch/master/graph/badge.svg?token=7DZTYSH7TY)](https://codecov.io/gh/pweigmann/HarmonicPowerFlow.jl)
[![CI](https://github.com/pweigmann/HarmonicPowerFlow.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/pweigmann/HarmonicPowerFlow.jl/actions/workflows/CI.yml)
[![Code on Github.](https://img.shields.io/badge/code%20on-github-blue.svg)](https://github.com/pweigmann/HarmonicPowerFlow.jl)

An open source module to solve the harmonic power flow problem in Julia. Nonlinear devices are represented by their Norton equivalents using either a coupled
or uncoupled approach. A system of nonlinear equations is assembled and solved using the Newton-Raphson algorithm.

This package was created as part of a master's thesis at TU Berlin, in collaboration with [DAI-Labor](https://dai-labor.de/) and [elena international](https://www.elena-international.com/).

Comments, questions and contributions are highly welcome.

# Usage

Install the module via the package manager:

```
]add HarmonicPowerFlow
```

or

```Julia
using Pkg; Pkg.add("HarmonicPowerFlow")
```

To run a harmonic power flow, the following objects are required:

- a `PowerGrid` struct, containing information on buses and lines which can either be imported or created manually,
- a `Settings` struct with simulation parameters and
- the Norton equivalent parameters of the nonlinear devices connected to the grid.

An example of a simple harmonic power flow can be found [here](https://github.com/pweigmann/HarmonicPowerFlow.jl/tree/main/examples/simple_coupled_hpf).


## Module index
```@index
Modules = [HarmonicPowerFlow]
Order   = [:constant, :type, :function, :macro]
```