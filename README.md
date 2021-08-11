# HarmonicPowerFlow.jl


[![codecov](https://codecov.io/gh/pweigmann/HarmonicPowerFlow.jl/branch/master/graph/badge.svg?token=7DZTYSH7TY)](https://codecov.io/gh/pweigmann/HarmonicPowerFlow.jl)
[![CI](https://github.com/pweigmann/HarmonicPowerFlow.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/pweigmann/HarmonicPowerFlow.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pweigmann.github.io/HarmonicPowerFlow.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pweigmann.github.io/HarmonicPowerFlow.jl/dev)

An open source module to solve the harmonic power flow problem using Norton representation of nonlinear devices.

Feel free to contribute and get in contact and check out the [documentation](https://pweigmann.github.io/HarmonicPowerFlow.jl/dev), which is still a work in progess.

## Disclaimer

In the current state, the module does not produce trustworthy results! Further tweaks and additional validation is required. Please check out the issues and feel free to add ideas and suggestions. Contributions are highly welcome!

## Future Features

The following list gives some ideas about which additions would be most beneficial to increase usability:

- Modeling transformers
- Support for nonlinear power sources (PV buses), so far only linear sources
- Adding more devices to the library of Norton equivalent parameters
- Alternative solvers for cases when Newton-Raphson doesn't converge
- Load aggregation
- ...