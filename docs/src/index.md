# NeutralAtoms.jl

[![Build Status](https://github.com/mgoloshchapov/NeutralAtoms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mgoloshchapov/NeutralAtoms.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://mgoloshchapov.github.io/NeutralAtoms.jl/dev/)

`NeutralAtoms.jl` simulates neutral-atom experiments in optical tweezers, with a
focus on two-photon Rydberg excitation and blockade-mediated controlled-phase
gates.

The package is organized around the same experimental workflow as the papers
that motivate it:

- trap characterization and thermal sampling of atom motion
- release-and-recapture thermometry
- effective two-photon Rydberg excitation with Doppler shifts
- stochastic laser phase noise and spontaneous decay channels
- blockade-mediated two-qubit phase-gate simulation and calibration

## Installation

Install the package directly from GitHub:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/mgoloshchapov/NeutralAtoms.jl"))
```

## Package Scope

The current API models:

- Gaussian tweezer optics and trap frequencies
- Monte Carlo sampling of thermal positions and velocities
- time-dependent master-equation simulation for a single atom
- configurable beam models, including Gaussian and Hermite-Gauss based profiles
- two-atom blockade simulations for CZ-style global-pulse protocols
- helper routines for parity scans and infidelity budgets

The single-atom modeling follows the physical picture emphasized in
[de Léséleuc et al. (2018)](https://arxiv.org/abs/1802.10424): Doppler
dephasing, spontaneous emission through the intermediate state, and laser phase
noise. The two-qubit modeling follows the blockade-based global-pulse gate
protocol of [Levine et al. (2019)](https://arxiv.org/abs/1908.06101).

## Start Here

- Read [Examples](examples.md) for tutorial-style workflows using the current
  dict-based laser parameter API.
- Use [API](library.md) for the exported functions and configuration types.

## References

- Sylvain de Léséleuc, Daniel Barredo, Vincent Lienhard, Antoine Browaeys, and
  Thierry Lahaye, [“Analysis of imperfections in the coherent optical excitation
  of single atoms to Rydberg states”](https://arxiv.org/abs/1802.10424).
- Harry Levine, Alexander Keesling, Giulia Semeghini, Ahmed Omran, Tout T.
  Wang, Sepehr Ebadi, Hannes Bernien, Markus Greiner, Vladan Vuletić, Hannes
  Pichler, and Mikhail D. Lukin, [“Parallel implementation of high-fidelity
  multi-qubit gates with neutral atoms”](https://arxiv.org/abs/1908.06101).
