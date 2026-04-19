# API Reference

```@meta
CurrentModule = NeutralAtoms
```

## Package

```@docs
NeutralAtoms
```

## Beam And Trap Optics

```@docs
w0_to_z0
trap_frequencies
E
I
get_trap_params
```

## Thermal Sampling And Release/Recapture

```@docs
H
samples_generate
R
V
release_recapture
```

## Laser Phase Noise

```@docs
Sϕ
ϕ_amplitudes
ϕ
```

## Single-Atom Rydberg Simulation

```@docs
RydbergConfig
ket_0
ket_1
ket_r
ket_p
ket_l
gauss_field
simulation
Ω_twophoton
T_twophoton
δ_twophoton
Ωr_required
get_rydberg_probs
plot_rydberg_probs
```

## Beam Shaping Helpers

```@docs
simple_flattopHG_field
simple_flattopLG_field
HG_coeff
HG_coefficients
decomposition_HG_2d
reconstruct_HG_field_2d
```

## Two-Qubit CZ Workflow

```@docs
CZLPConfig
simulation_czlp
get_two_qubit_probs
plot_two_qubit_probs
get_fidelity_with_rz_phi
CZ_calibration_by_fidelity_oscillation
get_parity_osc
get_rydberg_fidelity_configs
get_rydberg_infidelity
get_cz_infidelity
```

## Gate Utilities

```@docs
get_gate
project_on_qubit
```
