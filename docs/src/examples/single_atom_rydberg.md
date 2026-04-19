# Single-Atom Two-Photon Rydberg Dynamics

The single-atom simulation path combines thermal motion, Doppler shifts,
optional phase noise, and spontaneous decay in the style of
[arXiv:1802.10424](https://arxiv.org/abs/1802.10424).

```@example rydberg
using NeutralAtoms
using QuantumOptics
using Plots
using Random

Random.seed!(1234)

atom_params = [86.9091835, 20.0]
trap_params = [1000.0, 1.1, w0_to_z0(1.1, 0.813, 1.3)]

Ωr = 2π * 8.0
Ωb = 2π * 8.0
Δ0 = 2π * 250.0

first_laser_params = Dict(
    "Ω" => Ωr,
    "w0" => 10.0,
    "z0" => w0_to_z0(10.0, 0.795),
    "θ" => 0.0,
    "n_sg" => 1,
    "type" => "gauss",
)
second_laser_params = Dict(
    "Ω" => Ωb,
    "w0" => 3.5,
    "z0" => w0_to_z0(3.5, 0.475),
    "θ" => 0.0,
    "n_sg" => 1,
    "type" => "gauss",
)

f = collect(0.02:0.02:0.10)
phase_amplitudes = zeros(length(f))
detuning_params = [Δ0, -δ_twophoton(Ωr, Ωb, Δ0)]
decay_params = [2π * 0.2, 2π * 0.2, 0.0, 0.0]
error_options = Dict(
    "laser_noise" => false,
    "spontaneous_decay_intermediate" => true,
    "spontaneous_decay_rydberg" => false,
    "atom_motion" => true,
    "free_motion" => true,
    "xy_motion" => true,
    "z_motion" => true,
    "Doppler" => true,
)

tspan = collect(0.0:0.05:0.35)

cfg = RydbergConfig(
    tspan,
    ket_1,
    atom_params,
    trap_params,
    2,
    f,
    phase_amplitudes,
    phase_amplitudes,
    first_laser_params,
    second_laser_params,
    [0.0, 0.0, 0.0],
    detuning_params,
    decay_params,
    error_options,
)

ρ, ρ2 = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        simulation(cfg; reltol = 1e-6, abstol = 1e-8)
    end
end

rydberg_probs = get_rydberg_probs(ρ, ρ2)
plot_rydberg_probs(tspan, rydberg_probs)
```

Key helpers in this workflow:

- `samples_generate`, `R`, and `V` encode finite-temperature motion.
- `Δ` and `δ` convert that motion into one- and two-photon Doppler detunings.
- `Sϕ`, `ϕ_amplitudes`, and `ϕ` construct stochastic phase-noise trajectories.
- `simulation` averages the master-equation dynamics over those realizations.

