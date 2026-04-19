# Two-Qubit Blockade Simulation

`simulation_czlp` extends the same ingredients to a two-atom blockade model in
the spirit of [arXiv:1908.06101](https://arxiv.org/abs/1908.06101).

```@example cz
using NeutralAtoms
using QuantumOptics
using Plots
using Random

Random.seed!(1234)

atom_params = [86.9091835, 10.0]
trap_params = [1000.0, 1.1, w0_to_z0(1.1, 0.813, 1.3)]

Ωr = 2π * 6.0
Ωb = 2π * 6.0
Δ0 = 2π * 200.0

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
decay_params = [2π * 0.1, 2π * 0.1, 0.0, 0.0]
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

ket_plus = (ket_0 + ket_1) / sqrt(2)
tspan = [0.0, 0.05]

cfg = CZLPConfig(
    tspan,
    ket_plus ⊗ ket_plus,
    atom_params,
    trap_params,
    1,
    f,
    phase_amplitudes,
    phase_amplitudes,
    first_laser_params,
    second_laser_params,
    detuning_params,
    decay_params,
    error_options,
    [[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    2π * 500.0 * 4.0^6,
    1.0,
    1.0,
    π / 2,
    0.0,
)

ρ, ρ2 = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        simulation_czlp(cfg; reltol = 1e-6, abstol = 1e-8)
    end
end

two_qubit_probs = get_two_qubit_probs(ρ, ρ2)
plot_two_qubit_probs(tspan, two_qubit_probs)
```

For phase-gate calibration and fidelity analysis, the main helpers are
`get_fidelity_with_rz_phi`, `CZ_calibration_by_fidelity_oscillation`,
`get_parity_osc`, and `get_cz_infidelity`.

