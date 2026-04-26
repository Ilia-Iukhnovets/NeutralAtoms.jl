"""
    w(z, w0, z0)

Return the Gaussian beam radius at longitudinal coordinate `z`.
"""
function w(z, w0, z0)
    return w0 .* sqrt.(1.0 .+ (z ./z0) .^2);
end;

"""
    w0_to_z0(w0, λ, M2=1.0)

Return the Rayleigh length associated with a beam waist `w0`.

The conversion is

```math
z_0 = \\frac{\\pi w_0^2}{\\lambda M^2}.
```

# Arguments
- `w0`: beam waist radius in `μm`.
- `λ`: beam wavelength in `μm`.
- `M2 = 1.0`: beam-quality factor. `M2 = 1` corresponds to an ideal Gaussian beam.

# Returns
- Rayleigh length in `μm`.
"""
function w0_to_z0(w0, λ, M2=1.0)
    return π*w0^2/λ / M2;
end;

"""
    A(x, y, z, w0, z0; n=1, θ=0.0)

Return the normalized real envelope of the generalized Gaussian beam model.
"""
function A(x, y, z, w0, z0; n=1, θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ)-zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    return (w0 ./ w(zt, w0, z0)) .* exp.(- ((xt .^2 .+ yt .^2) ./ (w(zt, w0, z0) .^2)) .^ n)
end;

"""
    I(x, y, z, w0, z0; n=1, θ=0.0)

Return the normalized intensity profile associated with `A`.
"""
function I(x, y, z, w0, z0; n=1, θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ) - zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    return ((w0 ./ w(zt, w0, z0)) .* exp.(-((xt .^2 .+ yt .^2) ./ (w(zt, w0, z0) .^2)).^n)) .^2
end;

"""
    A_phase(x, y, z, w0, z0; θ=0.0)

Return the complex propagation phase of the Gaussian beam model, including the
wavefront-curvature and Gouy-phase terms used by the excitation model.
"""
function A_phase(x, y, z, w0, z0; θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ)-zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    k = 2.0 * z0 / w0^2;
    #return exp.(-1.0im * k * zt .* (0.5*(xt .^2 + yt .^ 2) ./ (zt .^2  + z0 .^2)) + 1.0im * atan.(zt ./ z0));
    return exp.(1.0im .* ( atan.(zt ./ z0) .- k * zt .*(1.0 .+ 0.5 .* (xt .^2 + yt .^ 2) ./ (zt .^2  + z0 .^2) ) ) );
end;

"""
    E(x, y, z, w0, z0; n=1, θ=0.0)

Return the normalized complex field profile `A(...) * A_phase(...)`.
"""
function E(x, y, z, w0, z0;n=1, θ=0.0)
    return A(x,y,z,w0,z0;n=n, θ=θ) .* A_phase(x,y,z,w0,z0; θ=θ)
end;

"""
    trap_frequencies(atom_params, trap_params)

Return the harmonic trap frequencies for a Gaussian optical tweezer.

# Arguments
- `atom_params`: `[m, T]`, with mass `m` in atomic mass units and temperature
  `T` in `μK`.
- `trap_params`: `[U0, w0, z0]`, with trap depth `U0` in `μK`, waist `w0` in
  `μm`, and Rayleigh length `z0` in `μm`.

# Returns
- `(ωr, ωz)`: radial and longitudinal trap frequencies in `MHz`.
"""
function trap_frequencies(atom_params, trap_params)
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    ωr = vconst/w0 * sqrt(4 * U0/m);
    ωz = vconst/z0 * sqrt(2 * U0/m);
    
    return ωr, ωz;
end;

"""
    get_rydberg_probs(ρ, ρ2, eps=1e-12)

Extract single-atom level populations and estimated sampling errors from the
first and second moments returned by `simulation`.
"""
function get_rydberg_probs(ρ, ρ2, eps=1e-12)
    probs_dict = OrderedCollections.OrderedDict{String, Vector{Float64}}();

    names = ["0", "1", "r", "p", "l"];
    states = [ket_0, ket_1, ket_r, ket_p, ket_l];
    for i in 1:5
        P = real(expect(states[i] ⊗ dagger(states[i]), ρ))
        P2 = real(expect(states[i] ⊗ dagger(states[i]), ρ2))
        S = @. sqrt(P2 - P^2 .+ eps) / length(ρ)
        probs_dict["P"*names[i]] = P
        probs_dict["S"*names[i]] = S 
    end 

    return probs_dict
end

"""
    get_two_qubit_probs(ρ, ρ2, eps=1e-12)

Extract computational-basis populations and estimated sampling errors from the
first and second moments returned by `simulation_czlp`.
"""
function get_two_qubit_probs(ρ, ρ2, eps=1e-12)
    probs_dict = OrderedCollections.OrderedDict{String, Vector{Float64}}();
    names = ["00", "01", "10", "11"];

    states = [
        ket_0 ⊗ ket_0, 
        ket_0 ⊗ ket_1, 
        ket_1 ⊗ ket_0, 
        ket_1 ⊗ ket_1
        ];
    for i in 1:4
        P = real(expect(states[i] ⊗ dagger(states[i]), ρ))
        P2 = real(expect(states[i] ⊗ dagger(states[i]), ρ2))
        S = @. sqrt(P2 - P^2 .+ eps) / length(ρ)
        probs_dict["P"*names[i]] = P
        probs_dict["S"*names[i]] = S 
    end 

    return probs_dict
end

"""
    plot_rydberg_probs(tspan, probs_dict)

Plot the single-atom populations stored in `probs_dict`, typically produced by
`get_rydberg_probs`.
"""
function plot_rydberg_probs(tspan, probs_dict)
    names = ["0", "1", "r", "p", "l"];
    colors = ["lightblue", "blue", "red", "orange", "green"];

    plt = Plots.plot()
    for i in 1:5
        P = probs_dict["P"*names[i]]
        S = probs_dict["S"*names[i]]
        plot!(
            tspan, [P P], fillrange=[P+S P-S], 
            ylim=(0.0, 1.0), xlim=(minimum(tspan), maximum(tspan)), 
            fillalpha=0.25, c=colors[i], 
            label=[nothing "P" * names[i]], linewidth=3
            )
    end
    xlabel!("Time, μs")
    ylabel!("Probability")
    title!("Rydberg Rabi oscillations")

    return plt
end

"""
    plot_two_qubit_probs(tspan, probs_dict)

Plot the computational-basis populations stored in `probs_dict`, typically
produced by `get_two_qubit_probs`.
"""
function plot_two_qubit_probs(tspan, probs_dict)
    names = ["00", "01", "10", "11"];
    colors = Plots.cgrad(:bam, 4, categorical = true)

    plt = Plots.plot()
    for i in 1:4
        P = probs_dict["P"*names[i]]
        S = probs_dict["S"*names[i]]
        plot!(
            tspan, [P P], fillrange=[P+S P-S], 
            ylim=(0.0, 1.0), xlim=(minimum(tspan), maximum(tspan)), 
            fillalpha=0.25, c=colors[i], 
            label=[nothing "P" * names[i]], linewidth=3
            )
    end
    xlabel!("Time, μs")
    ylabel!("Probability")
    title!("Two-qubit probabilities")

    return plt
end

"""
    RydbergConfig

Configuration for single-atom two-photon Rydberg simulations.

This type bundles the ingredients used by the single-atom model inspired by
[arXiv:1802.10424](https://arxiv.org/abs/1802.10424): thermal atom motion,
laser phase noise, detuning, and spontaneous decay channels.

# Fields
- `tspan::Vector{Float64}`: time points in `μs` used for solver output.
- `ψ0`: initial pure state in the five-level basis.
- `atom_params::Vector{Float64}`: `[m, T]`, with mass in atomic mass units and
  temperature in `μK`.
- `trap_params::Vector{Float64}`: `[U0, w0, z0]`, with `U0` in `μK` and lengths
  in `μm`.
- `n_samples::Int64`: number of Monte Carlo trajectories to average.
- `f::Vector{Float64}`: phase-noise frequency grid in `MHz`.
- `first_laser_phase_amplitudes::Vector{Float64}`: Fourier amplitudes for the
  first-laser phase-noise trace.
- `second_laser_phase_amplitudes::Vector{Float64}`: Fourier amplitudes for the
  second-laser phase-noise trace.
- `first_laser_params::Dict{String, Any}`: dict-based description of the first
  laser beam and coupling.
- `second_laser_params::Dict{String, Any}`: dict-based description of the second
  laser beam and coupling.
- `shift::Vector{Float64}`: static spatial shift applied to the sampled atom
  position before simulation.
- `detuning_params::Vector{Float64}`: `[Δ0, δ0]` single- and two-photon
  detunings in angular-frequency units.
- `decay_params::Vector{Float64}`: spontaneous decay rates from the intermediate
  and Rydberg states.
- `error_options::Dict{String, Any}`: switches controlling motion, Doppler
  shifts, phase noise, and spontaneous-decay channels.
"""
mutable struct RydbergConfig
    tspan::Vector{Float64}
    ψ0::Ket{NLevelBasis{Int64}, Vector{ComplexF64}}

    atom_params::Vector{Float64}
    trap_params::Vector{Float64}
    n_samples::Int64

    f::Vector{Float64}
    first_laser_phase_amplitudes::Vector{Float64}
    second_laser_phase_amplitudes::Vector{Float64}
    
    first_laser_params::Dict{String, Any}#Vector{Float64}
    second_laser_params::Dict{String, Any}#Vector{Float64}
    shift::Vector{Float64}
    
    detuning_params::Vector{Float64}
    decay_params::Vector{Float64}

    error_options::Dict{String, Any}
end

"""
    CZLPConfig

Configuration for the blockade-mediated controlled-phase simulation.

This extends `RydbergConfig` with the atom-pair geometry and phase-gate
calibration parameters used by `simulation_czlp` and the fidelity-analysis
helpers. The intended workflow follows the global-pulse CZ protocol discussed in
[arXiv:1908.06101](https://arxiv.org/abs/1908.06101).

# Fields
- `tspan`, `ψ0`, `atom_params`, `trap_params`, `n_samples`, `f`,
  `first_laser_phase_amplitudes`, `second_laser_phase_amplitudes`,
  `first_laser_params`, `second_laser_params`, `detuning_params`,
  `decay_params`, `error_options`: same meaning as in `RydbergConfig`.
- `atom_centers::Vector{Vector{Float64}}`: equilibrium positions of the two
  traps in `μm`.
- `c6::Float64`: van der Waals interaction coefficient used for blockade.
- `ΔtoΩ::Float64`: detuning-to-Rabi ratio used when calibrating the CZ pulse.
- `Ωτ::Float64`: pulse-area parameter used by downstream calibration utilities.
- `ξ::Float64`: phase step between the two global Rydberg pulses.
- `ϕ_RZ::Float64`: single-qubit `RZ` compensation phase used in parity analysis.
"""
mutable struct CZLPConfig
    tspan::Vector{Float64}
    ψ0

    atom_params::Vector{Float64}
    trap_params::Vector{Float64}
    n_samples::Int64

    f::Vector{Float64}
    first_laser_phase_amplitudes::Vector{Float64}
    second_laser_phase_amplitudes::Vector{Float64}
    
    first_laser_params::Dict{String, Any}#Vector{Float64}
    second_laser_params::Dict{String, Any}#Vector{Float64}
    
    detuning_params::Vector{Float64}
    decay_params::Vector{Float64}
    
    error_options::Dict{String, Any}

    atom_centers::Vector{Vector{Float64}}
    c6::Float64
    ΔtoΩ::Float64
    Ωτ::Float64
    ξ::Float64
    ϕ_RZ::Float64
end

# function calibrate_rabi(trap_params, atom_params, laser_params; n_samples=1000)
#     Ω, w, z, θ, n = laser_params
#     samples = samples_generate(trap_params, atom_params, n_samples)[1]
#     Ω_samples = [A(sample[1], sample[2], sample[3], w, z; n=n, θ=θ) for sample in samples];
#     Ω2_samples = [A(sample[1], sample[2], sample[3], w, z; n=n, θ=θ)^2 for sample in samples];
#     factor = sum(Ω_samples) / length(Ω_samples)
#     factor2 = sqrt(sum(Ω_samples^2) / length(Ω_samples))
#     return factor, factor2
# end

using JLD2

function save_QO_operator(filename::String, op::Operator; name::String="operator")
    jldopen(filename, "w") do file
        file[name] = op
    end
    println("Operator saved to $filename")
end

function load_QO_operator(filename::String; name::String="operator")
    op = jldopen(filename, "r") do file
        file[name]
    end
    return op
end