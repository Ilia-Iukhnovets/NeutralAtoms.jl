using Test
using NeutralAtoms
using QuantumOptics

@testset "Analytic helpers" begin
    @test w0_to_z0(1.0, 0.5) ≈ 2π

    ωr, ωz = trap_frequencies([86.9, 10.0], [500.0, 1.2, w0_to_z0(1.2, 0.813)])
    @test ωr > 0
    @test ωz > 0

    @test Ω_twophoton(4.0, 2.0, 8.0) ≈ 0.5
    @test δ_twophoton(4.0, 2.0, 8.0) ≈ 0.375
    @test Ωr_required(0.5, 2.0, 8.0) ≈ 4.0
    @test gauss_field(0.0, 0.0, 0.0, 1.0, 2.0) ≈ 1.0 + 0.0im
end

@testset "Configuration types" begin
    tspan = [0.0, 1e-3]
    ψ0 = ket_1

    atom_params = [86.9, 10.0]
    trap_params = [500.0, 1.2, w0_to_z0(1.2, 0.813)]
    f = [0.1, 0.2]
    phase_amplitudes = zeros(2)

    first_laser_params = Dict(
        "Ω" => 1.0,
        "w0" => 10.0,
        "z0" => w0_to_z0(10.0, 0.795),
        "θ" => 0.0,
        "n_sg" => 1,
        "type" => "gauss",
    )
    second_laser_params = Dict(
        "Ω" => 1.0,
        "w0" => 10.0,
        "z0" => w0_to_z0(10.0, 0.475),
        "θ" => 0.0,
        "n_sg" => 1,
        "type" => "gauss",
    )

    shift = [0.0, 0.0, 0.0]
    detuning_params = [2π * 1000, 0.0]
    decay_params = zeros(4)
    error_options = Dict(
        "laser_noise" => false,
        "spontaneous_decay_intermediate" => false,
        "spontaneous_decay_rydberg" => false,
        "atom_motion" => false,
        "free_motion" => false,
        "xy_motion" => false,
        "z_motion" => false,
        "Doppler" => false,
    )

    cfg = RydbergConfig(
        tspan,
        ψ0,
        atom_params,
        trap_params,
        1,
        f,
        phase_amplitudes,
        phase_amplitudes,
        first_laser_params,
        second_laser_params,
        shift,
        detuning_params,
        decay_params,
        error_options,
    )
    @test cfg.first_laser_params["Ω"] == 1.0
    @test cfg.error_options["laser_noise"] === false

    cfg_cz = CZLPConfig(
        tspan,
        ket_0 ⊗ ket_0,
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
        [[0.0, -1.5, 0.0], [0.0, 1.5, 0.0]],
        1.0,
        0.25,
        1.0,
        0.0,
        0.0,
    )
    @test length(cfg_cz.atom_centers) == 2
    @test cfg_cz.ξ == 0.0
end

@testset "Single-atom simulation smoke test" begin
    tspan = [0.0, 1e-3]
    atom_params = [86.9, 10.0]
    trap_params = [500.0, 1.2, w0_to_z0(1.2, 0.813)]

    first_laser_params = Dict(
        "Ω" => 1.0,
        "w0" => 10.0,
        "z0" => w0_to_z0(10.0, 0.795),
        "θ" => 0.0,
        "n_sg" => 1,
        "type" => "gauss",
    )
    second_laser_params = Dict(
        "Ω" => 1.0,
        "w0" => 10.0,
        "z0" => w0_to_z0(10.0, 0.475),
        "θ" => 0.0,
        "n_sg" => 1,
        "type" => "gauss",
    )
    error_options = Dict(
        "laser_noise" => false,
        "spontaneous_decay_intermediate" => false,
        "spontaneous_decay_rydberg" => false,
        "atom_motion" => false,
        "free_motion" => false,
        "xy_motion" => false,
        "z_motion" => false,
        "Doppler" => false,
    )

    cfg = RydbergConfig(
        tspan,
        ket_1,
        atom_params,
        trap_params,
        1,
        [0.1, 0.2],
        zeros(2),
        zeros(2),
        first_laser_params,
        second_laser_params,
        [0.0, 0.0, 0.0],
        [2π * 1000, 0.0],
        zeros(4),
        error_options,
    )

    ρ, ρ2 = simulation(cfg)
    @test length(ρ) == length(tspan)
    @test length(ρ2) == length(tspan)
    @test real(expect(ket_1 ⊗ dagger(ket_1), ρ[end])) ≤ 1.0
    @test real(expect(ket_1 ⊗ dagger(ket_1), ρ[end])) ≥ 0.0
end

include("motion_regressions.jl")
