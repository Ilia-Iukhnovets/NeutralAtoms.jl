using NeutralAtoms
using QuantumOptics
using Serialization
using DataFrames
using CSV
using DelimitedFiles 

function get_6P_config()
       m = 86.9091835;            # Atom params
       T = 1.0 # 70.0;
       atom_params = [m, T];

       U0 = 800 * 0.65; #797.7 changed 26.02.2026
       w0 = 1.42; # Trap params
       λ0 = 0.813;
       M2 = 1.0;
       z0 = w0_to_z0(w0, λ0, M2);
       trap_params = [U0, w0, z0];

       h0 = 13.0 * 1e-6;    #MHz^2/MHz #Params for laser phase noise
       hg1 = 25.0 * 1e-6;   #MHz^2/MHz
       hg2 = 10.0e3 * 1e-6; #MHz^2/MHz
       fg1 = 130.0 * 1e-3;  #MHz
       fg2 = 234.0 * 1e-3;  #MHz
       σg1 = 18.0 * 1e-3;   #MHz
       σg2 = 1.5 * 1e-3;    #MHz
       first_laser_phase_params  = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       second_laser_phase_params = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       f = [0.01:0.0025:1.0;];
       first_laser_phase_amplitudes = ϕ_amplitudes(f, first_laser_phase_params);
       second_laser_phase_amplitudes = ϕ_amplitudes(f, second_laser_phase_params);

       λ2 = 1.012;        ### Excitation beam parameters
       λ1 = 0.42;
       w2 = 20.0;
       w1 = 8.0 #10.0;
       z2 = w0_to_z0(w2, λ2);
       z1 = w0_to_z0(w1, λ1);

       Δ0 = 2.0*π * 2000 #870 #1600.0 #1000.0 #* 1600
       Ω = 2π * 2.3 # 2.0;
       a = 1 
       Ω1 = 1/a * sqrt(2* Δ0 * Ω);
       Ω2 = a * sqrt(2* Δ0 * Ω);
       θ2 = π;
       θ1 = 0.0;
       n2 = 1;
       n1 = 1;

       first_laser_params = Dict("Ω" => Ω2,"w0" => w2,"z0" => z2,
              "θ" => θ2,"n_sg" => n2,"type" => "gauss") #first_laser_params = [Ωr, wr, zr, θr] #, nr];
       second_laser_params = Dict("Ω" => Ω1,"w0" => w1,"z0" => z1,
              "θ" => θ1,"n_sg" => n1,"type" => "gauss")      #second_laser_params = [Ωb, wb, zb, θb] #, nb];

       detuning_params = [Δ0, -δ_twophoton(Ω2, Ω1, Δ0)]; #???

       Γ =  8.475 #2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       τr = 445.74; 
       Γr = 1/τr;
       decay_params = [Γ0, Γ1, Γl, Γr];
       
       #T0 = T_twophoton(Ω2, Ω1, Δ0) # Simulation params
       #tspan = [0.0:T0/20:3.0;]; #15.0
       T0 = 2.5 #T_twophoton(Ω1, Ω2, Δ0)
       tspan = [0.0:T0/200:T0;];
       ψ0 = ket_1;
       n_samples = 20;
       shift = [0.0,0.0,0.0]

       error_options = Dict("laser_noise" => false,
                        "spontaneous_decay_intermediate" => true,
                        "spontaneous_decay_rydberg" => true,
                        "atom_motion" => true,
                        "free_motion" => true,
                        "xy_motion" => true,
                        "z_motion" => true,
                        "Doppler" => true
                        )

       cfg = NeutralAtoms.RydbergConfig(
              tspan,
              ψ0,

              atom_params,
              trap_params,
              n_samples,

              f,
              first_laser_phase_amplitudes,
              second_laser_phase_amplitudes,

              first_laser_params,
              second_laser_params,
              shift,

              detuning_params,
              decay_params,

              error_options
              );
 
       d = 2.7 #2.0;
       #atom_centers = [[-d/2, 0.0, 0.0], [d/2, 0.0, 0.0]]
       atom_centers = [[0.0, 0.0,-d/2], [0.0, 0.0,d/2]]
       c6 = 2π * 135298
       
       ΔtoΩ = 0.377371 #идеальные
       Ωτ = 4.29268
       ξ = 3.90242
       
       ket_pos = (ket_0 + ket_1)/sqrt(2)
       ψ0_cz = ket_pos ⊗ ket_pos
       
       Ω_twophoton = (2π/T_twophoton(Ω1, Ω2, Δ0))
       τ = 2π / (Ω_twophoton * sqrt(ΔtoΩ^2 + 2.0))

       #Ω_twophoton = (2π/T0)
       #τ = Ωτ/Ω_twophoton #τ = 2π / (Ω_twophoton * sqrt(ΔtoΩ^2 + 2.0))
       VV = c6 / d^6 #3.4^6 #       Ω_twophoton / 2 / V
       detuning_params = [Δ0, ΔtoΩ * Ω_twophoton - δ_twophoton(Ω1, Ω2, Δ0)]; #ΔtoΩ = ΔtoΩ + Ω_twophoton/(2*VV)
       
       ϕ2 = 2*τ * ΔtoΩ * Ω_twophoton;
       ϕ1 = (ϕ2 - π)/2  
       tspan_cz = [0.0, 2*τ];        #detuning_params = [Δ0, ΔtoΩ*Ω_twophoton + δ_twophoton(Ωr, Ωb, Δ0)];
       
       cfg_czlp = NeutralAtoms.CZLPConfig(
              tspan_cz,
              ψ0_cz,

              atom_params,
              trap_params,
              n_samples,

              f,
              first_laser_phase_amplitudes,
              second_laser_phase_amplitudes,

              first_laser_params,
              second_laser_params,

              detuning_params,
              decay_params,

              error_options,

              atom_centers,
              c6,
              ΔtoΩ,
              Ωτ,
              ξ,
              ϕ1
       )
       return cfg, cfg_czlp
end

function get_5P_config()
       _, cfg_CZ_5P = get_default_config()
       ket_pos = (ket_0 + ket_1) / sqrt(2)
       cfg_CZ_5P.ψ0 = ket_pos ⊗ ket_pos 

       w = 2. ;
       bl_lsr_prms = cfg_CZ_5P.second_laser_params
       Ω0, w0, z0 = bl_lsr_prms["Ω"], bl_lsr_prms["w0"], bl_lsr_prms["z0"];
       z = z0*(w/w0)^2        #cfg_CZ_6P.second_laser_params["type"] = "gauss" 
       cfg_CZ_5P.second_laser_params["w0"] = w ; 
       cfg_CZ_5P.second_laser_params["z"] = z;

       Γ = 2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 111.74; 
       Γr = 1/τr;
       cfg_CZ_5P.decay_params = [Γ0, Γ1, Γl, Γr];

       return cfg_CZ_5P
end;