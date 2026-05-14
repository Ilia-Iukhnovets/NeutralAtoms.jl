using Distributed 
using NeutralAtoms  
using JLD2

num_proc = 5 #14
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../config5P_Flattop.jl") # include("../conf/default.jl")

@everywhere function compute_ρ(i)
    cfg_rd, cfg_CZ = get_5P_config();

    cfg_rd.n_samples = 40 #40

    ρ = NeutralAtoms.simulation(cfg_rd)[1] 
    
    #ρ = NeutralAtoms.simulation_czlp(cfg_CZ)[1][end]
    return ρ
end 

function main()
    ρ_computed = pmap(compute_ρ, 1:num_proc)

    ρ_mean = sum(ρ_computed) ./ num_proc
    
    Pr_r = real(expect(ket_r ⊗ dagger(ket_r), ρ_mean));  

    NeutralAtoms.save_QO_operator("data\\rydberg_rabi_1.jld2", Pr_r)

    #_, cfg_CZ = get_6P_config();
    #ϕs =  [0.0:0.001:2π;];
    #Fids = [NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg_CZ.ψ0, ϕ) for ϕ in ϕs];
    #println(maximum(Fids), " ", ϕs[argmax(Fids)])
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())