using Distributed 
using NeutralAtoms  
using JLD2

num_proc = 8
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../conf/configs_gen.jl") # include("../conf/default.jl")

@everywhere function compute_ρ(i)
    _, cfg_CZ = get_6P_config();
    
    #cfg_CZ.atom_params[2] = 70.0; #temperature #
    cfg_CZ.error_options["z_motion"] = false #true
    cfg_CZ.n_samples = 50 
    
    ρ_end = NeutralAtoms.simulation_czlp(cfg_CZ)[1][end]
    return ρ_end 
end 

function main()
    ρ_computed = pmap(compute_ρ, 1:num_proc)

    ρ = sum(ρ_computed) ./ num_proc
    
    save_QO_operator("data\\rho_distributed.jld2", ρ)

    _, cfg_CZ = get_6P_config();
    ϕs =  [0.0:0.001:2π;];
    Fids = [NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg_CZ.ψ0, ϕ) for ϕ in ϕs];
    println(maximum(Fids), " ", ϕs[argmax(Fids)])
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())