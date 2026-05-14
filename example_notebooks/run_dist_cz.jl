using Distributed 
using NeutralAtoms  
using JLD2

num_proc = 5 #14
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../config5P_EZ.jl") # include("../conf/default.jl")

@everywhere function compute_ρ(i)
    cfg_rd, cfg_CZ = get_5P_config();

    cfg_CZ.n_samples = 100 #40
    #cfg_CZ.error_options["atom_motion"] = false # true # 

    #ρ = NeutralAtoms.simulation(cfg_rd)[1] 
    ρ = NeutralAtoms.simulation_czlp(cfg_CZ)[1][end]
    return ρ
end 

function main()
    ρ_computed = pmap(compute_ρ, 1:num_proc)

    ρ_mean = sum(ρ_computed) ./ num_proc
    
    #Pr_r = real(expect(ket_r ⊗ dagger(ket_r), ρ_mean));  
    NeutralAtoms.save_QO_operator("data\\rho_distr.jld2", ρ_mean)

end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())