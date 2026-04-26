using Distributed 
using NeutralAtoms  

include("../conf/configs_gen.jl") 
_, cfg_CZ = get_6P_config();
cfg_CZ.n_samples=10
configs = get_rydberg_fidelity_configs(cfg_CZ, cfg_CZ.n_samples; int_prob=true) 
configs_list = [[key, configs[key]] for key in configs.keys];

num_proc = length(configs_list)
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../conf/configs_gen.jl") # include("../conf/default.jl")

@everywhere function compute_ρ(cnfg)
    cfg = deepcopy(cnfg[2])

    ρ = NeutralAtoms.simulation_czlp(cfg)[1][end]
    fid = NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg.ψ0, cfg.ϕ_RZ)
    
    return [cnfg[1], 1-fid]
end 

function main()
    err_computed = pmap(compute_ρ, configs_list) #1:num_proc)
    
    budget = Dict()
    for err in err_computed
        budget[err[1]] = err[2]
    end

    for key in keys(budget)
        println(key, " ", budget[key])
    end
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())