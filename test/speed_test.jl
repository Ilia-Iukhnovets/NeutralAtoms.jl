using NeutralAtoms 

include("../conf/configs_gen.jl")
cfg_CZ = get_6P_config();
cfg_CZ.n_samples = 1;

# Manual benchmark: both runs should finish in the same rough time range.
# `xy_motion=false` should remove only thermal x/y motion while keeping
# the static trap-center separation intact.
p = simulation_czlp(cfg_CZ)

cfg_CZ.error_options["xy_motion"] = false
p = simulation_czlp(cfg_CZ)
