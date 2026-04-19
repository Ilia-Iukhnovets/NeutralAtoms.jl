# Release And Recapture

Release-and-recapture is the package's temperature-estimation entry point for a
single atom in a Gaussian optical tweezer.

```@example recapture
using NeutralAtoms
using Plots
using Random

Random.seed!(1234)

trap_params = [1000.0, 1.1, w0_to_z0(1.1, 0.852, 1.3)]
atom_params = [86.9091835, 35.0]
tspan = collect(0.0:2.0:40.0)

recapture, acc_rate = release_recapture(
    tspan,
    trap_params,
    atom_params,
    200;
    harmonic = true,
)

plt = plot(tspan, recapture; lw = 3, label = "recapture")
xlabel!(plt, "Release time (μs)")
ylabel!(plt, "Recapture probability")
title!(plt, "Release-and-recapture curve")
plt
```

`harmonic = true` is the fast approximation used in the documentation and for
quick scans. For full-temperature studies in the Gaussian trap potential, use
`harmonic = false`.

