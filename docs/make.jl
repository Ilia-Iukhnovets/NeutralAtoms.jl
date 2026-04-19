using Documenter, NeutralAtoms

push!(LOAD_PATH,"../src/")

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    edit_link = "main",
    assets = [joinpath("assets", "logo.ico")],
    size_threshold_ignore = ["library.md"]
)

makedocs(
    sitename="NeutralAtoms.jl",
    modules=[NeutralAtoms],
    checkdocs=:exports,
    format=format,
    authors="M.Y. Goloshchapov",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Overview" => "examples/index.md",
            "Release And Recapture" => "examples/release_recapture.md",
            "Single-Atom Two-Photon Rydberg Dynamics" => "examples/single_atom_rydberg.md",
            "Two-Qubit Blockade Simulation" => "examples/two_qubit_cz.md",
        ],
        "API" => "library.md"
    ]
    )

deploydocs(
    repo = "https://github.com/mgoloshchapov/NeutralAtoms.jl",
    devbranch = "main",
)
