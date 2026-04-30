using Documenter, CalciumSpike

pages = [
    "Home"             => "index.md",
    "Forward Model"    => "forward_model.md",
    "Post-Processing"  => "postprocessing.md",
    "Model Comparison" => "comparison.md",
    "Noise Correction" => "noise_correction.md",
    "MLSpike"          => "mlspike.md",
    "Visualization"    => "visualization.md",
    # "API Reference"    => "api_reference.md",
]

makedocs(
    sitename = "CalciumSpike.jl",
    modules  = [CalciumSpike],
    warnonly = [:autodocs_block, :missing_docs, :cross_references],
    format   = Documenter.HTML(),
    pages    = pages,
)

deploydocs(repo = "github.com/aquaresi/CalciumSpike.jl.git")
