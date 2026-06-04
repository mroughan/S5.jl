using Documenter
using S5

makedocs(
    sitename = "S5.jl",
    modules  = [S5],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = nothing,
    ),
    pages = [
        "Home"          => "index.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)

# Uncomment and configure to deploy to GitHub Pages:
# deploydocs(
#     repo = "github.com/mroughan/S5.jl.git",
# )
