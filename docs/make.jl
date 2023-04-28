using GridHelpers
using Documenter

DocMeta.setdocmeta!(GridHelpers, :DocTestSetup, :(using GridHelpers); recursive=true)

makedocs(;
    modules=[GridHelpers],
    authors="Cédric BELMANT",
    repo="https://github.com/serenity4/GridHelpers.jl/blob/{commit}{path}#{line}",
    sitename="GridHelpers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://serenity4.github.io/GridHelpers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/serenity4/GridHelpers.jl",
    devbranch="main",
)
