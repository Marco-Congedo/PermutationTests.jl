push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, @__DIR__)

using Documenter
using PermutationTests

makedocs(
    #remotes = nothing, # ELIMINATE for deploying
    sitename = "PermutationTests",
    #format = Documenter.HTML((prettyurls = false)), # ELIMINATE pretty URL for deploying
    format = Documenter.HTML(),
    authors="Marco Congedo, CNRS, Grenoble, France; Livio Finos, Uni. Padova, Italia",
    modules = [PermutationTests],
    pages =  [
        "index.md",
        "About" => "about.md",
        "Main Module" => "PermutationTests.md",
        "Univariate tests" => "univariate tests.md",
        "Multiple comparisons tests" => "multiple comparisons tests.md",
        "Package tests" => "package tests.md",
        "Tools" => "tools.md",
        "Extras" => Any[
            "Statistics" => "statistics.md",
            "Test statistics" => "test statistics.md",
            "p-value combinations" => "pCombination.md",
            "Create your own test" => "create your own test.md",
            "Chose a Test" => "chose a test.md"
        ]
    ]
)

  
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    # root,
   # target = "build", # add this folder to .gitignore!
   repo = "github.com/Marco-Congedo/PermutationTests.jl.git",
   branch = "gh-pages",
   # osname = "linux",
   # deps = Deps.pip("pygments", "mkdocs"),
   devbranch = "dev",
   devurl = "dev",
   push_preview = true,
   # versions = ["stable" => "v^", "v#.#", devurl => devurl]
)
