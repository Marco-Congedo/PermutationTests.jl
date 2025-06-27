# Build the documenattion locally
# This si useful to check the documentation before deploying
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, @__DIR__)

using Documenter, DocumenterTools, DocumenterCitations, DocumenterInterLinks
using PermutationTests, Distributions

makedocs(
    remotes = nothing, # ELIMINATE for deploying
    sitename = "PermutationTests",
    format = Documenter.HTML((prettyurls = false)), # ELIMINATE pretty URL for deploying
    authors="Marco Congedo, CNRS, Grenoble, France; Livio Finos, Uni. Padova, Italia",
    modules = [PermutationTests],
    pages =  [
        "index.md",
        "Main Module" => "PermutationTests.md",
        "Univariate tests" => "univariate tests.md",
        "Multiple comparisons tests" => "multiple comparisons tests.md",
        "Statistics" => "statistics.md",
        "Test statistics" => "test statistics.md",
        "p-value combinations" => "pCombination.md",
        "Tools" => "tools.md",
        "Create your own test" => "create your own test.md",
        "Chose a Test" => "chose a test.md"
    ]
)
