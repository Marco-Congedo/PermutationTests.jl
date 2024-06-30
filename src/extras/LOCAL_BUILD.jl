#   This script is not part of the PermutatioTests.jl package.
#   It allows to build that package and its documentation locally from the source code,
#   without actually installing the package.
#   It is useful for developing purposes using the Julia
#   `Revise` package (that you need to have installed on your PC,
#   together with the `Documenter` package for building the documentation).
#   You won't need this script for using the package.
#
#   MIT License
#   Copyright (c) 2024, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the PermutationTests.jl from github or Julia registry, uninstall it.
#   2) Change the `juliaCodeDir` path here below to the path
#           where the PermutationTests.jl folder is located on your computer.
#   3) Run this sript (With VS code, click anywhere here and hit ALT+Enter)
#
begin
  juliaCodeDir = joinpath(homedir(),"Documents", "@ Documenti", "Code", "julia")
  scrDir       = joinpath(juliaCodeDir, "PermutationTests", "src")
  docsDir      = joinpath(juliaCodeDir, "PermutationTests", "docs")

  push!(LOAD_PATH, scrDir)
  using Documenter, Revise, PermutationTests

  cd(docsDir)
  clipboard("""makedocs(sitename="PermutationTests", modules=[PermutationTests], remotes = nothing)""")
  @info("\nRight-click and press ENTER on the REPL for building the documentation.");
end
