#   This script is not part of the PermutatioTests.jl package.
#   It allows to build that package from the source code,
#   without actually installing the package.
#   It is useful for developing purposes using the Julia
#   `Revise` package 
#   You won't need this script for using the package.
#
#   MIT License
#   Copyright (c) 2024-2025, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the PermutationTests.jl from github or Julia registry, uninstall it.
#   2) Run this sript (With VS code, click anywhere here and hit ALT+Enter)
#
begin
  push!(LOAD_PATH, abspath(joinpath(dirname(@__FILE__), "..")))
  using PermutationTests
end
