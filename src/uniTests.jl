#=
uniTests.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home

################################
# Univariate Permutation tests #
################################

=============================
  EXPORTED:
- _permTest! - ultimate function to perform all univariate tests 
- testStatistic - compute observed and permuted test statistic for univariate tests 

  UTILITIES:
- _randperm_uni!   - generate a random permutattion overwriting either the 𝐱 or 𝐲 vector
- _additional_kwargs_uni! - additional keyword arguments for _permTest! depending on the statistic
=#


# ----- #
# Generate a random permutation overwriting either vector 𝐱 or vector 𝐲, depending on stat and return them as a tuple.
# For both BivStatistic and IndSampStatistic this is a suffling of the elements of 𝐱 and ns is ignored 
# For RepMeasStatistic statistics this is a suffling of all elements of 𝐱 taken k at a time. ns is needed 
# For OneSampStatistic statistics, the sign of the elements in 𝐲 is flipped with probability 0.5. ns is ignored 
#   Note that instead for exact test only the 𝐱 vector is changed
_randperm_uni!(𝐱::UniData, 𝐲, stat::Union{BivStatistic, IndSampStatistic}, rng::MersenneTwister; ns::nsType = 0, ft = nothing) =
    return (shuffle!(rng, 𝐱), 𝐲)


function _randperm_uni!(𝐱::UniData, 𝐲, stat::RepMeasStatistic, rng::MersenneTwister; ns = @NamedTuple{n::Int, k::Int}, ft = nothing) 
    @simd for i=0:ns.n-1 
        @inbounds shuffle!(rng, view(𝐱, (i*ns.k)+1:(i*ns.k)+ns.k)) # shuffling in-place of 𝐱
    end
    return (𝐱, 𝐲)
end

function _randperm_uni!(𝐱::UniData, 𝐲, stat::OneSampStatistic, rng::MersenneTwister; ns::nsType = 0, ft::Tuple=(-1.0, 1.0)) 
    @simd for i ∈ eachindex(𝐲) 
        @inbounds 𝐲[i] *= rand(rng, ft) # flip the sign of the elements of 𝐲 with probability 0.5 (rand pick a number at random from `ft`)
    end
    return (𝐱, 𝐲)
end
# ----- #


# Add some keyword arguments to kwargs for some statistics.
# ----- #
function _additional_kwargs_uni!(kwargs, 𝐲, ns, stat, cpcd)
    stat isa AnovaF_RM && (kwargs = (kwargs..., ∑Y²kn=_∑Y²kn(𝐲, ns), ∑y²=_∑y²(𝐲), ∑S²k=_∑S²k(𝐲, ns)))
    stat isa StudentT_1S && (kwargs = (kwargs..., ∑y²=∑of²(𝐲)))
    cpcd === nothing || (kwargs = (kwargs..., cpcd=cpcd))
    return kwargs
end
# ----- #


# Compute observed and permuted statistics. Use an alias of `statistic` to allow a consistent
# sintax when creating custom univariate and multiple comparison tests
# ----- #
testStatistic(𝐱, 𝐲::UniData, stat::AllStatistics; kwargs...) = 
    statistic(𝐱, 𝐲, stat; kwargs...)
# ----- #



# ----- #
"""
```julia
function _permTest!(x, y, ns::nsType, stat::Stat, asStat::AsStat;
                    standardized::Bool=false, centered::Bool=false, 
                    nperm::Int = 20000, 
                    fstat::Function = abs,
                    compfunc::Function = >=,
                    switch2rand::Int = Int(1e8),
                    seed::Int = 1234,
                    verbose::Bool = true,
                    cpcd = nothing,
                    kwargs...) where {Stat<:Statistic, AsStat<:Statistic}
```

This function ultimately performs all **univariate permutation tests** implemented in *PermutationsTests.jl*, 
both *exact* and *approximate* (Monte Carlo). 

For running tests use the [univariate test functions](@ref "Univariate tests").
You need this function only for [creating your own tests](@ref "Create your own test").

Rewrite `x` and/or `y`, depending on the test performed.

For the `ns` argument see [ns](@ref).

`stat` can be a singleton of the [Statistic](@ref) type or a user-defined singleton of this type 
to be used as argument of a function implemented by the user
to compute both the observed and permuted test statistics, see [create your own test](@ref "Create your own test").

`asStat` is a singleton of the [Statistic](@ref) type. It is used to determine the permutation scheme 
and for this purpose it will be internally passed to [`genPerms`](@ref) and [`nrPerms`](@ref).

`asStat` determines also the input data format if you declare your own `stat` type. In this case
the function you write for computing the observed and permuted test statistic will take the `x` and 
`y` arguments as it follows: 

For `Stat` belonging to [group](@ref "Statistic groups")

 - `BivStatistic` : `x`, `y` are the two vectors of ``N`` elements each for which the bivariate statistic is to be tested. `ns` is ignored.
 - `IndSampStatistic` : we have ``K`` groups and ``N=N_1+...+N_K`` total observations; `y` holds all observations in a single vector such as `[y1;...;yK]` and `x` is the [`membership(::IndSampStatistic)`](@ref) vector. For example, for ``K=2``, ``N_1=2`` and ``N_2=3``, `x=[1, 1, 2, 2, 2]`.
 - `RepMeasStatistic` : we have ``K`` measures (*e.g.*, treatements) and ``N`` subjects; `y` holds the  ``K*N`` observations in a single vector such as `[y1;...;yN]`, where each vector ``y_i``, for ``i=1...N``, holds the observation at the ``K`` treatments and `x=collect(1:K*N)` (see [`membership(::RepMeasStatistic)`](@ref)).
 - `OneSampStatistic` : We have ``N`` observations (*e.g.*, subjects); `y` holds the ``N`` observations and `x=ones(Int, N)` (see [`membership(::OneSampStatistic)`](@ref)).

!!! note "Nota Bene"
    In all cases `x` is treated as the permutation vector that will be permuted before calling the [`testStatistic`](@ref)
    function.

For `length(x)>30` the approximate test is performed in all cases, 

otherwise,

if the number of systematic permutations exceeds `switch2rand` the approximate test is performed 
using `nperm` random permutations (default 20000), 

otherwise, 

the exact test is performed. 

`switch2rand` defaults to 1e8. 
To perform approximate tests in all cases, set `switch2rand`, for example, to 1. 

If `stat` is a `BivStatistic`, it optionally uses kwargs `standardized` or `centered` to compute them faster,
see for example [`correlationTest`](@ref).

`seed` is the initial seed for generating random permutations (not used for exact tests). 
To use a random seed, pass `seed=0`. For `seed` any natural number, the test will be reproducible.

`fstat` is a function applied to the test statistic. By default this is the julia `abs` function, 
which takes the absolute value, hence yieds a bi-directional test for a test statistic distributed symmetrically
around zero. For a right-directional test using such test statistics pass here `identity`. 
For a left-directional using such test statistics pass here [`flip`](@ref).

`compfunc` is the function to compare the observed statistics to the permuted statistics. 
The default function is `>=`. Don't change it unless you have studied the code of the function.

If `verbose` is true, print some information in the REPL while running the test. 
Set to false if you run benchmarks. The default is true.

For the `cpcd` and `kwargs...` argument, see [create you own test](@ref "Create your own test").

Return a [UniTest](@ref) structure.

*Examples*
```julia
using PermutationTests
x=randn(6)
y=randn(6) 
# bi-directional exact test of the correlation between x and y
t8 = _permTest!(μ0(x), μ0(y), length(x), Covariance(), Covariance(); centered=true) 
t8.p
t8.stat
#...

# make a left-directional test and standardize the data
t8_2 = _permTest!(μ0σ1(x), μ0σ1(y), length(x), Covariance(), Covariance(); 
        standardized=true, fstat=flip) 

# the same but force an approximate test
t8_3 = _permTest!(μ0σ1(x), μ0σ1(y), length(x), Covariance(), Covariance(); 
        standardized=true, fstat=flip, switch2rand=1) 

# the same using 5000 random permutations
t8_4 = _permTest!(μ0σ1(x), μ0σ1(y), length(x), Covariance(), Covariance(); 
        standardized=true, fstat=flip, switch2rand=1, nperm=5000) 
```

To check more examples, see the *uniTests_API.jl* unit located in the *src* github folder
and function `test_unitests()` in the *runtests.jl* unit located in the *test* github folder.
"""
function _permTest!(𝐱, 𝐲, ns::nsType, stat::Stat, asStat::AsStat;
                    standardized::Bool=false, centered::Bool=false, # means::Tuple=(), sds::Tuple=(), # optional kwa for correlation-like statistics
                    nperm::Int = 20000, 
                    fstat::Function = abs,
                    compfunc::Function = >=,
                    switch2rand::Int = Int(1e8),
                    seed::Int = 1234,
                    verbose::Bool = true,
                    cpcd = nothing,
                    kwargs...) where {Stat<:Statistic, AsStat<:Statistic}

    # check arguments and prepare test
    testtype, nperm, direction, design, mykwargs, rng = _prepare_permtest!(𝐱, 𝐲, ns, asStat, fstat, nperm, switch2rand, seed, standardized, centered)#, means, sds)      

    c, f, s, r! = compfunc, fstat, testStatistic, _randperm_uni! # aliases

    ####################### START TEST ##############################

    verbose && println("Performing a test using ", testtype == :exact ? "$nperm systematic permutations..." : "$nperm random permutations...")

    # observed statistic. eps is to avoid floating point arithmetic errors when comparing the permuted stats
    mykwargs = _additional_kwargs_uni!(mykwargs, 𝐲, ns, stat, cpcd) # more kwargs for some test statistics
    obsStat = f(testStatistic(𝐱, 𝐲, stat; mykwargs..., kwargs...)) - sqrt(eps()) 

    # get p-value
    #############
    #   The p-value is otained summing the number of time the function `fstat` applied to the permutated statistics satisfyies fucntionc `c` 
    #   (default: the permutation statistic is greater then or equal to the observed one) and dividing by the number of permutation.
    #   `fstat` is abs for two-tailes test, identity for right-tailed tests and flip for a left-tails test.
    #   For exact tests P is an itarator over all systematic permutations.
    #   Nota bene: for exact tests the permutation corresponding to the observed data is the first in P, 
    #       hence the sum across all permutations is considered. For approximate (monte carlo) tests 1 is added to the sum to occount 
    #       for the permutation corresponding to the observed data.
    if testtype == :exact
        P = genPerms(asStat, 𝐱, ns, direction, design) # generate systematic permutations

        # check to be removed later on
        #nperm ≠ length(P) && throw(ArgumentError("nperm is not equal to the length of P, $nperm, $(length(P))"))

        pvalue = sum(c(f(s(p, 𝐲, stat; mykwargs..., kwargs...)), obsStat) for p ∈ P) / nperm

    else    # Monte Carlo test
        # N.B. r! return 𝐱_, 𝐲_, one of which is modified depending on stat to generate a random data permutation 
        pvalue = (1+sum(c(f(s(r!(𝐱, 𝐲, asStat, rng; ns)..., stat; mykwargs..., kwargs...)), obsStat) for i ∈ 1:nperm-1)) / nperm  
    end
    #############

    return UniTest(pvalue, stat, obsStat, 1/nperm, nperm, testtype, direction, design)
    #return pvalue, stat, obsStat, 1/nperm, nperm, testtype, direction, design

end
# ----- #

# Note: for Onesampstat for exact tests x is changed, for approximate tests y is changed