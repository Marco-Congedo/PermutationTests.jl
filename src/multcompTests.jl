#=
multcompTests.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home

##########################################
# Multiple Comparisons Permutation tests #
##########################################

Run `test_multcompTests()` in unit _test.jl to test this unit

=============================
  EXPORTED:

_permMcTest! - ultimate function to perform all multiple comparison tests 
testStatistic - compute the ith observed and permuted statistic for multiple comparison tests

  UTILITIES:

_preComputedData - pre-computed data for rep meas ANOVA and one-sample t-test for multiple comparison tests
_randperm_multComp! - generate a random permutation for multiple comparison tests
=============================
=#

######################################################################################################################################


# ----- #
# Pre-compute data that is invariant to permutations for repeated measure ANOVA
_preComputedData(𝐘::UniDataVec, ns, stat::AnovaF_RM) = [_∑Y²kn_∑y²_∑S²k(𝐲, ns) for 𝐲 ∈ 𝐘]
_preComputedData(𝐘::UniDataVec, ns, stat::StudentT_1S) = [∑of²(𝐲) for 𝐲 ∈ 𝐘] 
_preComputedData(𝐘::UniDataVec, ns, stat::ParFreeStatistics) = nothing 
_preComputedData(𝐘, ns, stat) = nothing # allow calling with custom data types and or custum statistic
# ----- #
    


# ----- #
# Generate a random permutation overwriting either vector 𝐱 or all vectors in 𝐘, depending on stat. Does not return anything.
# For both BivStatistic and IndSampStatistic this is a suffling of the elements of 𝐱 and ns is ignored 
# For RepMeasStatistic statistics this is a suffling of all elements of 𝐱 taken k at a time. ns is needed 
# For OneSampStatistic statistics, the sign of the elements in the 𝐲 vectors of 𝐘 is flipped with probability 0.5. ns is ignored 
#    Only for this latter case this procedure is algorithmically different as compared to the univariate tests, but here never return anything.
# Note that instead for exact test only the 𝐱 vector is changed
_randperm_multComp!(𝐱::UniData, 𝐘, stat::Union{BivStatistic, IndSampStatistic}, rng::MersenneTwister; ns::nsType = 0, ft = nothing) = 
    shuffle!(rng, 𝐱)

function _randperm_multComp!(𝐱::UniData, 𝐘, stat::RepMeasStatistic, rng::MersenneTwister; ns = @NamedTuple{n::Int, k::Int}, ft = nothing) 
    @simd for i=0:ns.n-1 
        @inbounds shuffle!(rng, view(𝐱, (i*ns.k)+1:(i*ns.k)+ns.k)) # shuffling in-place of 𝐱
    end
end

function _randperm_multComp!(𝐱::UniData, 𝐘, stat::OneSampStatistic, rng::MersenneTwister; ns::nsType = 0, ft::Tuple=(-1.0, 1.0)) 
    signs = [rand(rng, ft) for i = 1:length(𝐘[1])] # vector of random signa with probability 0.5 (rand pick a number at random from `ft`)
    @simd for i ∈ eachindex(𝐘) 
        @inbounds 𝐘[i] .*= signs # the vector of signs is the same for all variables
    end
end
# ----- #



# ----- #
"""
```julia

# METHOD 1
function testStatistic(x, y, stat::mystat, fstat::Function; 
                        cpcd=nothing, kwargs...)

# METHOD 2
function testStatistic(x, Y, i::Int, stat::mystat, fstat::Function; 
                        cpcd=nothing, kwargs...)
                        
        where mystat<:Statistic 
```

Compute the observed and permuted test statistic for univariate tests (Method 1) or the ``i^{th}`` 
observed and permuted test statistic for the ``i^{th}`` hypothesis, with ``i=1...M``, 
for multiple comparisons permutation tests (Method 2).

If you [create your own test](@ref "Create your own test") you will write new methods
for these functions taking as `stat` a test statistic of type [Statistic](@ref) you have declared.

If not, you never need these functions.

`Y` is a vector of elements (typically, vectors themelves) and the test-statistic is to be computed on `Y[i]`,
using the permutation vector `x`.

For the `fstat` and `cpcd` argument, see how to [create your own test](@ref "Create your own test").

"""
testStatistic(𝐱, 𝐘::UniDataVec, i::Int, stat::AnovaF_RM; cpcd=nothing, kwargs...) = 
    statistic(𝐱, 𝐘[i], stat; ∑Y²kn=cpcd[i][1], ∑y²=cpcd[i][2], ∑S²k=cpcd[i][3], kwargs...)

testStatistic(𝐱, 𝐘::UniDataVec, i::Int, stat::StudentT_1S; cpcd=nothing, kwargs...) = 
    statistic(𝐱, 𝐘[i], stat; ∑y²=cpcd[i], kwargs...) 

# all other tests implemented in PermutationsTests.jl    
testStatistic(𝐱, 𝐘::UniDataVec, i::Int, stat::ParFreeStatistics; cpcd=nothing, kwargs...) =
    statistic(𝐱, 𝐘[i], stat; kwargs...) 
# ----- #


# ----- #
"""
```julia
function _permMcTest!(x, Y, ns::nsType, stat::Stat, asStat::AsStat;
            standardized::Bool=false, centered::Bool=false, 
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            nperm::Int = 20000, 
            fstat::Function = abs,
            compfunc::Function = >=,
            switch2rand::Int = Int(1e8),
            seed::Int = 1234,
            threaded::Bool = Threads.nthreads()>=4,
            verbose::Bool = true,
            cpcd = nothing) where {Stat<:Statistic, AsStat<:Statistic}
```

This function ultimately performs all **multiple comparisons permutation tests** implemented in *PermutationsTests.jl*, 
both *exact* and *approximate* (Monte Carlo). 

For running tests use the [multiple comparisons test functions](@ref "Multiple comparisons tests").
You need this function only for [creating your own tests](@ref "Create your own test").

The step-down version of the test is performed if `stepdown` is true (default).
In this case the `fwe` (family-wise error) rate is used for rejection at each step (default=0.05).

Rewrite `x` and/or `Y`, depending on the test performed.

For the `ns` argument see [ns](@ref).

`stat` can be a singleton of the [Statistic](@ref) type or a user-defined singleton of this type 
to be used as argument of two functions implemented by the user
to compute the observed and permuted test statistics, see [create your own test](@ref "Create your own test").

!!! warning
    In contrast to univariate tests, equivalent statistics are not possible for multiple comparison tests, 
    with the exception of `CrossProd()` if the data has been standardized and `Covariance()` if the data 
    has been centered and those only for correlation-like tests.
    
    The test statistics that must be used as `Stat` for the other kinds of test if a singleton of the 
    [Statistic](@ref) type is used are `AnovaF_IS()`, `StudentT_IS()`, `AnovaF_RM()`, and `StudentT_1S()`.

`asStat` is a singleton of the [Statistic](@ref) type. It is used to determine the permutation scheme 
and for this purpose it will be passed to [`genPerms`](@ref) and [`nrPerms`](@ref).

`asStat` determines also the input data format if you declare your own `stat` type. In this case
the two functions you write for computing the observed and permuted statistics will take the `x` and 
`Y` arguments as it follows. 

For `Stat` belonging to [group](@ref "Statistic groups")

 - `BivStatistic` : `x` is a fixed vector with ``N`` elements and `Y` an an ``M``-vector of ``N``-vectors. The ``M`` bivariate statistics between `x` and the ``y_m`` vectors of `Y` are tested simultaneously. `ns` is ignored.
 - `IndSampStatistic` : we have ``K`` groups and ``N=N_1+...+N_K`` total observations; `Y` is an ``M``-vector, each one holding all observations for the ``m^{th}`` hypothesis in a single vector. The ``m^{th}`` vector ``y_m`` concatenates the observations for all groups such as [y[m][1];...;y[m][K]]. `x` is the [`membership(::IndSampStatistic)`](@ref) vector, common to all hypotheses. For example, for ``K=2``, ``N_1=2`` and ``N_2=3``, `x=[1, 1, 2, 2, 2] `. 
 - `RepMeasStatistic` : we have ``K`` measures (*e.g.*, treatements) and ``N`` subjects; `Y` is an M-vector, each one holding the ``K*N`` observations for the ``m^{th}`` hypothesis in a single vector. The ``m^{th}`` vector ``y_m`` is such as `[Y[m][1];...;Y[m][N]]`, where each vector Y[1][n], for ``n=1…N``, holds the observations for the ``K`` treatments and `x=collect(1:K*N)` (see [`membership(::RepMeasStatistic)`](@ref)). 
 - `OneSampStatistic` : We have ``N`` observations (*e.g.*, subjects); `Y` is an ``M``-vector, each one holding the ``N`` observations and `x=ones(Int, N)` (see [`membership(::OneSampStatistic)`](@ref)).


!!! note "Nota Bene"
    In all cases `x` is treated as the permutation vector that will be permuted before calling the [`testStatistic`](@ref)
    function for each of the elements in `Y`.

Optional keyword arguments `switch2rand`, `nperm`, `standardized`, `centered`, `seed`, `fstat`, `compfunc` and `verbose`
have the same meaning as in the [`_permTest!`](@ref) function.

If `threaded` is true (default) the function is multi-threaded if the product of the number of hypotheses, 
    observations, and permutations exceed 500 millions.
If you have unexpected problems with the function, try setting `threaded` to false.

For the `cpcd` and `kargs...` arguments, see [create you own test](@ref "Create your own test").

Return a [MultcompTest](@ref) structure.

The number of executed steps ``S`` can be retrived as the length of the `.rejections` field 
of the returned structure. 

*Examples*
```julia 
using PermutationTests
N, M = 8, 100 # 100 hypotheses, N=8
x=randn(N)
Y=[randn(N) for m=1:M]
# bi-directional exact test of the correlation between x and 
# all the M vector in Y.
T12 = _permMcTest!(x, Y, N, PearsonR(), PearsonR())
T12.p
T12.stat
#...

# bi-directional exact test. Faster test by data standardization 
T12_2 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y in Y], N, CrossProd(), PearsonR(); 
                    standardized=true) 

# left-directional exact test.
T12_3 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y in Y], N, CrossProd(), PearsonR(); 
                    standardized=true, fstat=flip)

# as above, but force an approximate test
T12_4 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y in Y], N, CrossProd(), PearsonR(); 
                    standardized=true, fstat=flip, switch2rand=1)

# the same using 5000 random permutations
T12_4 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y in Y], N, CrossProd(), PearsonR(); 
                    standardized=true, fstat=flip, switch2rand=1, nperm=5000)
```

To check more examples, see the *multcompTests_API.jl* unit located in the *src* github folder
and function `test_multicompTests()` in the *runtests.jl* unit located in the *test* github folder.
"""
function _permMcTest!(𝐱, 𝐘, ns::nsType, stat::Stat, asStat::AsStat;
            standardized::Bool=false, centered::Bool=false, # means::Tuple=(), sds::Tuple=(), # optional kwa for correlation-like statistics
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            nperm::Int = 20000, 
            fstat::Function = abs,
            compfunc::Function = >=,
            switch2rand::Int = Int(1e8),
            seed::Int = 1234,
            threaded::Bool = Threads.nthreads()>=4,
            verbose::Bool = true,
            cpcd = nothing,
            kwargs...) where {Stat<:Statistic, AsStat<:Statistic}
    
    # check that all vectors in 𝐘 have the same length. This is special for multivariate tests
    length(unique(length(𝐲) for 𝐲 ∈ 𝐘)) ≠ 1 && throw(ArgumentError(📌*" Function _permMcTest!: The elements (typically, vectors) in 𝐘 do not have all the same length"))

    # check arguments and prepare test
    testtype, nperm, direction, design, mykwargs, rng = _prepare_permtest!(𝐱, 𝐘[1], ns, asStat, fstat, nperm, switch2rand, seed, standardized, false)#, (), ())      

    # pre-computed data in case of some statistics. This is special for multivariate tests. Set to nothing for custom statistics
    pcd = _preComputedData(𝐘, ns, stat) 

#    println("pcd: ", pcd)

    ####################### START TEST ##############################

    # observed statistics for all variables in 𝐘
    # eps is to avoid floating point arithmetic errors when comparing the permuted stats
    obsStats = [fstat(testStatistic(𝐱, 𝐘, i, stat; cpcd=pcd, mykwargs..., kwargs...)) - sqrt(eps()) for i ∈ eachindex(𝐘)]

    verbose && println("Performing a ", threaded ? "(multi-threaded) " : "", "test using ", testtype == :exact ? "$nperm systematic permutations..." : "$nperm random permutations...")
    
    # initialization and reserve required memory upfront to boost performance 
    accepted = trues(length(𝐘)) # BitArray, mask of accepted H0. Be careful, for multi-threading it is not thread safe. Use booleans instead
    rejections = Vector{Vector{Int64}}() # keep track of the rejected hypotheses at each step (a vector of indeces at each step)
    pvalues = Vector{Float64}(undef, length(𝐘)) # corrected p-values
    nullDistr = Vector{Float64}(undef, nperm) # Null distribution (only the one at last step will be returned)
    maxrejp = 0. # max rejected p-value at current step.
    step = 0
    nSig = 0
    hasRejected = true # set to false to exit the while loop and stop the step-down procedure


    # perform max-statistic or step-down max-statistic test, depending on whether `stepdown` is true.
    while hasRejected  

        # get p-values
        ############################################################
        #   For each variable the p-value is otained as the number of elements in the null distribution of max statistics satisfting fenction 
        #   `compfunc` when compared to the the observed statistic for that variable (>= by default, don't touch it) divided by the 
        #   number of permutations. `fstat` is abs for two-tailes test, identity for right-tailed tests and flip for a left-tails test.
        #   For exact tests P is a lazy itarator over all systematic permutations.
        #   Nota bene: for exact tests the permutation corresponding to the observed data is the first in P, thus the observed max statistic
        #   is naturally included in the null distribution. For the approximate (monte carlo) tests, the observed max statistic is added
        #   to the null distribution.

        tx = threaded && count(accepted)*length(𝐘[1])*nperm>5e8 # decide if using multithreading. count(accepted) considers the remaining Hyp.
        maxf = tx ? Folds.maximum : Base.maximum # the only multi-threading function is the maximum of the p-value

        if testtype == :exact
            P = genPerms(asStat, 𝐱, ns, direction, design) # generate systematic permutations as a lazy iterator
            if asStat isa StudentT_1S # bug fix: for StudentT_1S only, the iterator yields tuples and can be enumerated 
                for (j, p) in enumerate(P) #  and this iterator works only expliciting the for loop
                    nullDistr[j] = maxf(fstat(testStatistic(p, 𝐘, i, stat; cpcd=pcd, mykwargs..., kwargs...)) for i ∈ eachindex(𝐘, accepted) if accepted[i])
                end
            else # on the other hand the other iterators cannot be enumerated
                nullDistr = [maxf(fstat(testStatistic(p, 𝐘, i, stat; cpcd=pcd, mykwargs..., kwargs...)) for i ∈ eachindex(𝐘, accepted) if accepted[i]) for p ∈ P]
            end
        else    # Monte Carlo test
            for j ∈ 1:nperm-1
                _randperm_multComp!(𝐱, 𝐘, asStat, rng; ns) # one of 𝐱, 𝐘 is modified depending on stat to generate a random data permutation
                nullDistr[j] = maxf(fstat(testStatistic(𝐱, 𝐘, i, stat; cpcd=pcd, mykwargs..., kwargs...)) for i ∈ eachindex(𝐘, accepted) if accepted[i])
            end
            nullDistr[end] = _condMax(obsStats, accepted)
        end

        # find significance (p-values) and indices of rejected hypotheses, flag rejected hypotheses, update highest rejected p-value,
        # limit significance to the highest p-value that has been rejected at previous step
        sort!(nullDistr, rev=true) # descending order 
        highestRejected = 0.
        rejpos=Int[]   
        @simd for i ∈ eachindex(pvalues)
            @inbounds if accepted[i] # at first pass all hypothesis are accepted
                pos = findlast(x->compfunc(x, obsStats[i]), nullDistr) 
                pvalues[i] = pos===nothing ? 1/nperm : pos/nperm # find significance 
                if pvalues[i] <= fwe # if there is a rejection
                    push!(rejpos, i) # keep track of rejections for giving it as output
                    accepted[i] = false # flag this hypotheses as rejected  
                    pvalues[i] > highestRejected && (highestRejected = pvalues[i]) # find the highest rejected p-value
                    pvalues[i] = max(pvalues[i], maxrejp) # limit significance to the highest p-value that has been rejected at previous step 
                    nSig += 1 # number of rejected hypotheses
                end
            end
        end
        maxrejp = highestRejected # Limit significance hereafter to the highest rejected pvalue at this step        

        ############################################################

        if isempty(rejpos) # no rejeced hypothesis 
            verbose && print("Step ", step+1, ": no rejection; step-down process stops.\n")
        else
            push!(rejections, rejpos) # keep track of the rejected hypotheses at each step
            length(rejpos)<30 ?     verbose && println("Step ", step+1, ": rejected hypotheses ", rejpos) : 
                                    verbose && println("Step ", step+1, ": ", length(rejpos), " hypotheses have been rejected")
        end
    
        # exit while stepdown=false or if there were no rejeced hypothesis
        hasRejected = !stepdown || isempty(rejpos) || nSig==length(𝐘) ? false : true                 

        step += 1
    end # while
    ############################################################

    return MultcompTest(pvalues, stat, obsStats, 1/nperm, nperm, testtype, direction, design, nullDistr, rejections, stepdown, fwe)
end
# ----- #

# ----------------------------------------------------- #