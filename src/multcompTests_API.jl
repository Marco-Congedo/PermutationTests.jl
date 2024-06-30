#=
multcompTests_API.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home

#################################################
# API for Multiple Comparison Permutation tests #
#################################################

Run `test_mult_compTests_API()` in unit _test.jl to test this unit

========================================
  EXPORTED:

correlationMcTest!, rMcTest!,
correlationMcTest,  rMcTest,
trendMcTest!, trendMcTest,
pointBiSerialMcTest,
anovaMcTestIS, fMcTestIS,
studentMcTestIS, tMcTestIS, 
chiSquaredMcTest, Χ²McTest, 
fisherExactMcTest,
anovaMcTestRM, fMcTestRM,
cochranqMcTest, qMcTest,
mcNemarMcTest,
studentMcTestRM, McTestRM,
studentMcTestRM!, McTestRM!,
studentMcTest1S!, tMcTest1S!,
studentMcTest1S, tMcTest1S,
signMcTest!, signMcTest
========================================
=#

# GENERAL KEYWORD ARGUMENTS IN COMMON WITH UNIVARIATE TESTS
# ---------------------------------------------------------
# direction: test direction, either Right(), Left() or Both(). Default: Both()
# switch2rand (Int): the # of permutations starting from which the approximate test will be carried out. Default: 1e8
# nperm (Int): the # of permutations to be used if the test will be approximate (default: 20_000) 
# seed (Int): random # generator seed. Set to 0 to take a random seed. Any natural number allow a reproducible test.(Default: 1234) 
# verbose (Bool) print some information in the REPL while running the test. Set to false for running @benchmark (Default: true).

# GENERAL KEYWORD ARGUMENTS THAT ARE SPECIAL TO MULTIPLE COMPARISON TESTS (NB: no `equivalent` arg is provided for multiple comparisons)
# -----------------------------------------------------------------------
# stepdown: if true the step-down procedure is applied
# fwe: rejection threshold for step-down procedures
# threaded: if true the function will be run in multithreading for large problems
##################################################################################

### Correlation test

# Example scenario of a multiple comparison correlation test:
#   A neuroscientist carries out an fMRI experiment on N=16 subjects, obtaining a metabolic measure of neuronal workload 
#   at each of M=24000 brain voxels.
#   The experiment involes a cognitive task allowing a measure of latency to complete the task for each subject.
#   The scientist wish to know if there exist one ore more brain region which metabolism correlate or anticorrelate with the latency.
#   Let 𝐱 be the vector of N latencies for each subject and let 𝐘=[𝐲1, ..., 𝐲M] be the vector of M vectors, holding each one 
#   the N methabolic measures for each subject, in the same order as in 𝐱. The multiple comparion test is obtained as
#   rMc = rMcTest!(𝐱, 𝐘) # bi-directional by default
#   If only a negative correlation was expected, the test would have been
#   rMc = rMcTest!(𝐱, 𝐘; direction=(Left))
"""
```julia
function correlationMcTest(x::UniData, Y::UniDataVec;
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true, 
            #
            standardized::Bool = false,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection
```

Multiple comparisons [Pearson product-moment correlation test](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Testing_using_Student's_t-distribution) 
by data permutation. 

Run ``M`` correlation tests simultaneously. 
The input data are a fixed vector `x` and ``M`` vectors given as `Y`, a vector of ``M`` vectors ``y_1,...,y_M``. 

`x` and all vectors in `Y` must have equal length. 

The ``M`` null hypotheses have form 

``H_0(m):r_{(x, y_m)}=0, \\quad m=1...M``, 

where ``r_{x,y_m}`` is the Pearson correlation coefficient between vector `x` and the ``m^{th}`` vector of `Y`.

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded`, see [here](@ref "Common kwargs for multiple comparisons tests").

If `standardized` is true, both `x` and all vectors of `Y` are assumed standardized (zero mean and unit standard deviation).
With standardized input data the test can be executed faster as in this case the cross-product is actually the correlation. 
If `standardized` is false, the data will be standardized to execute a faster test.

*Directional tests, permutation scheme and number of permutations for exact tests:* as per 
*univariate version* [`correlationTest`](@ref)

*Aliases:* `rMcTest`, [`trendMcTest`](@ref) 

Return a [MultcompTest](@ref) structure.


*Examples*
```julia
using PermutationTests
N=10; # number of observations
M=100; # number of tests
x=randn(N);
Y=[randn(N) for i=1:M];
t=rMcTest(x, Y) # bi-directional test
```
```julia
tR=rMcTest(x, Y; direction=Right()) # right-directional test
tL=McTest(x, Y; direction=Left()) # left-directional test
tMC=rMcTest(x, Y; switch2rand=1) # Force a monte carlo test
# Force a monte carlo test and performs 50K permutations
t5K=rMcTest(x, Y; switch2rand=1, nperm= 50_000) 
tnoSD=rMcTest(x, Y; stepdown=false) # don't do stepdown
tnoMT=rMcTest(x, Y; threaded=false) # don't run it multithreaded
t001=rMcTest(x, Y; fwe=0.01) # stepdown rejects at 0.01 level instead of 0.05
```

**Similar tests**

See [correlationTest](@ref)

"""
function correlationMcTest(𝐱::UniData, 𝐘::UniDataVec;
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true, 
            #
            standardized::Bool = false,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘)==1 && (return correlationTest(𝐱, 𝐘[1]; direction, switch2rand, nperm, seed, standardized, verbose))            
    length(𝐱) == length(𝐘[1]) || throw(ArgumentError(📌*"Function correlationMcTest: the first argument (𝐱) and the vectors in the second (𝐘) must have equal length. Check the documentation"))
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function correlationMcTest: the vectors in the second argument (𝐘) must have all equal length. Check the documentation"))

    return standardized ?   _permMcTest!(copy(𝐱), 𝐘, length(𝐱), PearsonR(), PearsonR(); 
                        standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose) :
                            _permMcTest!(μ0σ1(𝐱), [μ0σ1(𝐲) for 𝐲 in 𝐘], length(𝐱), PearsonR(), PearsonR(); 
                        standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
end        

# alias. The trendTest is obtained passed the trend to be tested as vector `𝐱`. For example [1, 2,...] for an ascending linear trend                        
rMcTest = correlationMcTest


# same as correlationMcTest! but in no case vector 𝐱 will be overwritten.
"""
```julia
function correlationMcTest!(<same args and kwargs as `correlationMcTest`>)
```

As [`correlationMcTest`](@ref), but `x` is overwritten if `standardized` is true. 

*Aliases:* `rMcTest!`, [`trendMcTest!`](@ref)

*Univariate version:*  [`correlationTest!`](@ref)

"""
function correlationMcTest!(𝐱::UniData, 𝐘::UniDataVec;
    direction::TestDir = Both(),
    switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)), 
    nperm::Int = 20_000, 
    seed::Int = 1234, 
    verbose::Bool = true, 
    #
    standardized::Bool = false,
    #
    stepdown::Bool = true,
    fwe::Float64 = 0.05,
    threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

length(𝐘)==1 && (return correlationTest!(𝐱, 𝐘[1]; direction, switch2rand, nperm, seed, standardized, verbose))   
length(𝐱) == length(𝐘[1]) || throw(ArgumentError(📌*"Function correlationMcTest!: the first argument (𝐱) and the vectors in the second (𝐘) must have equal length. Check the documentation"))
length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function correlationMcTest!: the vectors in the second argument (𝐘) must have all equal length. Check the documentation"))

return standardized ?   _permMcTest!(𝐱, 𝐘, length(𝐱), PearsonR(), PearsonR(); 
                        standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose) :
                    _permMcTest!(μ0σ1(𝐱), [μ0σ1(𝐲) for 𝐲 in 𝐘], length(𝐱), PearsonR(), PearsonR(); 
                        standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
end

# alias
rMcTest! = correlationMcTest! 
 

"""
```julia
function trendMcTest(<same args and kwargs as correlationMcTest>)
```
"""
trendMcTest = correlationMcTest

"""
```julia
function trendMcTest!(<same args and kwargs as correlationMcTest!>)
```

Actually aliases of [`correlationMcTest`](@ref) and [`correlationMcTest!`](@ref), respectively. 

`x` is any specified trend (linear, polynomial, exponential, logarithmic, trigonometric, ...)
and `Y` holds the observed data. A multiple comparison Pearson product-moment correlation 
test between `x` and all ``M`` vectors in `Y` is then carried out. 

Directional tests, permutation scheme and number of permutations for exact tests: as per 
*univariate versions* [`trendTest`](@ref) or [`trendTest!`](@ref)

Both methods return a [MultcompTest](@ref) structure.

*Examples*

```julia
using PermutationTests
# We are goint to test an upward linear trend
N=10
M=100
x=Float64.(collect(Base.OneTo(N))) # [1, 2,..., N]
# out of the M vectors created here below, only one has a significant correlation
Y=[[1., 2., 4., 3., 5., 6., 7., 8., 10., 9.], ([randn(N) for m=1:M-1]...)];
# Since we expect an upward linear trend, the correlation is expected to be positive,
# hence we use a right-directional test to increase the power of the test.
t = trendMcTest(x, Y; direction=Right()) 
```
""" 
trendMcTest! = correlationMcTest!


### ANOVA for independent samples
"""
```julia
# METHOD (1)
function anovaMcTestIS(Y::UniDataVec, ns::IntVec;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection
```
"""
function anovaMcTestIS(𝐘::UniDataVec, ns::IntVec;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘)==1 && (return anovaMcTestIS(𝐘[1], ns; direction, switch2rand, nperm, seed, verbose))                            
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function anovaMcTestIS: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))
    !(direction isa Both) && throw(ArgumentError(📌*"Function anovaMcTestIS: The ANOVA test can only be bi-directional. Correct the `direction` keyword argument. Check the documentation")) 
    length(ns) < 2 && throw(ArgumentError(📌*"Function anovaMcTestIS: the second argument (ns) must be a vector of two or more integers. Check the documentation"))
    sum(ns) ≠ length(𝐘[1]) && throw(ArgumentError(📌*"Function anovaMcTestIS: the sum of the integers in second argument (ns) must be equal to the length of the vectprs in first argument. Check the documentation"))

    𝐱 = membership(AnovaF_IS(), ns)
    return length(ns)>2 ?   _permMcTest!(𝐱, 𝐘, ns, AnovaF_IS(), AnovaF_IS(); stepdown, fwe, nperm, switch2rand, seed, threaded, verbose) :
                            _permMcTest!(𝐱, 𝐘, ns, StudentT_IS(), StudentT_IS(); stepdown, fwe, nperm, switch2rand, seed, threaded, verbose)
end

# method 2: Input data is given as a multiple comparison vector of K vectors of observations for group 1...group K. Each vector has arbitrary length. 
"""
```julia
# METHOD (2)
function anovaMcTestIS(𝐘vec::UniDataVec²; <same kwargs>)
```

**METHOD (1)**

Multiple comparisons [1-way analysis of variance (ANOVA) for independent samples](https://en.wikipedia.org/wiki/One-way_analysis_of_variance) 
by data permutation. 

Run ``M`` ANOVAs simultaneously. The Input data is given as a vector `Y` holding
``M`` vectors ``𝐲1,...,𝐲M`` concatenating all observations, that is, holding each ``N=N_1+...+N_K`` observations 
for ``K>2`` independent samples (groups). The observations are ordered with group 1 first, then group 2,..., finally group K.
Note that the group numerosity ``N_1,...,N_K`` must be the same for all ``M`` hypotheses. 
The only check performed is that the first vector in `Y` contains `sum(ns)` elements.

The ``M`` null hypotheses have form

``H_0(m): μ_{m1}= \\ldots =μ_{mK}, \\quad m=1...M``,

where ``μ_{mk}`` is the mean of the ``k^{th}`` group for the ``m^{th}`` hypothesis.

`ns` is a vector of integers holding the group numerosity ``N_1,...,N_K`` (see examples below).

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests").

Directional tests, permutation scheme and number of permutations for exact tests: as per
*univariate version* [`anovaTestIS`](@ref)

*Alias:* `fMcTestIS`

Return a [MultcompTest](@ref) structure.


**METHOD (2)**

As method (1), but input data `Yvec` is a vector holding ``M`` vectors of K vectors of 
observations for group 1,..., group K. 

*Examples*

```julia
# method (1)
using PermutationTests
ns=[3, 4, 5] # number of observations in group 1, 2 and 3
M=10 # number of hypotheses
Yvec = [[randn(n) for n in ns] for m=1:M]; # some random Gaussian data for example 
t = fMcTestIS([vcat(y...) for y in Yvec], ns) # ANOVA tests are always bi-directional

# Force an approximate test with 5000 random permutations
tapprox = fMcTestIS([vcat(y...) for y in Yvec], ns; switch2rand=1, nperm=5000) 

# in method (2) only the way the input data is formatted is different 
t2 = fMcTestIS(Yvec)
# of course, method (1) and (2) give the same p-values
println(sum(abs.(t.p-t2.p))≈0. ? "OK" : "error") 
```

**Similar tests**

See [`anovaTestIS`](@ref)
"""
function anovaMcTestIS(𝐘vec::UniDataVec²;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘vec), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘vec)==1 && (return anovaMcTestIS(𝐘vec[1]; direction, switch2rand, nperm, seed, verbose))                                            
    !(direction isa Both) && throw(ArgumentError(📌*"Function anovaMcTestIS: The ANOVA test can only be bi-directional. Correct the `direction` keyword argument. Check the documentation")) 
    K = length(𝐘vec[1])
    K < 2 && throw(ArgumentError(📌*"Function anovaMcTestIS: the vectors in first argument (𝐘vec) must hold two or more vectors. Check the documentation"))
    ns=[length(𝐲) for 𝐲 in 𝐘vec[1]]
    unique(length(y) for y in 𝐘vec)[1]==length(ns) || throw(ArgumentError(📌*"Function anovaMcTestIS: All vectors in first arguments (𝐘vec) must contain the same number of vectors. Check the documentation")) 
    unique(length.(y) for y in 𝐘vec)[1]==ns || throw(ArgumentError(📌*"Function anovaMcTestIS: All vectors in first arguments (𝐘vec) must contain vectors of equal length in the same order. Check the documentation"))

    𝐱 = membership(AnovaF_IS(), ns)
    return K>2 ?    _permMcTest!(𝐱, [vcat(𝐲vec...) for 𝐲vec ∈ 𝐘vec], ns, AnovaF_IS(), AnovaF_IS(); stepdown, fwe, nperm, switch2rand, seed, threaded, verbose) :
                    _permMcTest!(𝐱, [vcat(𝐲vec...) for 𝐲vec ∈ 𝐘vec], ns, StudentT_IS(), StudentT_IS(); stepdown, fwe, nperm, switch2rand, seed, threaded, verbose)
end

# alias
fMcTestIS = anovaMcTestIS

### t-test for independent samples

"""
```julia
# METHOD (1)
function studentMcTestIS(Y::UniDataVec, ns::IntVec;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4,
                asPearson::Bool = true) where TestDir <: TestDirection
```
"""
function studentMcTestIS(𝐘::UniDataVec, ns::IntVec;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4,
                asPearson::Bool = true) where TestDir <: TestDirection

    length(𝐘)==1 && (return studentTestIS(𝐘[1], ns; direction, equivalent=true, switch2rand, nperm, seed, verbose, asPearson))
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function studentMcTestIS: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))
    length(ns) ≠ 2 && throw(ArgumentError(📌*"Function studentMcTestIS: the second argument (ns) must be a vector of two integer"))
    sum(ns) ≠ length(𝐘[1]) && throw(ArgumentError(📌*"Function studentMcTestIS: the sum of the two integers in second argument (ns) must be equal to the length of the vectors in first argument (𝐘)."))

    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱 = membership(StudentT_IS(), ns; rev=reverse)
        return _permMcTest!(μ0σ1(𝐱), [μ0σ1(𝐲) for 𝐲 in 𝐘], length(𝐱), CrossProd(), PearsonR();
                standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    else
        𝐱 = membership(StudentT_IS(), ns)
        return _permMcTest!(𝐱, 𝐘, ns, StudentT_IS(), StudentT_IS(); 
                stepdown, fwe, fstat=_fstat(direction), nperm, switch2rand, seed, threaded, verbose)
    end
end

# method 2: input data is given as a multiple comparison-vector of two vectors of arbitrary length, holding each the data of group 1 and group 2, 
# respectively.
"""
```julia
# METHOD (2)
function studentMcTestIS(Yvec::UniDataVec²; <same kwargs>)
```
"""
function studentMcTestIS(𝐘vec::UniDataVec²;
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘vec), Int(1e4)), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4,
            asPearson::Bool = true) where TestDir <: TestDirection

    length(𝐘vec)==1 && (return studentTestIS(𝐘vec[1]; direction, equivalent=true, switch2rand, nperm, seed, verbose, asPearson))
    K = length(𝐘vec[1])
    K ≠ 2 && throw(ArgumentError(📌*"Function studentMcTestIS: the vectors in first argument must hold two vectors"))
    ns=[length(𝐲) for 𝐲 in 𝐘vec[1]]
    unique(length(y) for y in 𝐘vec)[1]==length(ns) || throw(ArgumentError(📌*"Function studentMcTestIS: All vectors in first arguments (𝐘vec) must contain the same number of vectors. Check the documentation")) 
    unique(length.(y) for y in 𝐘vec)[1]==ns || throw(ArgumentError(📌*"Function studentMcTestIS: All vectors in first arguments (𝐘vec) must contain vectors of equal length in the same order. Check the documentation"))

    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱 = membership(StudentT_IS(), ns; rev=reverse)
        return _permMcTest!(μ0σ1(𝐱), [μ0σ1(vcat(𝐲vec...)) for 𝐲vec ∈ 𝐘vec], length(𝐱), CrossProd(), PearsonR(); 
                standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    else
        𝐱 = membership(StudentT_IS(), ns)
        return _permMcTest!(𝐱, [vcat(𝐲vec...) for 𝐲vec ∈ 𝐘vec], ns, StudentT_IS(), StudentT_IS(); 
                stepdown, fwe, fstat=_fstat(direction), nperm, switch2rand, seed, threaded, verbose)            
    end
end

# method 3: Input data is given as two multiple comparison-vectors of vectors of observations, 𝐗 for group 1 and 𝐘 for group 2. 
# The M vectors in 𝐗 have all the same length, so do the M vectors in 𝐘, but the vectors in 𝐗 and in 𝐘 don't have to be of the same length. 
"""
```julia
# METHOD (3)
function studentMcTestIS(X::UniDataVec, Y::UniDataVec; <same kwargs>)
```

**METHOD (1)**

Multiple comparisons [Student's t-test for independent samples](https://en.wikipedia.org/wiki/Student's_t-test#Independent_(unpaired)_samples) 
by data permutation. 

Run ``M`` t-tests simultaneously. Given ``M`` hypotheses with ``N=N_1+N_2`` observations for two groups each, 
the ``M`` null hypotheses have form

``H_0(m): μ_{m1}=μ_{m2}, \\quad m=1...M``,

where ``μ_{m1}`` and ``μ_{m2}`` are the mean of group 1 and group 2, respectively, for the ``m^{th}``hypothesis. 

For a bi-directional test, this t-test is equivalent to the 1-way ANOVA for two independent samples.
However, in contrast to the ANOVA, it can be directional. 

`ns` is a vector of integers holding the group numerosity ``N_1, N_2`` (see examples below).

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests"). 

If `asPearson` is true(default), the test is run as an equivalent version of a Pearson correlation test.
This is in general advantageous for multiple comparison tests, especially if approximate
(see the [benchmarks](@ref "Benchmarks")).
If you seek best performance for exact tests, benchmark the speed of the test
with `asPearson` set to true and to false to see what version is faster for your data. 

!!! note "nota"
    If `asPearson` is true, the `.stat` field of the test result will actually be `CrossProd()`,
    as the data will be standardized before running the test. See [`correlationTest`](@ref).

*Directional tests, permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`studentTestIS`](@ref)

*Aliases:* `tTestMcIS`, [`pointBiSerialMcTest`](@ref)

Return a [MultcompTest](@ref) structure.

**METHOD (2)**

As method (1), but input data `Yvec` is a vector containing ``M`` pairs of vector of arbitrary length, 
holding in the natural order the data corresponding to the ``m^{th}`` hypothesis for group 1 and group 2, respectively.

**METHOD (3)**

As method (1), but input data `X` and `Y` holds ``M`` vectors of observations each, `X` corresponding to data 
for group 1 and `Y` corresponding to data for group 2. 
The ``M`` vectors in `X` must have all the same length (``N_1``), so must the ``M`` vectors in `Y` (``N_2``). 

*Examples*

```julia
# (1)
using PermutationTests
M=100 # number of hypotheses
ns=[4, 5] # number of observations in group 1 and group 2 (N1 and N2)
N=sum(ns) # total number of observations
Yvec = [[randn(n) for n in ns] for m=1:M]; # some random Gaussian data for example 
Y=[vcat(yvec...) for yvec in Yvec];
t = tMcTestIS(Y, ns) # by default the test is bi-directional

# Force an approximate test with 10000 random permutations
tapprox = tMcTestIS(Y, ns; switch2rand=1, nperm=10000) 
tR=tMcTestIS(Y, ns; direction=Right()) # right-directional test
tL=tMcTestIS(Y, ns; direction=Left()) # left-directional test

# with a bi-directional test, t is equivalent to a 1-way ANOVA for independent samples
tanova= fMcTestIS(Y, ns) 
println(sum(abs.(t.p - tanova.p)) ≈ 0. ? "OK" : "error")

# do not run it using the CrossProd test statistic
tcor = tMcTestIS(Y, ns; asPearson=false) 

# in method (2) only the way the input data is formatted is different 
t2 = tMcTestIS(Yvec)
println(sum(abs.(t.p - t2.p)) ≈ 0. ? "OK" : "error")

# in method (3) also, only the way the input data is formatted is different 
t3 = tMcTestIS([Yvec[m][1] for m=1:M], [Yvec[m][2] for m=1:M])
println(sum(abs.(t.p - t3.p)) ≈ 0. ? "OK" : "error")
```

**Similar tests**

See [`studentTest1S`](@ref)
"""
function studentMcTestIS(𝐗::UniDataVec, 𝐘::UniDataVec;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)), 
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4,
                asPearson::Bool = true) where TestDir <: TestDirection

    length(𝐗)==1 && length(𝐘)==1 && (return studentTestIS(𝐗[1], 𝐘[1]; direction, equivalent=true, switch2rand, nperm, seed, verbose, asPearson)) 
    length(𝐗)==length(𝐘) || throw(ArgumentError(📌*"Function studentMcTestIS: the first two arguments must contain the same number of vectors"))
    length(unique(length.(𝐗)))==1 || throw(ArgumentError(📌*"Function studentMcTestIS: the vectors in the first argument (𝐗) must have all equal length. Check the documentation"))
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function studentMcTestIS: the vectors in the second argument (𝐘) must have all equal length. Check the documentation"))

    ns=[length(𝐗[1]), length(𝐘[1])]
    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱 = membership(StudentT_IS(), ns; rev=reverse)
        return _permMcTest!(μ0σ1(𝐱), [μ0σ1([𝐱; 𝐲]) for (𝐱, 𝐲) ∈ zip(𝐗, 𝐘)], length(𝐱), CrossProd(), PearsonR(); 
                standardized=true, stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    else
        𝐱 = membership(StudentT_IS(), ns)
        return _permMcTest!(𝐱, [[𝐱; 𝐲] for (𝐱, 𝐲) ∈ zip(𝐗, 𝐘)], ns, StudentT_IS(), StudentT_IS(); 
                stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    end
end

# alias
tMcTestIS = studentMcTestIS

"""
```julia
function pointBiSerialMcTest(<same args and kwargs as `studentMcTestIS`>)
```
Actually an alias for [`studentMcTestIS`](@ref).

Run ``M`` [point bi-serial correlation tests](https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient) simultaneously.
The correlations are between the ``M`` input data vectors ``y_1,...,y_M`` given as argument `Y`, 
all holding ``N=N_1+N_2`` elements, and a vector ``x``, internally created, with the first ``N_1`` elements equal to `1`
and the remaining ``N_2`` elements equal to `2`.

If you need to use other values for the dicothomous variable ``x`` or a different order for its elements, 
use [`correlationMcTest`](@ref) instead. 

The ``M`` null hypotheses have form 

``H_0(m): b_{(x,y_m)}=0, \\quad m=1...M``,

where ``b_{(x,y_m)}`` is the point bi-serial correlation between ``y_m`` (the ``m^{th}`` input data vectors in `Y`) 
and the internally created vector ``x``.

*Directional tests, permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`pointBiSerialTest`](@ref)

Return a [MultcompTest](@ref) structure.

*Examples*

```julia
using PermutationTests
ns=[4, 6]; # N1=4, N2=6
N=sum(ns); # number of observations

Y = [rand(N) for m=1:M]; # some Gaussian data as example
# implicitly, the point bi serial correlation is between 
# y1,...,yM and x=[1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
t=pointBiSerialMCTest(Y, ns) # by default the test is bi-directional
```
```julia
tR=pointBiSerialMCTest(Y, ns; direction=Right()) # right-directional test
tL=pointBiSerialMCTest(Y, ns; direction=Left()) # left-directional test
```
"""
# alias
pointBiSerialMcTest = studentMcTestIS



### Chi-Square and Fisher Exact Test
"""
```julia
function chiSquaredMcTest(tables::AbstractVector{Matrix{I}};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,            
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4,
            asPearson::Bool = true) 
                    where {I <: Int, TestDir <: TestDirection}
```

Multiple comparisons [chi-squared](https://en.wikipedia.org/wiki/Chi-squared_test) (``\\chi^2``) permutation test 
for ``2 \\cdot K`` contingency tables, where ``K`` is ≥2. 
It tests simultaneously ``M`` contingency tables, which must all have same dimension and same column sums.

The ``M`` null hypotheses have form 

``H_0(m): O_m=E_m``, \\quad m=1...M``,

where ``O_m`` and ``E_m`` are the observed and expected frequencies of the ``m^{th}``contingency table.

`tables` is a vector of ``M`` contingency tables. See [`chiSquaredTest`](@ref) for more explanations.

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests").

For ``K=2`` this function calls [`studentMcTestIS`](@ref) and pass to it also argument `asPearson`,
otherwise calls [`anovaMcTestIS`](@ref) and argument `asPearson` is ignored.

*Directional tests, permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`chiSquaredTest`](@ref)

*Aliases:* `Χ²McTest`, [`fisherExactMcTest`](@ref)

Return a [MultcompTest] structure.

!!! warning
    For ``K>2``, permutations of dicothomous tables may yield a null "sum of squares within", 
    thus an infinite *F* statistic for the ANOVA. In this case the `.nulldistr` field of the 
    returned structure will contain some julia `Inf` elements. This does not apply for the univariate 
    version of the test ([`chiSquaredTest`](@ref)), as in this case an equivalent statistic for the 
    ANOVA is used (see [Statistic](@ref)) and those statistics can never go to infinity. 

*Examples*
```julia
using PermutationTests
tables=[[3 3 2; 0 0 1], [1 2 2; 2 1 1], [3 1 2; 0 2 1]];
t=Χ²McTest(tables) # the test is bi-directional

tables=[[6 1; 2 5], [4 1; 4 5], [5 2; 3 4]]
tR=fisherExactMcTest(tables; direction=Right()) 
# or tR=Χ²Test(tables; direction=Right())

# do not use PearsonR statistic
tR_=fisherExactMcTest(tables; direction=Right(), asPearson=false) 

```
"""
function chiSquaredMcTest(tables::AbstractVector{Matrix{I}};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(tables), Int(1e4)), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,            
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4,
            asPearson::Bool = true) where {I <: Int, TestDir <: TestDirection}

    length(tables)==1 && (return chiSquaredTest(tables[1]; direction, equivalent=true, switch2rand, nperm, seed, verbose))  
    !(direction isa Both) && (size(tables[1], 2) > 2) && throw(ArgumentError(📌*"Function chiSquaredMcTest or fisherExactMcTest: For input data matrices with more than 2 columns the test can only be bi-directional. Correct or eliminate the `direction` keyword argument"))
    𝐘, ns = table2vec(tables, AnovaF_IS()) # anovaTestIS will switch to t-test if there are only two groups
    length(ns) < 2 && throw(ArgumentError(📌*"Function chiSquaredMcTest or fisherExactMcTest: the second argument must be a vector of two or more integers"))

    return length(ns)>2 ?   anovaMcTestIS(𝐘, ns; direction, switch2rand, nperm, seed, verbose, stepdown, fwe, threaded) :
                            studentMcTestIS(𝐘, ns; direction, switch2rand, nperm, seed, verbose, stepdown, fwe, threaded, asPearson)
end

# alias
Χ²McTest = chiSquaredMcTest


"""
```julia
function fisherExactMcTest(<same args and kwargs as `chiSquaredMcTest`>)
```
Actually an alias for [`chiSquaredMcTest`](@ref). It can be used for 2x2 contingency tables.
See [`chiSquaredTest`](@ref) for more explanations.

*Univariate version:* [`fisherExactTest`](@ref)

Return a [MultcompTest](@ref) structure.

*Examples*

```julia
using PermutationTests
tables=[[6 1; 2 5], [4 1; 4 5], [5 2; 3 4]];
tR=fisherExactMcTest(tables; direction=Right()) 
# or tR=Χ²Test(tables; direction=Right())
```
"""
# alias
fisherExactMcTest = chiSquaredMcTest


### ANOVA for Repeated Measures
"""
```julia
# METHOD (1)
function anovaMcTestRM(Y::UniDataVec, ns::@NamedTuple{n::Int, k::Int};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection
```
"""
function anovaMcTestRM(𝐘::UniDataVec, ns::@NamedTuple{n::Int, k::Int};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘) == 1 && (return anovaTestRM(𝐘[1], ns; direction, equivalent=true, switch2rand, nperm, seed, verbose))
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function anovaMcTestRM: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))
    ns.k < 2 && throw(ArgumentError(📌*"Function anovaMcTestRM: the `k` element of second argument (a named tuple) must be 2 or more"))
    ns.n < ns.k && throw(ArgumentError(📌*"Function anovaMcTestRM: the `n` element of second argument (a named tuple) must be greater than the `k` element"))
    #ns.n % ns.k ≠ 0 && throw(ArgumentError("Function anovaMcTestRM: impossible configuration for the `n` and `k` elements of second argument (a named tuple). The former is the number of subjects, the latter the treatments."))
    ns.n*ns.k ≠ length(𝐘[1]) && throw(ArgumentError("📌*Function anovaMcTestRM: the `n` and `k` elements of second argument (a named tuple) must be such that their product is equal to the length of the first argument (𝐲)"))
    
    if ns.k == 2 # do t-test if there are only two groups. Actually do one-sample test on the difference 
        𝐱 = membership(StudentT_1S(), ns.n)
        return _permMcTest!(𝐱, [𝐲[1:2:ns.n*2-1].-𝐲[2:2:ns.n*2] for 𝐲 ∈ 𝐘], ns.n, StudentT_1S(), StudentT_1S(); stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    else 
        𝐱 = membership(AnovaF_RM(), ns)
        !(direction isa Both) && throw(ArgumentError("📌*Function anovaMcTestRM: The ANOVA test can only be bi-directional. Correct the `direction` keyword argument")) 
        return _permMcTest!(𝐱, 𝐘, ns, AnovaF_RM(), AnovaF_RM(); stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    end
end

"""
```julia
# METHOD (2)
function anovaMcTestRM(Yvec::UniDataVec²; <same kwargs>)
```

**METHOD (1)**

Multiple comparison [1-way analysis of variance (ANOVA) for repeated measures](https://en.wikipedia.org/wiki/Repeated_measures_design#Repeated_measures_ANOVA) 
by data permutation. Given ``M`` hypotheses, each with ``N`` observation units (e.g., *subjects*, *blocks*, etc.) 
for each of ``K`` repeated measures (e.g., *treatments*, *time*, etc.), the null hypotheses have form

``H_0(m): μ_{m1}= \\ldots =μ_{mk}, \\quad m=1...M``,

where ``μ_{mk}`` is the mean of the ``k^{th}`` treatment for the ``m^{th}`` hypothesis. 

`Y` is a vector hoding ``M`` vectors, each one concatenaning the observations for the ``K`` treatments 
in the natural order, that is, the ``N`` observation for treatment 1, ..., 
the ``N`` observations for treatment ``K``. Thus, `Y` holds ``M`` vectors of ``N \\cdot K`` elements. 

`ns` is a julia [named tuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) 
with form `(n=N, k=K)` (see examples below).

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests"). 

*Directional tests, permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`anovaTestRM`](@ref)

*Alias:* `fTestMcRM`

Return a [MultcompTest](@ref) structure.

**METHOD (2)**

As (1), but `Yvec` is a vector of ``M`` vectors, each one holding the ``K`` vectors of ``N`` observations.

*Examples*

```julia
# method (1)
using PermutationTests
N=6 # number of observation units per treatment
K=3 # number of treatments
M=10 # number of hypotheses
Yvec = [[randn(N) for k=1:K] for m=1:M]; # some random Gaussian data for example 
t = fMcTestRM([vcat(yvec...) for yvec in Yvec], (n=N, k=K)) 
# ANOVA tests are always bi-directional
```
```julia
# Force an approximate test with 5000 random permutations
tapprox = fMcTestRM([vcat(yvec...) for yvec in Yvec], (n=N, k=K); 
            switch2rand=1, nperm=5000)

# in method (2) only the way the input data is formatted is different 
 
t2 = fMcTestRM(Yvec)
println(sum(abs.(t.p - t2.p)) ≈ 0. ? "OK" : "error")

```

**Similar tests**

See [`anovaTestRM`](@ref)
"""
function anovaMcTestRM(𝐘vec::UniDataVec²;
                direction::TestDir = Both(),
                switch2rand::Int = max(Int(1e8) ÷ length(𝐘vec), Int(1e4)),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                #
                stepdown::Bool = true,
                fwe::Float64 = 0.05,
                threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    K = length(𝐘vec[1])
    K < 2 && throw(ArgumentError(📌*"Function anovaMcTestRM: the first (𝐘vec) argument must be a vector of vector of two or more vectors"))
    N = length(𝐘vec[1][1])
    N < K && throw(ArgumentError(📌*"Function anovaMcTestRM: the length of the vectors in first argument must be greater than their number"))
    ns=(n=N, k=K)
    length(𝐘vec) == 1 && (return anovaTestRM(𝐘vec[1], ns; direction, equivalent=true, switch2rand, nperm, seed, verbose))
    unique(length(y) for y in 𝐘vec)[1]==K || throw(ArgumentError(📌*"Function anovaMcTestRM: All vectors in first arguments (𝐘vec) must contain the same number of vectors (K). Check the documentation")) 
    unique(length.(y) for y in 𝐘vec)[1]==repeat([N], K) || throw(ArgumentError(📌*"Function anovaMcTestIS: All vectors in first arguments (𝐘vec) must contain vectors of equal length (N). Check the documentation"))

    if K == 2 # do t-test if there are only two groups. Actually do one-sample test on the difference         
        𝐱 = membership(StudentT_1S(), ns.n)
        return _permMcTest!(𝐱, [𝐲[1].-𝐲[2] for 𝐲 ∈ 𝐘vec], ns.n, StudentT_1S(), StudentT_1S(); stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    else 
        !(direction isa Both) && throw(ArgumentError(📌*"Function anovaMcTestRM: The ANOVA test can only be bi-directional. Correct the `direction keyword argument`")) 
        𝐱 = membership(AnovaF_RM(), ns)
        return _permMcTest!(𝐱, [vcat(𝐲...) for 𝐲 ∈ 𝐘vec], ns, AnovaF_RM(), AnovaF_RM(); stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
    end
end

# alias
fMcTestRM = anovaMcTestRM 


### Cocharn Q test

# NB: not particular efficient. For an efficient alternative see https://dl.acm.org/doi/10.1145/3095076
"""
```julia
function cochranqMcTest(tables::AbstractVector{Matrix{I}};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) 
                where {I <: Int, TestDir <: TestDirection}
```

Multiple comparisons [Cochran Q](https://en.wikipedia.org/wiki/Cochran's_Q_test) by data permutation. 

The Cochran Q test is analogous to the 1-way ANOVA for repeated measures, but takes as input dicothomous 
data (zeros and ones). Given ``M`` hypotheses, consisting each in ``N`` observation units 
(e.g., *subjects*, *blocks*, etc.) and ``K`` repeated measures (e.g., *treatments*, *time*, etc.), 
the null hypotheses have form

``H_0(m): μ_{m1}= \\ldots =μ_{mk}, \\quad m=1...M``,

where ``μ_{mk}`` is the mean of the ``k^{th}`` treatment. 

Input `tables` is a vector of ``M`` tables of zeros and ones with size ``NxK``, where ``N`` is the number 
of observations and ``K`` the repeated measures. See [`cochranqTest`](@ref) for more explanations.

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests"). 

*Directional tests, Permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`cochranqTest`](@ref)

*Aliases:* `qMcTest`, [`mcNemarMcTest`](@ref)

Return a [MultcompTest](@ref) structure.

*Examples*
```julia
using PermutationTests
tables=[ [1 1 0; 1 0 0; 1 1 1; 1 1 0; 1 0 1; 1 1 0],
        [1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 0 0; 1 0 0],
        [1 0 0; 0 0 1; 1 0 1; 1 1 0; 1 0 1; 1 0 0]];
t=qMcTest(tables) # the test with K>2 can only be bi-directional
```
"""
function cochranqMcTest(tables::AbstractVector{Matrix{I}};
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(tables), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where {I <: Int, TestDir <: TestDirection}

    length(tables)==1 && (return cochranqTest(tables[1]; direction, equivalent=true, switch2rand, nperm, seed, verbose))
    𝐘, ns = table2vec(tables, AnovaF_RM()) # anovaTestRM will switch to t-test if there are only two groups
    return anovaMcTestRM(𝐘, ns; direction, switch2rand, nperm, seed, verbose, stepdown, fwe, threaded)
end

# alias
qMcTest = cochranqMcTest 

### McNemar test
"""
```julia
function mcNemarMcTest(same args and kwargs as `cochranqMcTest`>)
```
Actually an alias for [cochranqMcTest)(@ref).

Run ``M`` McNemar test simultaneously.

Input `tables` is a vector of ``M`` tables of zeros and ones with size ``Nx2``, where ``N`` is the number 
of observations and ``2`` the number of repeated measures. See [`cochranqTest`](@ref) for more explanations.

*Univariate version:* [`mcNemarTest`](@ref)

Return a [MultcompTest](@ref) structure.

*Examples*
```julia
using PermutationTests
tables=[[1 1; 1 0; 1 0; 0 0; 1 0; 1 0],
        [1 0; 1 1; 1 0; 0 1; 0 0; 1 0],
        [0 1; 0 0; 1 0; 1 0; 1 0; 1 1]];
t=mcNemarMcTest(tables) # by default the test is bi-directional

tR=mcNemarMcTest(tables; direction=Right()) # right-directional test
```
"""
mcNemarMcTest = cochranqMcTest # alias


### One-Sample t-test
"""
```julia
function studentMcTest1S(Y::UniDataVec;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection 
```

Multiple comparison [one-sample t-test](https://en.wikipedia.org/wiki/Student's_t-test#One-sample_t-test) 
by data permutation. 

Run ``M`` one-sample t-tests simultaneously. The null hypotheses have form 

``H_0(m): μ_m=μ_0, \\quad m=1...M``,

where ``μ_m`` is the mean of the observations for the ``m^{th}`` hypothesis and ``μ_{0}`` is a reference 
population mean.

`refmean` is the reference mean (``μ_0``) above. The default is `0.0`, which is the value needed in most situations.

`Y` is a vector of ``M`` vectors holding each the ``N`` observations which mean is to be compared to ``μ_0``.

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests"). 

*Directional tests, permutation scheme and number of permutations for exact tests:* as per
*univariate version* [`studentTest1S`](@ref)

*Alias:* `tTest1S` 

Return a [MultcompTest](@ref) structure.

*Examples*
```julia
using PermutationTests
N=20 # number of observations
M=100 # number of hypotheses
y = [randn(N) for m=1:M]; # some random Gaussian data for example
t = tMcTest1S(Y) # By deafult the test is bi-directional

tR = tMcTest1S(Y; direction=Right()) # right-directional test
tL = tMcTest1S(Y; direction=Left()) # Left-directional test

# Force an approximate test with 5000 random permutations
tapprox = tMcTest1S(Y; switch2rand=1, nperm=5000) 

# test H0m: μ(ym)=1.5: all will be rejected as the expected mean is 0.0
t1 = tMcTest1S(Y; refmean=1.5) 
```

**Similar tests**

See [`studentTest1S`](@ref)
"""
function studentMcTest1S(𝐘::UniDataVec;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘)==1 && (return stedentTest1S(𝐘[1]; refmean, direction, equivalent=true, switch2rand, nperm, seed, verbose))       
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function studentMcTest1S: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))

    _permMcTest!(membership(StudentT_1S(), length(𝐘[1])), refmean===nothing ? copy(𝐘) : [𝐲.-refmean for 𝐲 in 𝐘], length(𝐘[1]), StudentT_1S(), StudentT_1S();
    stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
end

tMcTest1S = studentMcTest1S # aliases    


# method (2) as studentMcTest1S!, but copy the 𝐘 vectors, so it does not overwrite them. 
# This method directly subtracts the reference mean    
"""
```julia
function studentMcTest1S!(<same args and kwargs as `studentMcTest1S`>)
```

As [`studentMcTest1S`](@ref), but `Y` is overwritten in the case of approximate (random permutations) tests.

*Alias:* `tTest1S!` 

Univariate version: [`studentTest1S!`](@ref)

"""
function studentMcTest1S!(𝐘::UniDataVec;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            switch2rand::Int = max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection 

    length(𝐘)==1 && (return stedentTest1S!(𝐘[1]; refmean, direction, equivalent=true, switch2rand, nperm, seed, verbose))       
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function studentMcTest1S!: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))

    _permMcTest!(membership(StudentT_1S(), length(𝐘[1])), refmean===nothing ? 𝐘 : [𝐲.-refmean for 𝐲 in 𝐘], length(𝐘[1]), StudentT_1S(), StudentT_1S();
    stepdown, fwe, nperm, fstat=_fstat(direction), switch2rand, seed, threaded, verbose)
end 

tMcTest1S! = studentMcTest1S!

### Sign Test
# Takes as input K vectors of vectors of N booleans, where N is the number of subjects.
"""
```julia
function signMcTest(Y::Union{AbstractVector{BitVector}, AbstractVector{Vector{Bool}}};
            direction::TestDir = Both(),
            switch2rand::Int= max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

```

Multiple comparisons [sign test](https://en.wikipedia.org/wiki/Sign_test) by data permutation. 

Run ``M`` sign tests simultaneously.The null hypotheses have form 

``H_0(m): E_m(true)=E_m(false), \\quad m=1...M``,

where ``E_m(true)`` and ``E_m(false)`` are the expected number of true and false occurrences, respectively,
in the ``m^{th}`` hypothesis.

`Y` ia a vector of ``M`` vectors holding each ``N`` booleans.

For optional keyword arguments, `direction`, `switch2rand`, `nperm`, `seed`, `verbose`, `stepdown`, `fwe`
and `threaded` see [here](@ref "Common kwargs for multiple comparisons tests"). 

*Directional tests, permutation scheme and number of permutations for exact tests:* as per 
*univariate version* [`signTest`](@ref)

Return a [MultcompTest](@ref) structure.

*Examples*

```julia
using PermutationTests
N=20; # number of observations
M=100; # number of hypotheses
Y = [rand(Bool, N) for m=1:M]; # some random Gaussian data for example
t = signMcTest(Y) # By deafult the test is bi-directional

tR = signMcTest(Y; direction=Right()) # right-directional test
tL = signMcTest(Y; direction=Left()) # Left-directional test

# Force an approximate test with 5000 random permutations
tapprox = signMcTest(Y; switch2rand=1, nperm=5000) 
```
"""
function signMcTest(𝐘::Union{AbstractVector{BitVector}, AbstractVector{Vector{Bool}}};
            direction::TestDir = Both(),
            switch2rand::Int= max(Int(1e8) ÷ length(𝐘), Int(1e4)),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            #
            stepdown::Bool = true,
            fwe::Float64 = 0.05,
            threaded::Bool = Threads.nthreads()>=4) where TestDir <: TestDirection

    length(𝐘)==1 && (return signTest(𝐘[1]; direction, equivalent=true, switch2rand, nperm, seed, verbose))       
    length(unique(length.(𝐘)))==1 || throw(ArgumentError(📌*"Function signMcTest: the vectors in the first argument (𝐘) must have all equal length. Check the documentation"))

    studentMcTest1S!([Float64.(𝐲).-0.5 for 𝐲 ∈ 𝐘]; refmean=0., direction, switch2rand, nperm, seed, verbose, stepdown, fwe, threaded)
end

"""
```julia
function studentMcTestRM(<same args and kwargs as `studentMcTest1S`>)
```

Actually an alias for [`studentMcTest1S`](@ref)
    
In order to run a multiple comparisons t-test for repeated measure, use as data input the vector of 
``M`` vectors of differences across the two measurements (or treatments, time, etc.).

Do not change the `refmean` default value. See [`studentMcTest1S`](@ref) for more details.

*Directional tests, permutation scheme and number of permutations for exact tests:* as per [`studentTest1S`](@ref)

*Univariate version:* [`studentTestRM`](@ref)

*Alias:* `tMcTestRM`

Return a [MultcompTest](@ref) structure.

*Examples*

```julia
using PermutationTests
N=10 # number of observation per treatment
M=100 # number of hypotheses

# suppose you have data as
Y1=[randn(N) for m=1:M]; # measurement 1
Y2=[randn(N) for m=1:M]; # measurement 2

# Let us compute the differences as
Y=[y1-y2 for (y1, y2) in zip(Y1, Y2)];
t=tMcTestRM(Y) # by default the test is bi-directional

tR=tMcTestRM(Y; direction=Both()) #  right-directional test
# if test tR is significant for some hypotheses, 
# for these hypotheses the mean of measurement 1 
# exceeds the mean of measurement 2.
```
"""
studentMcTestRM = studentMcTest1S

tMcTestRM = studentMcTest1S


"""
(2)
function studentMcTestRM!(<same args and kwargs as `studentMcTest1S!`>)

Actually an alias for [`studentMcTest1S!`](@ref).

`Y` is overwritten in the case of approximate (random permutations) tests.

*Alias:* `tMcTestRM!`

"""
studentMcTestRM! = studentMcTest1S!

tMcTestRM! = studentMcTest1S!

# ----------------------------------------------------- #
