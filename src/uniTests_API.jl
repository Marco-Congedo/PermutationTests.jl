#=
uniTests_API.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home

########################################
# API for Univariate Permutation tests #
########################################

Run `test_uniTests_API()` in unit _test.jl to test this unit

========================================
  EXPORTED:

correlationTest!, rTest!,
correlationTest,  rTest,
trendTest!, trendTest,
pointBiSerialTest,
anovaTestIS, fTestIS,
studentTestIS, tTestIS, 
chiSquaredTest, Χ²Test, 
fisherExactTest,
anovaTestRM, fTestRM,
cochranqTest, qTest,
mcNemarTest,
studentTestRM, tTestRM,
studentTestRM!, tTestRM!,
studentTest1S!, tTest1S!,
studentTest1S, tTest1S,
signTest!, signTest,

  UTILITIES:
- _permutationTest
- _permutationTest!
========================================
=#

# GENERAL KEYWORD ARGUMENTS
# direction: test direction, either Right(), Left() or Both(). Default: Both()
# equivalent (Bool): if true, the fastest equivalent statistic will be used. Default: true
# switch2rand (Int): the # of permutations starting from which the approximate test will be carried out. Default: 1e8
# nperm (Int): the # of permutations to be used if the test will be approximate (default: 20_000) 
# seed (Int): random # generator seed. Set to 0 to take a random seed. Any natural number allow a reproducible test.(Default: 1234) 
# verbose (Bool) print some information in the REPL while running the test. Set to false for running @benchmark (Default: true).


#####
## This function ultimately run all tests. No check on data input is performed

# This method run all tests BUT the correlation and trend tests and the 1 sample tests - 9 arguments        
function _permutationTest(𝐲::UniData, ns::Union{Vector{Int}, @NamedTuple{n::Int, k::Int}}, stat, direction, equivalent, switch2rand, nperm, seed, verbose)
    𝐱 = membership(stat, ns)
    design = assignment(stat, ns)
    eqstat = equivalent ? eqStat(stat, direction, design) : stat
    return _permTest!(𝐱, 𝐲, ns, eqstat, eqstat; 
                fstat=_fstat(direction), switch2rand, nperm, seed, verbose)      
end

# This run only correlation and trend tests. Don't use for other tests! - 10 arguments
function _permutationTest!(𝐱::UniData, 𝐲::UniData, standardized, centered, direction, equivalent, switch2rand, nperm, seed, verbose)
    centered && direction==Both() && (equivalent=true)
    if standardized 
        eqstat = CrossProd()
    else
        eqstat = equivalent ? eqStat(PearsonR(), direction, Balanced()) : PearsonR()
    end
    return _permTest!(𝐱, 𝐲, length(𝐱), eqstat, eqstat; 
                standardized, centered, fstat=_fstat(direction), switch2rand, nperm, seed, verbose)
end

# This run only one-sample t-tests, hence repeated-measures t-tests. Don't use for other tests!
# NB: for exact tests 𝐱 is rewritten, but for monte carlo tests input data vector 𝐲 is rewritten - 8 arguments 
function _permutationTest!(𝐲::UniData, ns::Int, direction, equivalent, switch2rand, nperm, seed, verbose)
    𝐱 = membership(StudentT_1S(), ns)
    eqstat = equivalent ? eqStat(StudentT_1S(), direction, Balanced()) : StudentT_1S()
    return _permTest!(𝐱, 𝐲, ns, eqstat, eqstat; fstat=_fstat(direction), switch2rand, nperm, seed, verbose)     
end


### Correlation test 
"""
```julia
function correlationTest(x::UniData, y::UniData;
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            standardized::Bool = false, 
            centerd::Bool = false,
            verbose::Bool = true) where TestDir <: TestDirection
```

Univariate [Pearson product-moment correlation test](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Testing_using_Student's_t-distribution) 
by data permutation. The null hypothesis has form 

``H_0: r_{(x,y)}=0``,

where ``r_{(x,y)}`` is the correlation between the two input data vectors, `x` and `y`, typically real, both holding ``N`` observations.

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests").

If `standardized` is true, both `x` and `y` are assumed standardized (zero mean and unit standard deviation).
Provided that the input data is standardized, the test provides the same p-value, however it can be executed faster as in this case the cross-product
is equivalent to the Pearon *r* statistic (see [Statistic](@ref)).

If `centered` is true, both `x` and `y` are assumed centered (zero mean).
The test provides the same p-value, however it can be executed faster if the test is bi-directional 
as in this case the equivalent statistic, the covariance, reduces to the cross-product divided by N.

If neither `standardized` nor `centered` is true, the data will be standardized to execute a faster test using the
cross-product as test-statistic.

*Directional tests*

 - For a right-directional test, the correlation is expected to be positive. A negative correlation will result in a p-value higehr then 0.5.
 - For a left-directional test, the correlation is expected to be negative. A positive correlation will result in a p-value higehr then 0.5.

*Permutation scheme:* under the null hypothesis, the position of the observations in the data input vectors bears no meaning.
The exchangeability scheme consists then in shuffling the observations of vector `x` or vector `y`.
*PermutationTests.jl* shuffles the observations in `x`. 

*Number of permutations for exact tests:* there are ``N!`` possible ways of reordering the ``N`` observations in `x`.

*Aliases:* `rTest!`, [`trendTest!`](@ref)

*Multiple comparisons version:* [correlationMcTest!](@ref)

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
N=10 # number of observations
x, y = randn(N), randn(N) # some random Gaussian data for example
t = rTest(x, y) # by deafult the test is bi-directional
```

```@julia
tR = rTest(x, y; direction=Right()) # right-directional test
tL = rTest(x, y; direction=Left()) # Left-directional test
# Force an approximate test with 5000 random permutations
tapprox = rTest(x, y; switch2rand=1, nperm=5000) 
```

**Similar tests**

Typically, the input data is real, but can also be of type integer or boolean. If either `x` or `y` is a vector 
of booleans or a vector of dicothomous data (only 0 and 1), this function will actually perform a permutation-based 
version of the [point bi-serial correlation test](https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient).
However, as shown in the preceeding link, the point bi-serial correlation test is equivalent to 
the t-test for independent sample, thus it can be tested using the t-test for independent samples,
which will need many less permutations as compared to a correlation test for an exact test 
(see examples below). A dedicated function in available with name [`pointBiSerialTest`](@ref),
which is an alias for [`studentTestIS`](@ref) and allowa the choice to run the test using a 
correlation- or t-test statistic.

If `x` or `y` represent a *trend*, for example a linear trend given by `[1, 2,...N]`, 
we otain the permutation-based trend correlation test, which can be used to test the fit of any type of regression
of `y` on `x` - see [trendTest](@ref).    

if `y` is a shifted version of `x` with a lag ``l``, this function will test the significance of the 
*autocorrelation at lag ``l``, see the page [Create your own test](@ref).

*Examples*
```julia
# Point bi-serial correlation test
using PermutationTests
N=10 # number of observations
x=[0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
y = rand(N)
t = rTest(x, y) 

# Exactly the same test can be obtained as a t-test for independent sample,
# but much faster as for an exact test the latter needs only 210 permutations 
# while the former needs 3628800 permutations.
# This is available with a dedicated function
t2=pointBiSerialTest(y, [4, 6])
println(t.p ≈ t2.p ? "OK" : "error")
```
"""
correlationTest(𝐱::UniData, 𝐲::UniData;
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            standardized::Bool = false, 
            centered::Bool = false,
            verbose::Bool = true) where TestDir <: TestDirection =
    correlationTest!((standardized || centered) ? copy(𝐱) : 𝐱, 𝐲; 
                    direction, equivalent, switch2rand, nperm, seed, standardized, centered, verbose)

# alias. The trendTest is obtained passed the trend to be tested as vector `𝐱`. For example [1, 2,...] for an ascending linear trend                        
rTest = correlationTest


"""
```julia
function correlationTest!(<same args and kwargs as `correlationTest`>)
```

As [`correlationTest`](@ref), but `x` is overwritten if neither `standardized` nor `centered` is true.

*Aliases:* `rTest!`, [`trendTest!`](@ref) 

*Multiple comparisons version:* [correlationMcTest!](@ref)
"""
function correlationTest!(𝐱::UniData, 𝐲::UniData;
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8), 
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            standardized::Bool = false, 
            centered::Bool = false,
            verbose::Bool = true) where TestDir <: TestDirection
 
    length(𝐱) == length(𝐲) || throw(ArgumentError(📌*"Function correlationTest!: the first two arguments are the vectors of the two variables and must have equal length."))
    if standardized 
        return _permutationTest!(𝐱, 𝐲, true, centered, direction, equivalent, switch2rand, nperm, seed, verbose)
    elseif centered
        return _permutationTest!(𝐱, 𝐲, false, true, direction, equivalent, switch2rand, nperm, seed, verbose)
    else
        _permutationTest!(μ0σ1(𝐱), μ0σ1(𝐲), true, false, direction, equivalent, switch2rand, nperm, seed, verbose)
    end
end

# alias
rTest! = correlationTest! 

# same as correlationTest! but copy vector 𝐱, thus this will not be overwritten.
"""
```julia
function trendTest(<same args and kwargs as `correlationTest`>)
```
"""
trendTest = correlationTest

"""
```julia
function trendTest!(<same args and kwargs as `correlationTest!`>)
```
Actually aliases for [`correlationTest`](@ref) and [`correlationTest!`](@ref), respectively.

The two vectors `x` and `y` of ``N`` elements each are provided as data input. 
`x` is a specified trend and `y` holds the observed data. A Pearson product-moment correlation test between 
`x` and `y` is then carried out. 

`x` can hold any trend, such as linear, polynomial, exponential, logarithmic, power, trigonometric...

*Directional tests, permutation scheme and number of permutations for exact tests:* as per [`correlationTest`](@ref)

Multiple comparisons versions: [`trendMcTest`](@ref) and [`trendMcTest!`](@ref)

Both methods return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
# We are goint to test an upward linear trend
N=10
x=Float64.(collect(Base.OneTo(N))) # [1, 2, ..., N]
y=[1., 2., 4., 3., 5., 6., 7., 8., 10., 9.]
# Supposing we expect an upward linear trend, 
# hence the correlation is expected to be positive,
# we can use a right-directional test to increase the power of the test.
t = trendTest(x, y; direction=Right()) 
```
""" 
trendTest! = correlationTest!


### ANOVA for independent samples
"""
```julia
# METHOD (1)
function anovaTestIS(y::UniData, ns::IntVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true) where TestDir <: TestDirection
```
"""
function anovaTestIS(𝐲::UniData, ns::IntVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true) where TestDir <: TestDirection

    length(ns) < 2 && throw(ArgumentError(📌*"Function anovaTestIS: the second argument (ns) must be a vector of two or more integers"))
    sum(ns) ≠ length(𝐲) && throw(ArgumentError(📌*"Function anovaTestIS: the sum of the integers in second argument (ns) must be equal to the length of the first argument (𝐲)"))
    stat = length(ns) == 2 ? StudentT_IS() : AnovaF_IS() # do t-test if there are only two groups
    stat isa AnovaF_IS && !(direction isa Both) && throw(ArgumentError(📌*"Function anovaTestIS: The ANOVA test can only be bi-directional. Correct the `direction` keyword argument")) 
    return _permutationTest(𝐲, ns, stat, direction, equivalent, switch2rand, nperm, seed, verbose)
end

"""
```julia
# METHOD (2)
function anovaTestIS(yvec::UniDataVec; <same kwargs>)
```

**METHOD (1)**

Univariate [1-way analysis of variance (ANOVA) for independent samples](https://en.wikipedia.org/wiki/One-way_analysis_of_variance) 
by data permutation. 
Given ``N=N1+...+N_K`` observations in ``K`` groups, the null hypothesis has form

``H_0: μ_1= \\ldots =μ_K``,

where ``μ_k`` is the mean of the ``k^{th}`` group. 

`y` is a vector concatenaning the vector of observations in each group, in the natural order. Thus, it holds ``N`` elements. 

`ns` is a vector of integers holding the group numerosity (see examples below).

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests"). 

*Directional tests*

Possible only for ``𝐾=2``, in which case the test reduces to a Student' t-test for independent samples
and the test directionality is given by keyword arguement `direction`. 
See function [`studentTestIS`](@ref) and its multiple comparisons version [`studentMcTestIS`](@ref). 

*Permutation scheme:* under the null hypothesis, the group membership of the observations bears no meaning.
The exchangeability scheme consists then in reassigning the ``N`` observations in the ``K`` groups 
respecting the original group numerosity. 

*Number of permutations for exact tests:* there are ``\\frac{N!}{N_1 \\cdot \\ldots \\cdot N_K}`` possible ways 
of reassigning the ``N`` observations in the ``K`` groups.

**Alias:** `fTestIS`

*Multiple comparisons version:* [anovaMcTestIS](@ref)

Both methods return a [UniTest](@ref) structure.

**METHOD (2)**

As (1), but `yvec` is a vector of K vectors of observations, of for each group.

*Examples*

```julia
# (1)
using PermutationTests
ns=[4, 5, 6] # number of observations in group 1, 2 and 3
yvec = [randn(n) for n in ns] # some random Gaussian data for example 
t = fTestIS(vcat(yvec...), ns) # ANOVA tests are always bi-directional
```
```julia
# Force an approximate test with 5000 random permutations
tapprox = fTestIS(vcat(yvec...), ns; switch2rand=1, nperm=5000) 
```

```julia
# in method (2) only the way the input data is formatted is different 
t2 = fTestIS(yvec)
println(t.p ≈ t2.p ? "OK" : "error")
```

**Similar tests**

Typically, for ANOVA the input data is real, but can also be of type integer or boolean. For dicothomous data, 
with this function one can obtain a permutation-based version of the [Χ² test](https://en.wikipedia.org/wiki/Chi-squared_test) 
for ``K \\cdot 2`` contingency tables, which has the ability to give exact p-values. 
For ``2 \\cdot 2`` contingency tables it yields exactly the same p-value of
the [Fisher exact test](https://en.wikipedia.org/wiki/Fisher's_exact_test), which is also exact, as the name suggests.
In these cases it is more convenient to use the [chiSquaredTest](@ref) and [fisherExactTest](@ref) functions though, 
which accept contingency tables as data input. 

"""
function anovaTestIS(𝐲vec::UniDataVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true) where TestDir <: TestDirection

    K = length(𝐲vec)
    K < 2 && throw(ArgumentError(📌*"Function anovaTestIS: the first argument (𝐲vec) must be a vector of two or more vectors"))
    #N = sum(length(𝐲) for 𝐲 ∈ 𝐲vec)
    stat = K == 2 ? StudentT_IS() : AnovaF_IS() # do t-test if there are only two groups
    stat isa AnovaF_IS && !(direction isa Both) && throw(ArgumentError(📌*"Function anovaTestIS: The ANOVA test can only be bi-directional. Correct the `direction keyword argument`")) 
    ns=[length(𝐲) for 𝐲 in 𝐲vec] 
    return _permutationTest(vcat(𝐲vec...), ns, stat, direction, equivalent, switch2rand, nperm, seed, verbose)
end

# alias
fTestIS = anovaTestIS


### t-test for independent samples
"""
```julia
# METHOD (1)
function studentTestIS(y::UniData, ns::IntVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                asPearson::Bool = true) where TestDir <: TestDirection
```
"""
function studentTestIS(𝐲::UniData, ns::IntVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                asPearson::Bool = true) where TestDir <: TestDirection
    length(ns) ≠ 2 && throw(ArgumentError(📌*"Function studentTestIS: the second argument (ns) must be a vector of two integer"))
    sum(ns) ≠ length(𝐲) && throw(ArgumentError(📌*"Function studentTestIS: the sum of the two integers in second argument (ns) must be equal to the length of the first argument (𝐲)"))
    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱 = membership(StudentT_IS(), ns; rev=reverse)
        return correlationTest(𝐱, 𝐲;
            direction, equivalent, switch2rand, nperm, seed, standardized=false, centered=false, verbose)
    else
        return _permutationTest(𝐲, ns, StudentT_IS(), direction, equivalent, switch2rand, nperm, seed, verbose)
    end
end

"""
```julia
# METHOD (2)
function studentTestIS(yvec::UniDataVec; <same kwargs>)
```
"""
function studentTestIS(𝐲vec::UniDataVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true,
                asPearson::Bool = true) where TestDir <: TestDirection
    K = length(𝐲vec)
    K ≠ 2 && throw(ArgumentError(📌*"Function studentTestIS: the first argument (𝐲vec) must be a vector of two vectors"))
    ns=[length(𝐲) for 𝐲 in 𝐲vec] 
            
    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱 = membership(StudentT_IS(), ns; rev=reverse)
        return correlationTest(𝐱, vcat(𝐲vec...);
            direction, equivalent, switch2rand, nperm, seed, standardized=false, centered=false, verbose)
    else
        return _permutationTest(vcat(𝐲vec...), ns, StudentT_IS(), direction, equivalent, switch2rand, 
                        nperm, seed, verbose)
    end
end


"""
```julia
# METHOD (3)
function studentTestIS(y1::UniData, y2::UniData; <same kwargs>)
```

**METHOD (1)**

Univariate [Student's t-test for independent samples](https://en.wikipedia.org/wiki/Student's_t-test#Independent_(unpaired)_samples) 
by data permutation. Given ``N=N1+N_2`` observations in two groups, 
the null hypothesis has form

``H_0: μ_1=μ_2``,

where ``μ_1`` and ``μ_1`` are the mean for group 1 and group 2, respectively. 

For a bi-directional test, this t-test is equivalent to a 1-way ANOVA for two independent samples.
However, in contrast to the ANOVA, it can be directional. 

`y` is a vector concatenaning the vector of observations in the two groups. Thus, it holds ``N`` elements.

`ns` is a vector of integers holding the group numerosity (see examples below).

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests"). 

If `asPearson` is true (default), the test is run as an equivalent version of a Pearson correlation test.
This is not faster in general for exact univariate tests, since the t-test needs 
less permutations, but is in general advantageous for approximate tests (see the [benchmarks](@ref "Benchmarks")).
If your need to perform exact tests, you may want to set `asPearson` to false.

!!! note "nota"
    If `asPearson` is true, the `.stat` field of the test result will actually be `CrossProd()`,
    as the data will be standardized before running the test. See [`correlationTest`](@ref).

*Directional tests*

 - For a right-directional test, ``μ_1`` is expected to exceed ``μ_2``. If the opposite is true, the test will result in a p-value higehr then 0.5.
 - For a left-directional test, ``μ_2`` is expected to exceed ``μ_1``. If the opposite is true, the test will result a p-value higehr then 0.5.

*Permutation scheme:* under the null hypothesis, the group membership of the observations bears no meaning.
The exchangeability scheme consists then in reassigning the ``N`` observations in the two groups 
respecting the original group numerosity. 

*Number of permutations for exact tests:* there are ``\\frac{N!}{N_1 \\cdot N_2}`` possible reassigments 
of the ``N`` observations in the two groups.

*Aliases:* `tTestIS`, [`pointBiSerialTest`](@ref)

*Multiple comparisons version:* [`studentMcTestIS`](@ref)

Return a [UniTest](@ref) structure.

**METHOD (2)**

As (1), but `yvec` is a vector of 2-vectors of observations for group 1 and group 2 (see examples below).

**METHOD (3)**

As (1), but the observations are given separatedly for the two groups as two vectors `y1` and `y2` (see examples below).

*Examples*

```julia
# (1)
using PermutationTests
ns=[4, 5]; # number of observations in group 1 and group 2 (N1 and N2)
y=[randn(n) for n∈ns]; # some Gaussian data as example
t = tTestIS(vcat(y...), ns) # by default the test is bi-directional

# with a bi-directional test, t is equivalent to a 1-way ANOVA for independent samples
tanova= fTestIS(vcat(y...), ns) 
println(t.p ≈ tanova.p ? "OK" : "error")

# do not run it using the CrossProd test statistic
tcor = tTestIS(vcat(y...), ns; asPearson=false) 

# Force an approximate test with 10000 random permutations
tapprox = fTestIS(vcat(y...), ns; switch2rand=1, nperm=10000) 

tR=tTestIS(vcat(y...), ns; direction=Right()) # right-directional test
tL=tTestIS(vcat(y...), ns; direction=Left()) # left-directional test

# in method (2) only the way the input data is formatted is different 
t2 = tTestIS(y)
println(t.p ≈ t2.p ? "OK" : "error")

# in method (3) also, only the way the input data is formatted is different 
t3 = tTestIS(y[1], y[2])
println(t.p ≈ t3.p ? "OK" : "error")

```

**Similar tests**

Typically, the input data is real, but can also be of type integer or boolean. 

For dicothomous data, with this function one can obtain the same p-value as the one given by the 
Fisher exact test, however in this case it is more convenient to use the [`fisherExactTest`](@ref) function, 
since it accepts contingency tables as input. 

This function can also be used to perform a permutation-based
point-biserial correlation test. See the dedicated function [`pointBiSerialTest`](@ref). 

"""
function studentTestIS(𝐱::UniData, 𝐲::UniData;
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            asPearson::Bool = true) where TestDir <: TestDirection
    # no check is needed when data is input this way
    ns=[length(𝐱), length(𝐲)]
    if asPearson # run t-test as a correlation test with reversed membership vector
        𝐱_ = membership(StudentT_IS(), ns; rev=reverse)
        return correlationTest(𝐱_, [𝐱; 𝐲];
            direction, equivalent, switch2rand, nperm, seed, standardized=false, centered=false, verbose)
    else
        return _permutationTest([𝐱; 𝐲], ns, StudentT_IS(), direction, equivalent, switch2rand, nperm, seed, verbose)

    end
end

# alias
tTestIS = studentTestIS


"""
```julia
function pointBiSerialTest(<same args and kwargs as `studentTestIS`>)
```

Actually an alias for [`studentTestIS`](@ref).

Run a [point bi-serial correlation test](https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient) 
between an input vector `y` of ``N=N_1+N_2`` elements and a vector ``x``, internally created, with the first ``N_1`` elements equal to `1`
and the remaining ``N_2`` elements equal to `2`.
If you need to use other values for the dicothomous variable ``x`` or a different order for its elements, 
use [`correlationTest`](@ref) instead. 

The null hypothesis has form 

``H_0: b_{(x,y)}=0``,

where ``b_{(x,y)}`` is the point bi-serial correlation between input data vectors `y` and the internally created vector ``x``.

Directional tests, permutation scheme and number of permutations for exact tests: as per [`studentTestIS`](@ref)

*Multiple comparisons version:* [`pointBiSerialMcTest`](@ref)

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
ns=[4, 6] # number of observations in group 1 and group 2 (N1 and N2)
N=sum(ns) # total number of observations

y = rand(N) # some Gaussian data as example
# implicitly, the point bi serial correlation is 
# between y and x=[1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
t=pointBiSerialTest(y, ns) # by default the test is bi-directional

tR=pointBiSerialTest(y, ns; direction=Right()) # right-directional test
tL=pointBiSerialTest(y, ns; direction=Left()) # left-directional test
```
"""
pointBiSerialTest = studentTestIS # alias


### Chi-Squared and Fisher Exact Test
# See `table2vec` for explanation of the table format.
"""
```julia
function chiSquaredTest(table::Matrix{I};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            asPearson::Bool = true) where {I <: Int, TestDir <: TestDirection}

```

Univariate [chi-squared](https://en.wikipedia.org/wiki/Chi-squared_test) (``\\chi^2``) permutation test 
for ``2 \\cdot K`` contingency tables, where ``K`` is ≥2. 
The null hypothesis has form 

``H_0: O=E``,

where ``O`` and ``E`` are the observed and expected frequencies of the contingency table.

`table` is a contingency table given in the form of a matrix of integers. For example, the contingency table


| 0 | 2 | 3 | Failures 

| 3 | 1 | 0 | Successes


will be given as

```julia
table=[0 2 3; 3 1 0]
```

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests").

For ``K=2`` this function calls [`studentTestIS`](@ref) and pass to it also argument `asPearson`,
otherwise calls [`anovaTestIS`](@ref) and argument `asPearson` is ignored.

In contrast to Pearson's asymptotic ``\\chi^2``, with permutation tests the sample size does not have to be large.
Actually, for large sample sizes Pearson's test is more efficient. For small sample sizes the p-value can be obtained 
using all possible permutations, thus being exact and will be the same as the p-value obtained using 
the [Fisher exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test).

This permutation test is therefore not particularly useful when compared to the standard ``\\chi^2`` and Fisher exact test, 
however, its multiple comparison version allows the control of the family-wise error rate through data permutation. 

*Directional tests*

Possible only for ``𝐾=2``, in which case the test reduces to a Fisher exact test 
and the test directionality is given by keyword arguement `direction`. 
See function [`fisherExactTest`](@ref) and its multiple comparisons version [`fisherExactMcTest`](@ref). 

*Permutation scheme:* the contingency table is converted to ``K`` vectors holding each as many observations as the corresponding column sum. 
The conversion is operated internally by function [`table2vec`](@ref).
The elements of the vectors are as many zeros and ones as the counts of the two cells of the correspondind column.
The F-statistic of the 1-way ANOVA for indepedent samples is then an equivalent test-statistic for the 
``\\chi^2`` and the permutation scheme of that ANOVA applies (see [`anovaTestIS`](@ref)).

*Number of permutations for exact tests:* there are ``\\frac{N!}{N_1 \\cdot\\ldots\\cdot N_K}`` possible permutations, where ``K`` is the number of columns
in the contingency table and ``N_k`` is the ``k^{th}`` column sum.


*Aliases:* `Χ²Test`, [`fisherExactTest`](@ref)

*Multiple comparisons version:* [chiSquaredMcTest](@ref)

Return a [UniTest](@ref) structure.

*Examples*
```julia
using PermutationTests
table=[0 2 2; 3 1 0]
t=Χ²Test(table) # the test is bi-directional
```
```julia
table=[6 1; 2 5]
tR=fisherExactTest(table; direction=Right()) 
# or tR=Χ²Test(table; direction=Right())
```
"""
function chiSquaredTest(table::Matrix{I};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true,
            asPearson::Bool = true) where {I <: Int, TestDir <: TestDirection}
            
    𝐲, ns = table2vec(table, AnovaF_IS()) # anovaTestIS will switch to t-test if there are only two groups
    k = length(ns)
    k < 2 && throw(ArgumentError(📌*"Function chiSquaredTest or fisherExactTest: error reading first argument (table)"))
    !(direction isa Both) && (k > 2) && throw(ArgumentError(📌*"Function chiSquaredTest or fisherExactTest: For input data matrices with more than 2 columns the test can only be bi-directional. Correct or eliminate the `direction` keyword argument")) 
    k==2 ?  tTestIS(𝐲, ns; direction, equivalent, switch2rand, nperm, seed, verbose, asPearson) :
            fTestIS(𝐲, ns; direction, equivalent, switch2rand, nperm, seed, verbose)
end

# alias
Χ²Test = chiSquaredTest

"""
```julia
function fisherExactTest(<same args and kwargs as `chiSquaredTest`>)
```
Alias for [chiSquaredTest](@ref). It can be used for ``2 \\cdot 2`` contingency tables. 
The contingency table in this case has form:

| a | b |

| c | d |

 - For a right-directional test, ``a/c`` is expected to exceed ``b/d``. If the opposite is true, the test will result in a p-value higehr then 0.5.

 - For a left-directional test, ``b/d`` is expected to exceed ``a/c``. If the opposite is true, the test will result in a p-value higehr then 0.5.

!!! tip "Nota Bene" 
    For ``K=2``, any input data matrix gives the same p-value as its transpose.

*Multiple comparisons version:* [fisherExactMcTest](@ref)

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
table=[6 1; 2 5]
t=fisherExactTestTest(table) 
# or t=Χ²Test(table) # bi-directional test
```
```julia
tR=fisherExactTest(table; direction=Right())  # right-directional test
```
"""
fisherExactTest = chiSquaredTest # alias


### ANOVA for Repeated Measures
"""
```julia
# METHOD (1)
function anovaTestRM(y::UniData, ns::@NamedTuple{n::Int, k::Int};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection
```
"""
function anovaTestRM(𝐲::UniData, ns::@NamedTuple{n::Int, k::Int};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection

    ns.k < 2 && throw(ArgumentError(📌*"Function anovaTestRM: the `k` element of second argument (a named tuple) must be 2 or more"))
    ns.n < ns.k && throw(ArgumentError(📌*"Function anovaTestRM: the `n` element of second argument (a named tuple) must be greater than the `k` element"))
    ns.n*ns.k ≠ length(𝐲) && throw(ArgumentError(📌*"Function anovaTestRM: the `n` and `k` elements of second argument (a named tuple) must be such that their product is equal to the length of the first argument (𝐲)"))
    if ns.k == 2
        stat = StudentT_1S() # do t-test if there are only two groups. Actually do one-sample test on the difference
        return _permutationTest!(𝐲[1:2:ns.n*2-1].-𝐲[2:2:ns.n*2], ns.n, direction, equivalent, switch2rand, nperm, seed, verbose)
    else
        stat = AnovaF_RM() 
        !(direction isa Both) && throw(ArgumentError(📌*"Function anovaTestRM: The ANOVA test can only be di-directional. Correct the `direction` keyword argument")) 
    return _permutationTest(𝐲, ns, stat, direction, equivalent, switch2rand, nperm, seed, verbose)
    end
end


"""
```julia
# METHOD (2)
function anovaTestRM(yvec::UniDataVec; <same kwargs>)
```
METHOD (1)

Univariate [1-way analysis of variance (ANOVA) for repeated measures](https://en.wikipedia.org/wiki/Repeated_measures_design#Repeated_measures_ANOVA) 
by data permutation. Given ``N`` observation units (e.g., subjects, blocks, etc.) for each of ``K`` repeated measures 
(e.g., treatments, time, etc.), the null hypothesis has form

``H_0: μ_1= \\ldots =μ_K``,

where ``μ_k`` is the mean of the ``k^{th}`` treatment. 

`y` is a vector concatenaning the observations for the ``K`` treatments in the natural order, that is,
the ``N`` observation for treatment 1, ..., the ``N`` observations for treatment ``K``.
Thus, `y` holds ``N \\cdot K`` elements. 

`ns` is a julia [named tuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) 
with form `(n=N, k=K)` (see examples below).

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests"). 

*Directional tests*

Possible only for ``𝐾=2``, in which case the test reduces to a one-sample Student' t-test on the 
differences of the two treatments and the test directionality is given by keyword arguement `direction`. 
See function [`studentTest1S`](@ref) and its multiple comparisons version [`studentMcTest1S`](@ref). 

*Permutation scheme:* under the null hypothesis, the order of the ``K`` measurements bears no meaning.
The exchangeability scheme consists then in reordering the ``K`` measurements within the ``N`` observation units. 

*Number of permutations for exact tests:* there are ``K!^N`` possible ways of reordering the ``K`` measurements in all ``N`` observation units.

**Alias:** `fTestRM`

*Multiple comparisons version:* [anovaMcTestRM](@ref)

Both methods return a [UniTest](@ref) structure.


METHOD (2)

As (1), but `yvec` is a vector of K vectors holding the observations for the ``k^{th}`` 
measurement (see examples below).

*Examples*

```julia
# (1)
using PermutationTests
N=6; # number of observation units
K=3; # number of measurements
y = [randn(N) for k=1:K] # some random Gaussian data for example 
t = fTestRM(vcat(y...), (n=N, k=K)) # ANOVA tests are always bi-directional
```

```julia
# Force an approximate test with 5000 random permutations
tapprox = fTestRM(vcat(y...), (n=N, k=K); switch2rand=1, nperm=5000)
```

```julia
# in method (2) only the way the input data is formatted is different 
t2 = fTestRM(y)
println(t.p ≈ t2.p ? "OK" : "error")
```

**Similar tests**

Typically, the input data is real, but can also be of type integer or boolean. For dicothomous data, 
with this function one can obtain the permutation-based Cochran Q test for ``K>2`` and the permutation-based McNemar test 
for ``K=2``, but with the ability to give exact p-values.
For these two tests the dedicated functions [`cochranqTest`](@ref) and [`mcNemarTest`](@ref) is available.
"""
function anovaTestRM(𝐲vec::UniDataVec;
                direction::TestDir = Both(),
                equivalent::Bool = true,
                switch2rand::Int = Int(1e8),
                nperm::Int = 20_000, 
                seed::Int = 1234, 
                verbose::Bool = true) where TestDir <: TestDirection

    K = length(𝐲vec)
    K < 2 && throw(ArgumentError(📌*"Function anovaTestRM: the first argument (𝐲vec) must be a vector of two or more vectors"))
    nn = length(unique(length(𝐲vec)))
    nn ≠ 1 && throw(ArgumentError(📌*"Function anovaTestRM: For a repeated-measure ANOVA all vectors in argument `𝐲vec` must be of the same length"))
    N = length(𝐲vec[1])
    N < K && throw(ArgumentError(📌*"Function anovaTestRM: the length of the vectors in first argument (𝐲vec) must be greater than their number"))
    if K == 2
        stat = StudentT_1S() # do t-test if there are only two groups. Actually do one-sample test on the difference
        return _permutationTest!(𝐲vec[1].-𝐲vec[2], N, direction, equivalent, switch2rand, nperm, seed, verbose)
    else
        stat = AnovaF_RM() 
        !(direction isa Both) && throw(ArgumentError(📌*"Function anovaTestRM: The ANOVA test can only be bi-directional. Correct the `direction keyword argument`")) 
        ns=(n=N, k=K) 
        return _permutationTest(vcat(𝐲vec...), ns, stat, direction, equivalent, switch2rand, nperm, seed, verbose)
    end
end

fTestRM = anovaTestRM # alias

### Cocharn Q test
# NB: not particular efficient. For an efficient alternative see https://dl.acm.org/doi/10.1145/3095076
"""
```julia
function cochranqTest(table::Matrix{I};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where {I <: Int, TestDir <: TestDirection}

```

Univariate [Cochran Q](https://en.wikipedia.org/wiki/Cochran's_Q_test) test by data permutation. The Cochran Q test 
is analogous to the 1-way ANOVA for repeated measures, but takes as input dicothomous data (zeros and ones).
Given ``N`` observation units (e.g., subjects, blocks, etc.) and ``K`` repeated measures 
(e.g., treatments, time, etc.), the null hypothesis has form

``H_0: μ_1= \\ldots =μ_K``,

where ``μ_k`` is the mean of the ``k^{th}`` measure. 

When ``K=2`` the test reduces to the [McNemar test](https://en.wikipedia.org/wiki/McNemar's_test),
which is the analogous to the Student's t-test for repeated measures taking as input dicothomous data (zeros and ones). 

Input `table` is a ``N \\cdot K`` table of zeros and ones, where ``N`` is the number of observations and 
``K`` the repeated measures. Transposed, one such data would look like


| 1 | 1 | 1 | 1 | 1 | 1 | Measure 1

| 1 | 0 | 1 | 1 | 0 | 1 | Measure 2

| 0 | 0 | 1 | 0 | 1 | 0 | Measure 3

This table shall be given as input such as 
    ```table=[1 1 0; 1 0 0; 1 1 1; 1 1 0; 1 0 1; 1 1 0]```

and internally it will be converted to the appropriate format by function [`table2vec`](@ref).

!!! note "Redundant permutations"
    Adding any number of vectors `[0 0 0]` or `[1 1 1]` in any combination to the table here above, 
    yields exactly the same p-value using systematic permutations, but increase the number of permutations
    to be listed. If such vector exist in your data, you can delete them to obatin a faster test.

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests").

In contrast to the [Cochran Q](https://en.wikipedia.org/wiki/Cochran's_Q_test) and 
[McNemar](https://en.wikipedia.org/wiki/McNemar's_test) asymptotic tests, 
with permutation tests the sample size does not have to be large.
Actually, for large sample sizes the Cochran Q and McNemar tests are more efficient, although they do not provide 
an exact p-value (an exact test for the case ``K=2`` can be derived though).
Overall, this permutation test is not very useful when compared to the standard Cochran Q and McNemar tests, 
however, its multiple comparison version allows the control of the family-wise error rate by means of data permutation. 

*Directional tests*

Possible only for ``𝐾=2``, in which case the test reduces to a MnNemar test 
and the test directionality is given by keyword arguement `direction`. 
See function [`mcNemarTest`](@ref) and its multiple comparisons version [`mcNemarMcTest`](@ref). 

*Permutation scheme:* the F-statistic of the 1-way ANOVA for repeated measure is an equivalent test-statistic for the 
Cochran Q test and the permutation scheme of that ANOVA applies (see [`anovaTestRM`](@ref)).
The permutation scheme is the same for tha case of ``K=2`` (McNemar test).

*Number of permutations for exact tests:* there are ``K!^N`` possible ways of reordering the ``K`` measurements in all the ``N`` observation units.

*Aliases:* `qTest`, `mcNemarTest`

*Multiple comparisons versions:* [`cochranqMcTest`](@ref), [`mcNemarMcTest`](@ref) 

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
table=[1 1 0; 1 0 0; 1 1 1; 1 1 0; 1 0 1; 1 1 0]
t=qTest(table) # the test is bi-directional
```
"""
function cochranqTest(table::Matrix{I};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where {I <: Int, TestDir <: TestDirection}
    𝐲, ns = table2vec(table, AnovaF_RM()) # anovaTestRM will switch to t-test if there are only two groups
    ns.k == 2 ? studentTest1S!(𝐲[1:2:(ns.n)*2-1].-𝐲[2:2:(ns.n)*2]; refmean = nothing, direction, equivalent, switch2rand, nperm, seed, verbose) :
                anovaTestRM(𝐲, ns; direction, equivalent, switch2rand, nperm, seed, verbose)
end

qTest = cochranqTest # alias


"""
```julia
function mcNemarTest(same args and kwargs as `cochranqTest`>)
```
Alias for [`cochranqTest`](@ref). It can be used for ``2 \\cdot 2`` contingency tables.

Notice that [`cochranqTest`](@ref) does not accept data input in the 
form of a contingency table.
If your data is in the form of a contingency table, here is how you can convert it:

given the contingency table

| a | b |

| c | d |

you will create a vector holding:

 - as many vectors `[0, 1]` as the ``b`` frequency.
 - as many vectors `[1, 0]` as the ``c`` frequency.
 
 For example, the contingency table

| 1 | 2 |

| 3 | 4 |

will be given as input such as as

```julia
table=[0 1; 0 1; 1 0; 1 0; 1 0]
```

!!! note "Redundant permutations"
    Adding any number of vectors `[0 0]` or `[1 1]` in any combination to the table here above, 
    yields exactly the same p-value using systematic permutations.
    Like in the asymptotic McNemar test, these vector correspond to elements `a` and `d` of the 
    contingency table and have no effect.

*Directional tests*

 - For a right-directional test, ``c`` is expected to exceed ``b``. If the opposite is true, the test will result in a p-value higher then 0.5.
 - For a left-directional test, ``b`` is expected to exceed ``c``. If the opposite is true, the test will result in a p-value higher then 0.5.

*Multiple comparisons version:* [`mcNemarMcTest`](@ref)

Return a [UniTest](@ref) structure.

*Examples*
```julia
using PermutationTests
table=[1 0; 1 0; 1 0; 1 0; 0 1]
t=mcNemarTest(table) # by default the test is bi-directional

tR=mcNemarTest(table; direction=Right()) # right-directional test
```
"""
mcNemarTest = cochranqTest # alias

# https://en.wikipedia.org/wiki/McNemar%27s_test


### One-Sample t-test
"""
```julia
function studentTest1S(𝐲::UniData;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection =
```

Univariate [one-sample t-test](https://en.wikipedia.org/wiki/Student's_t-test#One-sample_t-test) 
by data permutation. The null hypothesis has form 

``H_0: μ=μ_0``,

where ``μ`` is the mean of the observations and ``μ_0`` is a reference population mean.

`refmean` is the reference mean (``μ_0``) above. The default is ``0``, which is the value used for most tests.

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests").

*Directional tests*

 - For a right-directional test, ``μ`` is expected to exceeds ``μ_0``. If the opposite is true the test will result in a p-value higehr then 0.5.
 - For a left-directional test, ``μ_0`` is expected to exceeds ``μ``. If the opposite is true the test will result in a p-value higehr then 0.5.

*Permutation scheme:* the one-sample t-test test is equivalent to a t-test for two repeated measures, the mean of the difference of which 
is tested by the one-sample t-test. Under the null hypothesis, the order of the two measurements bears no meaning.
The exchangeability scheme consists then in reordering the two measurements in the ``N`` observation units. 
For a one-sample t-test, this amounts to considering the observations with the observed and switched sign.

*Number of permutations for exact tests:* there are ``2^N`` possible ways of reordering the two measurements in the ``N`` observations of `y`.
For the one-sample t-test, this is the number of possible configurations of the signs of the observations.

**Alias:** `tTest1S` 

Multiple comparisons version: [`studentMcTest1S`](@ref)

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
N=20 # number of observations
y = randn(N) # some random Gaussian data for example
t = tTest1S(y) # test H0: μ(y)=0. By deafult the test is bi-directional

t = tTest1S(y; refmean=1.) # test H0: μ(y)=1
tR = tTest1S(y; direction=Right()) # right-directional test
tL = tTest1S(y; direction=Left()) # Left-directional test
# Force an approximate test with 5000 random permutations
tapprox = tTest1S(y; switch2rand=1, nperm=5000) 
```

**Similar tests**

Typically, the input data is real, but can also be of type integer or boolean. If either `y` is a vector 
of booleans or a vector of dicothomous data (only 0 and 1), this function will actually perform a permutation-based 
version of the [sign test](https://en.wikipedia.org/wiki/Sign_test).
With boolean input, use the [signTest](@ref) function.

Passing as data input the vector of differences of two repeated measurements, this function carries out the 
Student's t-test for repeated measurements. If you need such a test you may want to use
the [`studentTestRM`](@ref) alias.
"""
studentTest1S(𝐲::UniData;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection =
    _permutationTest!(refmean===nothing ? copy(𝐲) : 𝐲.-refmean, length(𝐲), direction, equivalent, switch2rand, nperm, seed, verbose)

tTest1S = studentTest1S # aliases    


"""
```julia
function studentTest1S!(y::UniData; <same args and kwargs as `studentTest1S`>)
```

As [`studentTest1S`](@ref), but `y` is overwritten in the case of approximate (random permutations) tests.

**Alias:** `tTest1S!`

**Multiple Comparison version:** [`studentMcTest1S!`](@ref)

"""
studentTest1S!(𝐲::UniData;
            refmean::Realo = nothing,
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection =
    _permutationTest!(refmean===nothing ? 𝐲 : 𝐲.-refmean, length(𝐲), direction, equivalent, switch2rand, nperm, seed, verbose)

# alias
tTest1S! = studentTest1S!
# as studentTest1S!, but copy 𝐲, so it does not overwrite it. This method directly subtracts the reference mean    


### Sign Test
"""
```julia
signTest(𝐲::Union{BitVector, Vector{Bool}};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection =

```

**Multiple comparisons version**

[`signMcTest`](@ref)

Univariate [sign test](https://en.wikipedia.org/wiki/Sign_test) by data permutation. The null hypothesis has form 

``H_0: E(true)=E(false)``,

where ``E(true)`` and ``E(false)`` are the expected number of true and false occurrences, respectively.

`y` ia a vector of ``N`` booleans.

For optional keyword arguments, `direction`, `equivalent`, `switch2rand`, `nperm`, `seed` and `verbose`, 
see [here](@ref "Common kwargs for univariate tests").

*Directional tests*

- For a right-directional test, ``E(true)`` is expected to exceeds ``E(false)``. If the opposite is true the test will result in a p-value higehr then 0.5.
- For a left-directional test, ``E(false)`` is expected to exceeds ``E(true)``. If the opposite is true the test will result in a p-value higehr then 0.5.

Permutation scheme and number of permutations for exact tests: as per [`studentTest1S`](@ref). 

The significance of the univariate sign test can be obtained much more efficiently using the binomial distribution.
This permutation test is therefore not useful at all in the univariate case, 
however, its multiple comparison version allows the control of the family-wise error rate by data permutations.

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
N=20; # number of observations
y = rand(Bool, N); # some random Gaussian data for example
t = signTest(y) # By deafult the test is bi-directional

tR = signTest(y; direction=Right()) # right-directional test
tL = signTest(y; direction=Left()) # Left-directional test
# Force an approximate test with 5000 random permutations
tapprox = signTest(y; switch2rand=1, nperm=5000) 
```
"""
signTest(𝐲::Union{BitVector, Vector{Bool}};
            direction::TestDir = Both(),
            equivalent::Bool = true,
            switch2rand::Int = Int(1e8),
            nperm::Int = 20_000, 
            seed::Int = 1234, 
            verbose::Bool = true) where TestDir <: TestDirection =
    studentTest1S!((Float64.(𝐲)).-0.5; refmean=nothing, direction, equivalent, switch2rand, nperm, seed, verbose)


"""
```julia
function studentTestRM(<same args and kwargs as `studentTest1S`>)
```

Actually an alias for [`studentTest1S`](@ref).
   
In order to run a t-test for repeated measure, use as data input the vector of differences across measurements.

Do not change the `refmean` default value. See [`studentTest1S`](@ref) for more details.

**Alias:** `tTestRM`

*Multiple comparisons version:* [`studentMcTestRM`](@ref)

Return a [UniTest](@ref) structure.

*Examples*

```julia
using PermutationTests
y1=randn(10) # measurement 1
y2=randn(10) # measurement 2
t=tTestRM(y1.-y2) #  # by default the test is bi-directional

tR=tTestRM(y1.-y2; direction=Both()) #  right-directional test
# if test tR is significant, the mean of measurement 1 exceeds the mean of measurement 2.
```
"""
studentTestRM = studentTest1S

tTestRM = studentTest1S

"""

```julia
function studentTestRM!(<same args and kwargs as `studentTestRM`>)
```    
Actually an alias for [`studentTest1S!`](@ref).

See [`studentTestRM`](@ref) for the usage of this function.

**Alias:** `tTestRM!`

*Multiple comparisons version:* [`studentMcTestRM!`](@ref)
"""
studentTestRM! = studentTest1S!

tTestRM! = studentTest1S!

# ----------------------------------------------------- #

