# Create your own test

This page illustrates how you can create *univariate* and *multiple comparisons* permutation tests
besides those already supported by *PermutationsTests.jl*.

Good knowledge and understanding of the whole documentation is required for continuing reading this page,
including the **Extras** pages.

There are two ways you can extend the capability of *PermutationTests.jl*.

 - The first, very simple, but limited, is by using existing tests with particular data input.

 - The second, which needs some development, but is much more general, allows you to create new tests defining your own test statistics.

## Using existing tests

The approach is illustrated by creating tests for the autocorrelation, which are not in the API of the package.

**Example of a univariate permutation test**

As an example, let us detect the presence of autocorrelation at lag 1, similar to what the 
[Durbin-Watson test](https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic) does.

*PermutationTests.jl* implements the Pearson product moment correlation statistic and a test for it.
We can then test the null hypothesis that the autocorrelation at lag 1 is equal to zero
by performing a correlation test on a data vector `x` and a lagged version of itself `y`.

This is simply achieved by calling the [`_permTest!`](@ref), which performs all 
[univariate tests](@ref "Univariate tests") implemented in the package:

```@example
using PermutationTests
# First, let us create some data autorrelated at lag 1 as example
N=300
data=randn(N) 
for i=2:length(data)
  data[i-1]=data[i-1]+data[i]
end

# Data and lagged version of the data
x=data[1:end-1]
y=data[2:end]

# Our Univariate test on autocorrelation
t=_permTest!(copy(x), y, N-1, PearsonR(), PearsonR())
# NB: we pass a copy of x as this vector will be permuted. 
```

The testing procedure can be further improved and refined, see for example [`correlationTest`](@ref),
but this should suffice to give the idea.

---

**Example of a multiple comparison permutation test**

Similarly, we can easily test simultaneously the autocorrelation at several lags with a 
multiple comparison version of the above test.

This is achieved by calling the [`_permMcTest!`](@ref), which performs all 
[multiple comparisons tests](@ref "Multiple comparisons tests") implemented in the package.
Now we craete a vector `x` and a vector `Y` holding as many lagged versions of `x` as we wish to test:

```@example
using PermutationTests
# Let us create some data autorrelated at lag 1 and 2 as example
N=300
data=randn(N) 
for j=1:2, i=2:length(data)
  data[i-1]=data[i-1]+data[i]
end

# We will test lags 1 to 30
lags=1:30 
x=data[1:end-length(lags)]
Y=[data[l+1:end-length(lags)+l] for l in lags]

# Our muliple comparisons test on autocorrelations
tMc=_permMcTest!(copy(x), Y, N-length(lags), PearsonR(), PearsonR())
```

The p-values are retrived by 

```julia
tMc.p # p-values for lags 1...30 
```

Notice that the test above is different from the 
[Liung-Box test](https://en.wikipedia.org/wiki/Ljung%E2%80%93Box_test),
which is a multivariate test for a group of lags:
with a multiple comparison tests we have tested each lag
in the group *simultaneously*, controlling the family-wise error (FWE) rate at the nominal level
(0.05 by default). Thus on the base of the p-values 
*we can take a decision on the null hypothesis at each lag*.

You can check yourself that the test controls the FWE rate at the nominal level
by running the test above with random Gaussian data many times and controlling 
that the proportion of rejections, which are false by definition since random Gaussian data
is not autocorrelated, does not exceed the nominal level.

---

## Defining new test statistics

You can create your own permutation tests based on a *custom test statistic* as long as the permutation scheme 
for your test is supported by *PermutationTests.jl*. Four permutation schemes are implemented in the package, corresponding to four implemented groups of test-statistics.
Before continuing, read the documentation of [`genPerms`](@ref) and the documentation linked therein.

In order to create a test you will call function [`_permTest!`](@ref)(univariate) 
or [`_permMcTest!`](@ref)(multiple comaprisons). You will pass, among others, argument `asStat`,
which determines the permutation scheme of the test.

!!! tip "keep in mind"
    The `asStat` argument can be one of the following built-in test statistics or any test statistic belonging to the same [group](@ref "Statistic groups"):

    - `PearsonR()`
    - `AnovaF_IS()`
    - `AnovaF_RM()`
    - `StudentT_1S()`

The procedure to follow is slightly different depending on whether you want 
to create your own univariate or multiple comparisons test.

---

#### UNIVARIATE TESTS

Creating your own univariate test involves two preliminary steps:

1) define a new [`Statistic`](@ref) type to name the test statistic you want to use
2) write a method for the [`statistic`](@ref) function to compute the test statistic for both the observed and permuted data

You then obtain the test just calling the [`_permTest!`](@ref) functions, as we will show.

!!! tip "keep in mind"
    Argument `fstat` of `_permTest!` by default is the `abs` julia function, which yields a bi-directional test. 
    If appropriate for the test statistic at hand, use `fstat=identity` for a right-directional test and 
    `fstat=`[`flip`](@ref) for a left-directional test. 

---

**Example of a univariate permutation test**

We create a new test for the difference of the mean between two independent samples.
Such a test is equivalent to the *Student's t-test for independent samples*, since the mean difference
is the nominator of the t statistic and the denominator is invariant with respect to data permutation.

The permutations scheme is the same as the permutation scheme of the `AnovaF_IS()` (or `StudentT_IS`) 
test statistic, hence it is supported by *PermutationTests.jl* and we can reuse it.

```@example
# First, let us import `statistic` from PermutationTests.jl
using PermutationTests
import PermutationTests: statistic

# 1) define a new `Statistic` type to name the test statistic you want to use
struct MeanDiff <: Statistic end

#=
2) write a method for the `statistic` function to compute the mean difference, 
it does not matter if the data is observed or permuted. 
We must input the data exactly as for the `AnovaF_IS()` statistic. 

Our observations are divided in two groups. In this case,
ror this and equivalent statistics, data input is given as a single vector 
y of length N1+N2, with the N1 observations followed by the N2 observations. 

The `statistic` method also takes as input a 'membership' vector x, 
holding 1 for subjects of group 1 and 2 for subjects of group 2,
that is, x=[repeat ([1], N1); repeat([2], N2)].
When given as input, x will reflect the positioning of the observed data 
and for permuted data it will be a permutation thereof. 

A simple implementation: 
=#
function statistic(x, y, stat::MeanDiff; kwargs...)
    # get the number of subjects in each gruoup
    ns=[count(j->j==i, x) for i=1:length(unique(x))]

    s=[0., 0.]
    @simd for i ∈ eachindex(x, y) 
        @inbounds (s[x[i]] += y[i])
    end
    return (s[1]/ns[1]) - (s[2]/ns[2]) # mean difference
end

#=
The following code block is not necessary. 
We report it to illustrate a nice opportunity you have to write 
efficient code to compute test statistics;
the above implementation computes the group numerosity (ns=[N1, N2]) 
at each permutation, however ns is invariant to permutations. 

Let us rewrite the function then using the `cpcd` argument, 
which will be passed to the `_permTest!` functions in order to 
(optionally) pre-compute `ns`.
=#
function statistic(x, y, stat::MeanDiff; cpcd=nothing, kwargs...)
    if cpcd===nothing 
        ns=[count(j->j==i, x) for i=1:length(unique(x))]
    else
        ns=cpcd
    end

    s=[0., 0.]
    @simd for i ∈ eachindex(x, y) 
        @inbounds (s[x[i]] += y[i])
    end
    return (s[1]/ns[1]) - (s[2]/ns[2]) # mean difference
end

# That's it. Let's test our own test:

# create some data
y=[4., 5., 6., 7., 1., 2., 3]
ns=[4, 3] # N1=4, N2=3

# membership vector (the dedicated function `membership` is available)
x=membership(StudentT_IS(), ns) # = [repeat ([1], N1); repeat([2], N2)]

println("t-test independent samples (bi-directional) using PermutationsTest.jl: ");
# two-directional exact test. Equivalent to API function studentTestIS(y, ns)
t1 = _permTest!(copy(x), y, ns, StudentT_IS(), StudentT_IS()) 
# Note: we pass a copy of x as this vector will be permuted. 

# This runs the test we have created. Note that argument `stat` of `_permTest!`
# is `MeanDiff()`, while argument `asStat` is `AnovaF_IS()` 
println("t-test independent samples (bi-directional) using our own statistic: ");
t2 = _permTest!(copy(x), y, ns, MeanDiff(), AnovaF_IS())
```

As before, but using the `cpcd` argument to avoid redundant computations:

```julia
t2_ = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); cpcd=ns) 
```

You can check that the p-values obtained by the two tests is identical:
```julia
t1.p≈t2.p ? (@info "my test works!") : 
            (@error "my test does not work!" t1.p t2.p)
```

To run a right-directional test:
```julia
t3 = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); fstat = identity)
```

To run a left-directional approximate test (see [`_permTest!`](@ref)):
```julia
t4 = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); fstat = flip, switch2rand=1)
```

etc.

---

#### MULTIPLE COMPARISON TESTS

Creating your own multiple comparison test involves three preliminary steps:

1) define a new [`Statistic`](@ref) type to name the test statistic you want to use
2) write a method of the [`_observedStats`](@ref) function to compute function `fstat` of all observed test-statistics at once as a vector
2) write a method of the [`_permutedStat`](@ref) function to compute function `fstat` of each permutation statistics separately 

You then obtain the test just calling the [`_permMcTest!`](@ref) functions.

!!! tip "keep in mind"
    Argument `fstat` of `_permMcTest!` by default is the `abs` julia function, which yields a bi-directional test. 
    If appropriate for the test statistic at hand, use `fstat=identity` for a right-directional test and 
    `fstat=`[`flip`](@ref) for a left-directional test. 

---

**Example of a multiple comparison permutation test**

We create another version of the Pearson product-moment correlation test between a fixed variable given as a vector `x` and 
``M`` variables given as a vector of ``M`` vectors `Y`. We want to test simultaneously 
all correlations between `x` and `Y[m]`, for ``m=1...M``.

The permutations scheme is the one of the `PearsonR()` statistic, hence it is supported by *PermutationTests.jl*
and we can reuse it (as `asStat` argument of `_permMcTest!`).

```@example
# First, let us import some names from PermutationTests.jl
using PermutationTests
import PermutationTests: statistic, _observedStats, _permutedStat

using LinearAlgebra: dot, ⋅ # we are going to use the dot product

# 1) Define a new `Statistic` type to name the test statistic you want to use
struct MyPearsonR <: Statistic end

#= 
2) Write a function to compute function `fstat` of all observed 
test-statistics at once as a vector. 

In order to obtain a faster test we will work with standardized 
variables, thus the Pearson correlation is given simply by the 
cross-product divided by N. 

Note that we will subtract a small constant to all 
obseved statistics (sqrt(eps())) to avoid numerical errors when 
comparing the observed statistics to the permuted statistics. 
=#
_observedStats(x, Y, stat::MyPearsonR, fstat::Function=abs; kwargs...) = 
    [fstat((x ⋅ y)/length(x)) - sqrt(eps()) for y ∈ Y]

#=
3) Write a function to compute function `fstat` of the mth 
permuted test-statistic. 

Again, we will work with standardized variables, thus the Pearson 
correlation is given simply by the cross-product divided by N.
=#
_permutedStat(x, Y, m::Int, stat::MyPearsonR, fstat::Function=abs; kwargs...) = 
    fstat((x ⋅ Y[m])/length(x))


# That's it. Let's test our own test:

# Let's create some data as example
N, M = 7, 100 # N=7 observations, M=100 hypotheses
x=randn(N);
Y=[randn(N) for m=1:M];

println("correlation test using PermutationsTest.jl: ");
# Exact test. Equivalent to API function rMcTest!(x, Y)
# NB: the standardization function μ0σ1() is exported by PermutationTests.jl
t1 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y ∈ Y], length(x), PearsonR(), PearsonR()) 

println("correlation test using our own test: ");
t2 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y ∈ Y], length(x), MyPearsonR(), PearsonR()) 

```

Check that the p-values obtained by the two tests is identical:
```julia
abs(sum(t1.p-t2.p))<1e-8 ?  (@info "my test works!") : 
                            (@error "my test does not work!")
```

Check also the observed statistics:

```julia
abs(sum(t1.obsstat-t2.obsstat))<1e-8 ?  (@info "my test works!") : 
                                        (@error "my test does not work!")
```

To run a right-directional test:
```julia
t3 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y ∈ Y], length(x), MyPearsonR(), PearsonR(); 
        fstat = identity) 
```

To run a left-directional approximate test:
```julia
t4 = _permMcTest!(μ0σ1(x), [μ0σ1(y) for y ∈ Y], length(x), MyPearsonR(), PearsonR(); 
        fstat = identity, switch2rand=1)
```

---

We can also pre-compute data that is invariant to permutations in order to obtain 
a faster procedure. 

For example, suppose we didn't know that the cross-product 
divided by `N` is equal to the Pearson statistic when the data is standardized. 
Instead of computing the Pearson statistics each time from scratch at each permutation, 
we could still pre-compute the means and standard deviations, which are invariant 
to pemutations, using the `cpcd` (custom pre-computed data) keyword argument, such as:

```julia
# Pre-compute means and standard deviations. Here you can compute anything you need
# NB: the mean μ() and st. deviation function σ(y) are exported by PermutationTests.jl 
cpcd=[μ(x), [μ(y) for y in Y], σ(x), [σ(y) for y in Y]]
```

Declare functions for observed and permuted data using argument cpcd. Note that

 - `_observedStats` computes all ``M`` statistics, while `_permutedStat` computes the ``m^{th}`` statistic.
 - not the statistics, but function `fstat` of the statistics are returned

```julia
_observedStats(x, Y, stat::MyPearsonR, fstat::Function=abs; cpcd=cpcd, kwargs...) = 
    [fstat( (((x.-cpcd[1])./cpcd[3]) ⋅ ((y.-cpcd[2][m])./cpcd[4][m])) / length(x) ) - 
        sqrt(eps()) for (m, y) ∈ enumerate(Y)]

_permutedStat(x, Y, m::Int, stat::MyPearsonR, fstat::Function=abs; cpcd=cpcd, kwargs...) = 
    fstat((((x.-cpcd[1])./cpcd[3]) ⋅ ((Y[m].-cpcd[2][m])./cpcd[4][m])) / length(x))
```

```julia
println("correlation test using my own test and pre-computed data: ");
t2_ = _permMcTest!(copy(x), Y, length(x), MyPearsonR(), PearsonR(); cpcd=cpcd) 
```

Note that we copy `x` because otherwise, using this permutation scheme, it is overwritten.
If you don't care, you can simply pass `x`. 

!!! note "nota bene"
    If you reuse the permutations scheme given with  `asStat=StudentT_1S()`, `x` or `Y` are overwritten, depending on whether the test is exact or approximated. 

Note also that now running the `_permMcTest!` function without the `cpcd` kwarg gives an error,
since we did not implement versions of `_observedStats` and of `_permutedStat` 
that take it as an optional argument (although we could, as we have done
for the univariate test). 

---

## Useful functions for creating your own tests

---

```@docs
_permTest!
_permMcTest!
flip
membership
_observedStats
_permutedStat
_∑y²
_∑Y²kn_∑y²_∑S²k
```
