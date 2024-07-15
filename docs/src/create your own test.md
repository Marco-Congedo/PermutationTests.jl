# Create your own test

This page illustrates how you can create **univariate** and **multiple comparisons** permutation tests
besides those already supported by *PermutationsTests.jl*.

Good knowledge and understanding of the whole documentation is required for continuing reading this page,
including the **Extras** pages.

There are two ways you can extend the capability of *PermutationTests.jl*.

 - The first, very simple, but limited, is by using existing tests with particular data input.

 - The second, which needs a little development, but is much more general, allows you to create new tests defining your own test statistics.

The illustration proceeds by examples of increasing complexity.

## Index of examples

| Serial | Examples | 
|:----------:|:----------|
| 1 | [Univariate autocorrelation test](@ref "Example 1: univariate autocorrelation test") |
| 2 | [Multiple comparison autocorrelation test](@ref "Example 2: multiple comparison autocorrelation test") |
| 3 | [Univariate t-test for independent samples](@ref "Example 3: univariate t-test for independent samples") |
| 4 | [Multiple comparison correlation test](@ref "Example 4: multiple comparison correlation test") |
| 5 | [Univariate Chatterjee correlation](@ref "Example 5: univariate Chatterjee correlation") |
| 6 | [Multiple comparisons Chatterjee correlation](@ref "Example 6: multiple comparisons Chatterjee correlation") |
| 7 | [Univariate distance correlation](@ref "Example 7: univariate distance correlation") |
| 8 | [Multiple comparison distance correlation](@ref "Example 8: multiple comparison distance correlation") |


## Using existing tests

The approach for constructing univariate and multiple comparisons tests is illustrated by creating tests for the **autocorrelation**, which are not in the API of the package.

---

### Example 1: univariate autocorrelation test

As an example, let us detect the presence of autocorrelation at lag 1, similar to what the 
[Durbin-Watson test](https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic) does.

*PermutationTests.jl* implements the Pearson product moment correlation statistic and a test for it.
We can then test the null hypothesis that the autocorrelation at lag 1 is equal to zero
by performing a correlation test on a data vector `x` and a lagged version of itself `y`.

This is simply achieved by calling the [`_permTest!`](@ref) function, which ultimately performs all 
[univariate tests](@ref "Univariate tests") implemented in the package.

First, let us create some data autorrelated at lag 1 to get a working example:

```julia
using PermutationTests

N=300
data=randn(N) 
for i=2:length(data)
  data[i-1]=data[i-1]+data[i]
end

# Data and lagged version of the data
x=data[1:end-1]
y=data[2:end]
```

Our univariate test on autocorrelation is obtained with

```julia
t = _permTest!(copy(x), y, N-1, PearsonR(), PearsonR())
```

NB: we pass a copy of `x` as this vector will be permuted. 

If we wish to wrap the test in a friendly function, we can write

```julia
autocorTest(x; kwargs...) = 
  _permTest!(x[1:end-1], x[2:end], length(x)-1, PearsonR(), PearsonR(); kwargs...)
```

The test will be executed then simply as

```julia
t2 = autocorTest(data)
```

By default the test is bi-directional, thus the above call is a test for the autocorrelation regardless its sign.
Our function `autocorTest` accepts all optional keyword arguments of the [`_permTest!`](@ref) function,
thus we can easily obtain exact or approximate tests with a specified number of permutations, left- or right-directional tests, etc. 

---

### Example 2: multiple comparisons autocorrelation test

With the same ease, we can test simultaneously the autocorrelation at several lags using the 
multiple comparison version of the correlation test.

This is achieved by calling the [`_permMcTest!`](@ref), which ultimately performs all 
[multiple comparisons tests](@ref "Multiple comparisons tests") implemented in the package.
Now we craete a vector `x` and a vector `Y` holding as many lagged versions of `x` as we wish to test:

Let us create some data autorrelated at lag 1 and 2 to get a working example:

```julia
using PermutationTests

N=300
data=randn(N) 
for j=1:2, i=2:length(data)
  data[i-1]=data[i-1]+data[i]
end

# We will test lags 1 to 30
lags=1:30 
x=data[1:end-length(lags)]
Y=[data[l+1:end-length(lags)+l] for l in lags]
```

Our muliple comparisons test on autocorrelations is

```julia
tMc = _permMcTest!(copy(x), Y, N-length(lags), PearsonR(), PearsonR())
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

As we have done for the univariate test, we may want to wrap the call in a friendly function.
This is left as an exercise.

---

## Defining new test statistics

You can create your own permutation tests based on a *custom test statistic* as long as  
the permutation scheme for your test is supported by *PermutationTests.jl*. 
The ability to create **a test for virtually whatever test-statistic** is a fine characteristics
of *PermutationTests.jl* and will be illustrated below with examples of ascending complexity.

The procedure to follow is the same whether you want 
to create your own **univariate** or **multiple comparisons test**.
It involves two preliminary steps:

1) define a new [`Statistic`](@ref) type to name the test statistic you want to use
2) write a method for the [`testStatistic`](@ref) function to compute the test statistic for both the observed and permuted data

You then obtain the test just calling the [`_permTest!`](@ref) function for an univariate test
and the [`_permMcTest!`](@ref) for a multiple comparisons test, as we will show.
Two arguments of these functions are particularly important:

- Argument `asStat`, which determines the permutation scheme of the test. 

!!! tip "keep in mind"
    Four permutation schemes are implemented in the package, corresponding to four implemented groups of test-statistics.
    The `asStat` argument can be one of the following built-in test statistics or any test statistic belonging to the same [group](@ref "Statistic groups"):

    - `PearsonR()`
    - `AnovaF_IS()`
    - `AnovaF_RM()`
    - `StudentT_1S()`

    Refer to the documentation of [`genPerms`](@ref) and the documentation linked therein.

- Argument `fstat`, which is a function applied to all observed and permuted statistic.

!!! tip "keep in mind"
    Argument `fstat` of `_permTest!` by default is the `abs` julia function, which yields a bi-directional test for
    a test-statistic that is distributed symmetrically around zero. For such test-statistics you will use 
    `fstat=identity` for a right-directional test and `fstat=`[`flip`](@ref) for a left-directional test. Note that this may not be appropriate for your own test statistic. For example, if your test statistic is not symmetric and a right-tailed test is of interest, you must use `fstat=identity`. See below for examples.


!!! warning "important"
    The permutation vector, which is always the first argument of the [`_permTest!`](@ref) function for an univariate test
    and of the [`_permMcTest!`](@ref) for a multiple comparisons test must correspond to the permutation
    yielding the observed statistic (*i.e.*, no data permutation).

---

### Example 3: univariate t-test for independent samples

We create a new test for the difference of the mean between two independent samples.
Such a test is equivalent to the *Student's t-test for independent samples*, since the mean difference
is the nominator of the t statistic and the denominator is invariant with respect to data permutation.

The permutations scheme is the same as the permutation scheme of the `AnovaF_IS()` (or `StudentT_IS`) 
test statistic, hence it is supported by *PermutationTests.jl* and we can reuse it.

First, let us import `testStatistic` from *PermutationTests.jl* since we will write a new method for it:

```julia
using PermutationTests
import PermutationTests: testStatistic
```
(1) Define a new `Statistic` type to name the test statistic you want to use

```julia
struct MeanDiff <: Statistic end
```

(2) Write a method for the `testStatistic` function to compute the mean difference, it does not matter if the data is observed or permuted. 

We will input the data exactly as for the `AnovaF_IS()` (or `StudentT_IS`) statistic. 
Our observations are divided in two groups. In this case,
for this and equivalent statistics, data input is given as a single vector 
`y` of length ``N1+N2``, with the N1 observations followed by the N2 observations
(see [`studentTestIS`](@ref)). 

The first argument of the `testStatistic` method is **always** the permutation vector x, 
holding in this case 1 for subjects of group 1 and 2 for subjects of group 2,
that is, `x=[repeat ([1], N1); repeat([2], N2)]`.
When given as input, `x` will reflect the positioning of the observed data 
and for permuted data it will be a permutation thereof. 

A simple implementation could be: 

```julia
function testStatistic(x, y, stat::MeanDiff; kwargs...)
    # get the number of subjects in each gruoup
    ns=[count(j->j==i, x) for i=1:length(unique(x))]

    s=[0., 0.]
    @simd for i ∈ eachindex(x, y) 
        @inbounds (s[x[i]] += y[i])
    end
    return (s[1]/ns[1]) - (s[2]/ns[2]) # mean difference
end
```

That's it. 

Let's test our own test. First create some data

```julia
y=[4., 5., 6., 7., 1., 2., 3]
ns=[4, 3] # N1=4, N2=3
```

For the permuttaion vector the dedicated function `membership` is available:
```julia
x=membership(StudentT_IS(), ns) # = [repeat ([1], N1); repeat([2], N2)]
```

A two-directional exact test, as performed by PermutationTests.jl is 

```julia
println("t-test independent samples (bi-directional) using PermutationsTest.jl: ");
t1 = studentTestIS(y, ns)
```

Let us run the test we have created. Notice that argument `stat` of `_permTest!`
is `MeanDiff()`, while argument `asStat` is `AnovaF_IS()` 
The default `fstat` is julia abs() function, thus the test will be bi-directional
as for the test `t1` here above.

```julia
println("t-test independent samples (bi-directional) using our own statistic: ");
t2 = _permTest!(copy(x), y, ns, MeanDiff(), AnovaF_IS())
```
NB: we pass a copy of `x` as this vector will be permuted. 

---

We now illustrate with a simple example a nice opportunity we have to write 
efficient code to compute test statistics;
the above implementation of the `testStatistic` function computes the group numerosity (`ns=[N1, N2]`) 
at each permutation, however `ns` is invariant to permutations. 

Let us rewrite the function then using the `cpcd` argument, 
which will be passed to the `_permTest!` functions in order to 
(optionally) pre-compute `ns`:

```julia
function testStatistic(x, y, stat::MeanDiff; cpcd=nothing, kwargs...)
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
```

Running the test now can use the `cpcd` argument to avoid redundant computations:

```julia
t3 = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); cpcd=ns) 
```

You can check that the p-values obtained by all the above tests is identical:
```julia
t1.p≈t2.p≈t3.p ? (@info "It works!") : 
                 (@error "It does not work!")
```

To run a right-directional test:

```julia
t4 = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); fstat = identity)
```

To run a left-directional approximate test (see [`_permTest!`](@ref)):
```julia
t5 = _permTest!(copy(x), y, ns, MeanDiff(), StudentT_IS(); fstat = flip, switch2rand=1)
```

Etc. 
   
As an exercise, wrap the test we have created in a friendly function.

!!! tip "Nota Bene"
    We are not restricted to the use of a `cpcd` argument. In later examples we will show how to use 
    an arbitrary number of specially defined keywords in order to obtain efficient tests.

---


### Example 4: multiple comparisons correlation test

We create another version of the Pearson product-moment correlation test between a fixed variable given as a vector `x` and ``M`` variables given as a vector of ``M`` vectors `Y`. We want to test simultaneously 
all correlations between `x` and `Y[m]`, for ``m=1...M``.

The procedure to follow is the same with two differences: the function [`testStatistic`](@ref) will now compute 
each one of the ```M`` test statistics and we will finally call the [`_permMcTest!`](@ref) function
instead of the `_permMcTest!` function.

The permutations scheme for the correlation test is the one of the `PearsonR()` statistic, hence it is supported by *PermutationTests.jl* and we can reuse it (using the appropriate `asStat` argument of `_permMcTest!`).

First, let us import `testStatistic` from *PermutationTests.jl* and what we need for the computations:

```julia
using PermutationTests
import PermutationTests: testStatistic

using LinearAlgebra: dot, ⋅ # we are going to use the dot product
```

(1) Define a new `Statistic` type to name the test statistic you want to use, such as

```julia
struct MyPearsonR <: Statistic end
```

(2) Write a function to compute the ``i{th}`` observed and permuted test statistic. In order to obtain a faster test we will work with standardized variables, thus the Pearson correlation is given simply by the 
cross-product divided by N. We can write then

```julia
testStatistic(x, Y, i, stat::MyPearsonR; kwargs...) = (x ⋅ Y[i])/length(x)
```

That's it. 

Let's test our own test. First create some data as example:

```julia
N, M = 7, 100 # N=7 observations, M=100 hypotheses
x=randn(N);
Y=[randn(N) for m=1:M];
```

A two-directional exact test, as performed by PermutationTests.jl is obtained by

```julia
println("correlation test using PermutationsTest.jl: ");
t1 = rMcTest!(x, Y)
```

Let's run our own test. Remind that our test works only for standardized variables, thus we have to input standardized variables. We can use function [`μ0σ1`](@ref):

```julia
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

As an excercise, wrap our test in a friendly function that will standardize the input data before
running the test. 

---

### Example 5: univariate Chatterjee correlation

Examples (4) and (5) are not very interesting per se, as they are just equivalent forms
of tests already implemented in the API of *PermutationTests.jl*. Creating genuinely new tests
is not necesserely more complicated thought, as it is shown in this example.

The Chatterjee (2021) correlation ``ξ_{x→y}`` is an asymmetric, asymptotycally consistent 
estimator of how much a random variable ``y`` is a measurable function of ``x``. Thus,
we can use it to test a null hypothesis of the form 

``H_0: y≠f(x)``. 

The function ``f`` does not have to be monotonic, nor linear.
In the following, ``x`` and ``y`` are assumed continuous with no ties.

Let ``r`` be the ranks of variable ``y`` computed after sorting ``y`` by ascending 
order of ``x``. The coefficient is given by

``ξ_{x→y}= 1-\frac{3δ_{x→y}(r)}{n²-1}``, 

where

``δ_{x→y}(r)=\sum_{1}^{n-1}\left| r_{i+1} - r_i \right|``

``ξ_{x→y}`` takes value in interval ``[0, 1]`` asymptotically. The limits with finite sampling
are ``[-1/2+O(1/n), (n-2)/(n+1)]``. 

``ξ`` has a simple asymptotic Gaussian distribution given by ``√n*ξ \sim N(0, 2/5)``
and is computed in ``O(n\textrm{log}n)``.

**Reference**: S. Chatterjee (2021) a new coefficient of correlation, journal of the american statistical association 
116(536), 506-15 

---

Let us now construct a univariate permutation tests for ``ξ_{x→y}``. The following code block defines a function to create it:

```julia
using PermutationTests
using StatsBase: ordinalrank # install StatsBase.jl if you don't have it
import PermutationTests: testStatistic

struct Chatterjeeξ <: Statistic end

testStatistic(r, y, stat::Chatterjeeξ; kwargs...) = sum(abs, diff(r))

xiCorPermTest(x, y; kwargs...) =
   _permTest!(ordinalrank(y[sortperm(x)]), [], length(x), Chatterjeeξ(), PearsonR(); 
                fstat = flip, kwargs...)

```
As in [example 3](@ref "Example 3: univariate t-test for independent samples"), we import `testStatistic` and we define a type for the test statistic of interest as a [`Statistic`](@ref) type, which we have named `Chatterjeeξ`. Notice that we use ``δ_{x→y}(r)`` as an 
equivalent test statistic for ``ξ_{x→y}``. Notice also that while this is possible because
``δ_{x→y}(r)`` and ``ξ_{x→y}`` are a monotonic function of each other, their relationship is
inverted, which will be taken care of later.

Then, we declare a [`testStatistic`](@ref) method for computing the ``δ_{x→y}(r)`` quantity defined above for the observed and permuted data.

Finally, we invoke the [`_permTest!`](@ref) function. 

The test works this way: `_permTest!` takes as input the ``n`` ranks and list all ``n!``
permutations of the rank indices for an exact test, or a random number of them for an approximate test.
As the ranks are a permutation of ``(1.,...,n)`` themselves, at each permutation the `testStatistic`
function is invoked and the test statistic is computed directly on the permutation vector,
whcih is always passed by `_permTest!` to `testStatistic` as the first argument. 

Notice that the test of interest is meaningful only as a right-directional test for ``ξ_{x→y}``, but this is a 
left-directional test for ``δ_{x→y}(r)``, therefore we have used keyword argument `fstat = flip`.
This takes care of the aforementioned inverse relationship of our equivalents statistic.

Let's use the test we have created:

```julia
n = 100
x = randn(n)
y = randn(n)
t = xiCorPermTest(x, y) # is not signifant

y = x.^2+randn(length(x))/5
t1 = xiCorPermTest(x, y) # y is a function of x with little noise

```


---

### Example 6: multiple comparisons Chatterjee correlation

---

Next, let us define the multiple comparisons permutation test for ``ξ_{x→y}``. We wish to test simultaneously

``H_0(m): y_m≠f(x), m=1...M``

The code could be:

```julia
using PermutationTests
using StatsBase: ordinalrank # install StatsBase.jl if you don't have it
import PermutationTests: testStatistic

struct Chatterjeeξ <: Statistic end

testStatistic(p, Y, m::Int, stat::Chatterjeeξ; kwargs...) = 
    1.0-((3*sum(abs, diff(Y[m][p])))/(abs2(length(x))-1))

# Take as input x and Y=y1,...,yM
function xiCorPermMTest(x, Y; kwargs...) 
   p = sortperm(x)
   n = length(x)
   x_ = collect(Base.OneTo(n))
   return _permMcTest!(x_, [ordinalrank(y[p]) for y ∈ Y], n, Chatterjeeξ(), PearsonR(); 
                  fstat = identity, kwargs...)
end                
```
The first four lines are to be omitted if we continue writing on the same unit where we have written 
the univariate test. 

Again, we define a method for the `testStatistic` function and we invoke the the [`_permMcTest!`](@ref) function.

Here we have to adopt a slightly different strategy: `testStatistic` takes as input the permutation vector
and the ranks of the ``y_1,..,y_M`` vectors. It computes the ``ξ_{x→y}`` statistic (note: not an equivalent statistic)
for these vectors after sorting them all by **the same** sorting key at each permutation, which is the permutation vector passed by `_permMcTest!` to `testStatistic` as the `x_` vector.

Note that ``ξ_{x→y}`` can take negative values, thus we pass to `_permMcTest!` keyword argument `fstat = identity`;
this will result in the correct right-directional test we sought.

Let's use the test we have created:

```julia
n = 100
m = 10
x = randn(n)
Y = [randn(n) for i=1:m]
t = xiCorPermMTest(x, Y) # is not signifant

Y = [x.^2+randn(length(x))/5  for i=1:m]
t1 = xiCorPermMTest(x, Y) # y1,...ym are a function of x with little noise

```


---



### Example 7: univariate distance correlation

This and next example illustrate how we can crete a test for a complex test statistic,
still avoiding redundant computations, as if we wrote the code for a permutation test
from scratch.

The **distance correlation** (Skézely, Rizzo and Bakirov, 2007) extends the concept of correlation
to distances among objects living in metric spaces, yielding a measure of dependence that is somehow sensitive to 
linear and non-linear dependence. We will here limit ourselves to
scalars, vectors and matrices, it does not matter if real or complex.

Let ``{(X_k, Y_k): k=1...K}`` be ``K`` pairs of
observations (may be scalars, vectors, matrices). The elements within ``X`` and ``Y`` must have the same size, 
but the size of the elements in ``X`` may be different from the size of those in ``Y``.

Let ``D_x`` and ``D_y`` be the distance matrices with elements

``D_x^{(ij)}=\left| X_i - X_j \right|_p : i,j={1,...,K}``

and

``D_y^{(ij)}=\left| Y_i - Y_j \right|_p : i,j={1,...,K}``,

where ``\left|\cdot\right|_p`` is the p-norm, with ``p∈(0, 2]``.

Let 

``H=I-K^{-1}\textbf{1}\textbf{1}^T``

be the centering matrix, where ``I`` is the identity matrix and 
``\textbf{1}`` is the vector of ones.

Let

``A=HD_xH'``

and

``B=HD_yH'``

be the double-centered distance matrices. Then the **distance variance** is defined such as

``\nu(X)=\frac{1}{K^2}\sum_{i,j=1}^{K}a_{ij}^2``

and the **distance covariance** such as

``\nu(X, Y)=\frac{1}{K^2}\sum_{i,j=1}^{K}a_{ij}^2b_{ij}^2``.

Finally, the **distance correlation** is defined such as

``\rho(X, Y)=\frac{\nu(X, Y)}{\sqrt{\nu(X)\nu(Y)}}``,

where by convention it takes value zero if the denominator is zero.

**Reference**: G.J. Székely, M.L. Rizzo, N.K. Bakirov (2007). Measuring and testing dependence by correlation of distances Ann. Statist. 35(6): 2769-2794.

---
After giving the definition, let us see how we can create a univariate permutation tests
for the distance correlation.

First, let us define some functions we will need.

```julia

using LinearAlgebra: I, norm, Hermitian, LinearAlgebra

# Distance matrix for an array `X` of scalars, vectors, matrices, etc.
# The `n` elements of `X` can be real or complex.
# `pnorm` is the norm used to get distaces d_ij=norm(X[i]-X[j], pnorm),
# see julia LinearAlgenbra.norm function.
# Return the distance matrix filled only in the upper triangle.
# Note: this is always real.
function dm(x, n; pnorm=2)
    D = Matrix{Float64}(undef, n, n)
    for i=1:n-1 
        @simd for j=i+1:n
            @inbounds D[i, j] = norm(x[i]-x[j], pnorm)
        end
    end
    @simd for i=1:n 
        @inbounds D[i, i] = 0. 
    end
    return D
end    

# Distance variance. Eq. (2.9) in Székemy, Rizzo and Bakirov(2007). 
# D is a distance matrix
dVar(D) = sum(x->abs2(x), D)/size(D, 1)^2  

# Distance covariance. Eq. (2.8) in Székemy, Rizzo and Bakirov(2007)
# Dx and Dy are two distance matrices of equal size
dCov(Dx, Dy) = sum(Dx.*Dy)/size(Dx, 1)^2  

# Distance correlation. Eq. (2.10) in Székemy, Rizzo and Bakirov(2007)
# Dx and Dy are two distance matrices of equal size
function dCor(Dx, Dy)
    den = dVar(Dx) * dVar(Dy)
    return den > 0 ? sqrt(dCov(Dx, Dy)/sqrt(den)) : 0.  
end
```

---
We are now ready to create the univariate test.

First, let us import `testStatistic`:

```julia
using PermutationTests
import PermutationTests: testStatistic
```

Given two vectors of elements `x` and `y` with distance matrices `Dx` and `Dy`, the test is 
obtained permuting the indices of `x`. Instead of recomputing distance matrix `Dx` at each permutation, 
the strategy is to permute instead the rows and columns of `Dx` by a permutation 
matrix `P` that is changed at each permutation according to the permutation vector `x`.

Let us then write the function to do this:

```julia
# Overwrite P with a permutation matrix given a vector x holding a 
# permutation of n indices. P must exist, be square and have same dimension as x
function permMatrix!(P, x)
   fill!(P, 0.)
   @simd for i∈axes(P, 1)
         @inbounds P[i, x[i]] = 1. 
   end
   return P
end
```

As in the previous examples, we need to define a type for the distance correlation statistic:
```julia
struct Dcor <: Statistic end
```

Next, as in the previous example, we need to define a method for `testStatistic`.
This function takes as arguments a permutation vector `p` and `Dy`, the distance
matrces of `y`, which is fixed across permutations.

Furthermore, we will use four specially defined keyword arguments:
- `H` : the centering matrix
- `P` : the permutation matrix (to be updated at each call of the function)
- `Dx`: the original distance matrix for data `x`
- `dVarDy`: the distance variance of `y`

At each permutation then, we will permute and double-center `Dx` in order to 
compute the distance correlation. The quantities that are invariant by permutation
are passed to the function so that we do not need to recompute them.

```julia
function testStatistic(x, HDyH, stat::Dcor; H, P, Dx, dVarHDyH, kwargs...)
   permMatrix!(P, x)
   # Hx * P * Dx * P' * Hx', Dx with rows and columns permuted
   HP = H*P
   HxPDxPHx = HP * Dx * HP' # permuted and double centered Dx
   den = dVar(HxPDxPHx) * dVarHDyH # dVar(Dx) * dVar(Dy), square of the denominator of dCor 
   return den > 0 ? sqrt(dCov(HxPDxPHx, HDyH)/sqrt(den)) : 0. # dCor
end

```
Finally, we write a function preparing the data and calling the [`_permTest!`](@ref) function.
The preparation involves computing the distance matrices (once and for all),
making some checks, creating the permutation vector `p` and the keyword arguments
that will be passed internally to the function `testStatistic` we have created by `_permTest!`.

Importantly, the vector `p` is created so as to correspond to the permutation
yielding the observed statistic (``i.e.``, no data permutation), that is, in this case,
the vector with the indices in the natural order ``1...n``.

```julia
# The test takes as input two vectors of elements, which may be scalars, vectors or matrices.
function dCorPermTest(x, y; pnorm=2, kwargs...)
    n = length(x)
    Dx = Hermitian(dm(x, n; pnorm)) 
    Dy = Hermitian(dm(y, n; pnorm)) 
   
    size(Dx, 1) == size(Dx, 2) || throw(ArgumentError("Function dCorPermTest: Did you want to call dCorPermMTest intead? The `Dx` and `Dy` distance matrices that have been computed are not square"))
    size(Dx) == size(Dy) || throw(ArgumentError("Function dCorPermTest: Did you want to call dCorPermMTest intead? The `Dx` and `Dy` distance matrices that have been computed are not square or do not have the same size"))
    p = collect(Base.OneTo(n)) # permutation vector
   
    #  the centering matrix, the identity matrix, HDyH, dVar(HDyH)
    Id = Matrix{Float64}(I, n, n)
    H = Id - fill(1/n, n, n) # the centering matrix
    P = copy(Id) # the permutation matrix
    HDyH = H * Dy * H'
    dVarHDyH = dVar(HDyH) # the distance variace of HDyH (invariant to permutations)

    return _permTest!(p, HDyH, n, Dcor(), PearsonR(); fstat=identity, H, P, Dx, dVarHDyH, kwargs...)
end
```

We are done.

Let us use the function. We make an example where the elements of `x` and `y` are vectors.
You can verify yourself that the test works in the same way if they are scalars or matrices.

```julia
# Vectors
n=10
l=20
x=[randn(l) for i=1:n]
y=[randn(l) for i=1:n]
perm = dCorPermTest(x, y; pnorm=2, switch2rand=1)

# non-linear relationship (could be detected)
y=[x_.^2+randn(length(x_))./10 for x_ in x]
perm = dCorPermTest(x, y; switch2rand=1)
```

---

### Example 8: multiple comparisons distance correlation

This example reuses the code we have already written for the previous
[example 7](@ref "Example 7: univariate distance correlation").

The strategy here is the same as in that example, however here we have
to be more careful as there are several keyword arguments that will be
updated at each permutation.

The additional code we need to create a multiple comparisons permutation test
is here below. As compared to the previous example, note that:
- we are asking function `_permMcTest!` to pass internally to the function `testStatistic` the product `HPDxPH` as keyword argument because we do not want to compute it for every call of the `testStatistic` function, but only when it is called for the first hypothesis.
- we are similarly also passing the distance variance of Dx `dVarDx` as a keyword argument for the same reason.

Note also
- the syntax `[:]` to update variables passed as keyword 
- the fact that `dVarHDxH` is passed as a vector of one element so as to be possible to update it thanks to the syntax `[:]`.


```julia
function testStatistic(x, HDYH, m::Int, stat::Dcor; H=H, P=P, Dx, dVarHDYH, HxPDxPHx, dVarHDxH, kwargs...)
    if m==1
        permMatrix!(P, x)
        # Hx * P * Dx * P' * Hx', Dx with rows and columns permuted 
        HP = H * P
        HxPDxPHx[:] = HP * Dx * HP'  # notice the [:] syntax; this kwarg is updated when m=1
        dVarHDxH[:] = [dVar(HxPDxPHx)] # notice the [:] syntax; this kwarg is updated when m=1
    end
    
    den = dVarHDxH[1] * dVarHDYH[m] # dVar(Dx) * dVar(Dy), square of the denominator of dCor 
    return den > 0 ? sqrt(dCov(HxPDxPHx, HDYH[m])/sqrt(den)) : 0. # dCor
end
 

function dCorPermMTest(x, Y; pnorm=2, kwargs...)
    n = length(x)
    Dx = Hermitian(dm(x, n; pnorm)) 
    DY = [Hermitian(dm(y, n; pnorm)) for y ∈ Y] 
    size(Dx, 1) == size(Dx, 2) || throw(ArgumentError("Function dCorPermMTest: The `Dx` and `Dy` distance matrices that have been computed are not square"))
    size(Dx) == size(DY[1]) || throw(ArgumentError("Function dCorPermMTest: The `Dx` and `Dy` distance matrices that have been computed are not square or do not have the same size"))
    length(unique(size.(DY, 1))) == 1 || throw(ArgumentError("Function dCorPermMTest: All elements of second data input must have equal size"))
    n = size(Dx, 1)
    p = collect(Base.OneTo(n)) # permutation vector

    # keyword arguments that are not updated
    Id = Matrix{Float64}(I, n, n)
    H = Id - fill(1/n, n, n)
    P = copy(Id) 
    HDYH = [H*Dy*H' for Dy ∈ DY]
    dVarHDYH = dVar.(HDYH)

    # Initialize keyword arguments that will be updated
    HxPDxPHx = Matrix{Float64}(undef, size(Dx)...)
    dVarHDxH = [0.0] # NB, cannot pass a scalar as kwargs if it is to be updated!

    return _permMcTest!(p, HDYH, n, Dcor(), PearsonR(); 
                        fstat=identity, threaded=false, H, P, Dx, dVarHDYH, HxPDxPHx, dVarHDxH, kwargs...)
end

```
---

Let us use the multiple comaparions test we have just created.

**Example where `x` and `y` hold vectors:**
```julia
# Vectors
n=10
l=20
m=30
x=[randn(l) for i=1:n]
Y=[[randn(l) for i=1:n] for j=1:m]
# random data. No hypothesis should be significant
perm = dCorPermMTest(x, Y)

# The y1...ym variables are noisy copy of x.  
# The p-value should be significant for about the first M/2 hypothesis out of M.
# The actual number of rejected hypotheses could be slighty different then M/2. 
for i=1:m÷2
    Y[i]=[x[j].+randn(l)./2 for j=1:n]
end
perm = dCorPermMTest(x, Y)
```


**Example where `x` and `y` hold matrices:**
```julia
# Matrices
n=10
l=20
m=30
x=[randn(l, l) for i=1:n]
Y=[[randn(l, l) for i=1:n] for j=1:m]
# random data. No hypothesis should be significant
perm = dCorPermMTest(x, Y)

# The y1...ym variables are noisy copy of x.  
# The p-value should be significant for about the first M/2 hypothesis out of M
# The actual number of rejected hypotheses could be slighty different then M/2. 
for i=1:m÷2
    Y[i]=[x[j].+randn(l, l)./2 for j=1:n]
end
perm = dCorPermMTest(x, Y)

```
---

In conclusion, in the above examples we have illustrated how to create permutation tests of increasing complexity.
The last two examples concerns tests that are pretty complex and very different from any of the test
implemented in *Permutationtests.jl*, where even the permutation scheme had to be defined differently.
The diverse procedures exposed here above can be adapted to new problems in order to
create a countless number of new permutation tests.

## Useful functions for creating your own tests

---

```@docs
_permTest!
_permMcTest!
flip
membership
testStatistic
```
