# PermutationTests.jl

---

## Installation

Execute the following command in Julia's REPL:

    ]add PermutationTests

---
#### Requirements 

**Julia**: version ≥ 1.10

---
#### Dependencies

| standard Julia packages |     external packages    |
|:-----------------------:|:-----------------------|
| [Random](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Random) |  [Combinatorics](https://github.com/JuliaMath/Combinatorics.jl)|
| [Statistics](https://bit.ly/2Oem3li) | [Folds](https://github.com/JuliaFolds/Folds.jl) |
| [Test](https://docs.julialang.org/en/v1/stdlib/Test/#Test.@test) |  |



---
## Quick start

Here are some quick examples to show how *Permutationtests.jl* works:

---
**Univariate test for correlation**

Given two vectors of ``N`` observations, ``x`` and ``y``, test the null hypothesis 

``H_0: r_{(x,y)}=0``, 

where ``r_{(x,y)}`` is the [Pearson product-moment correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between ``x`` and ``y``. 

First, we chose a [Type I error](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors) ``α`` , typically fixed to Ronald Fisher's classical 0.05 level.

```@example
using PermutationTests
N=10 # number of observations
x, y = randn(N), randn(N) # some random Gaussian data for example
t = rTest(x, y)
```
We reject the null hypothesis if the p-value is less then or equal to ``α``, otherwise we suspend any judgement.

The result of the test is a structure, which fields are printed in yellow. For example:

```julia
t.p # p-value
t.nperm # number of permutations used by the test
```

---
**Multiple comparison test for correlation**

Given a vector of ``N`` observations ``x`` and ``M`` vectors of ``N`` observations ``y_m``, test simutaneously the ``M`` null hypotheses 

``H_0:r_{(x, y_m)}=0, m=1...M``, 

where ``r_{(x,y_m)}`` is the Pearson product-moment correlation between ``x`` and the ``m^{th}`` vector ``y_m``. 

First, we fix the family-wise error [(FWE)](https://en.wikipedia.org/wiki/Family-wise_error_rate) rate we are willing to tolerate, say, at level 0.05.

```@example
using PermutationTests
N=10 # number of observations
M=100 # number of hypotheses
x, Y = randn(N), [randn(N) for m=1:M] # random Gaussian data
t = rMcTest(x, Y) # by default the FWE is 0.05
```

For each one of the ``m`` hypotheses, we reject the null hypothesis if its p-value is smaller than the nominal level (0.05), otherwise we suspend any judgement. 

The result of the test is a structure, which fields are printed in yellow. For example:

```julia
t.p # p-values
t.obsstat # observed test statistics
```

---
## Preamble

If you have no idea what a [`statstical hypothesis test`](https://en.wikipedia.org/wiki/Statistical_hypothesis_test) is, most probably you are on the wrong web page.

If you do, but you have no idea what a [permutation test](https://en.wikipedia.org/wiki/Permutation_test) is, check first my [introduction to permutation tests](https://sites.google.com/site/marcocongedo/science/tutorials?authuser=0).

If you need help to establish what hypothesis test is appropriate for your data, 
check [this page](@ref "Chose a test").

Permutation tests offer a different way to obtain the p-value usually obtained with well-known parametric tests, such as the *t-test for independent samples*, the *ANOVA for repeated mesures*, etc. In contrast to parametric tests, permutation tests may provide the *exact* p-value for the test, not an approximated value based on probability distribution theory, and make use of less stringent assumptions. 

Moreover, using permutation tests it is possible to use any test-statistic, not just those with known distribution under the null hypothesis, as such distribution is obtained by data permutation.

Permutation tests have been introduced by none other than R.A. Fisher and E.J.G Pitman in the late '30 
(see the [references](@ref "References")), but has become feasable only thanks to the advent of modern computers.

When multiple hypotheses are to be tested simultaneously, the [multiple comparisons problem](https://en.wikipedia.org/wiki/Multiple_comparisons_problem) arises: statistical tests perforormed on each hypothesis separatedly cannot control anymore the probability to falsely reject the null hypothesis ([Type I error](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors)). Using the *max-statistic* (also known as *min-p*) approach, permutation tests control the probabilty to commit one or more Type I error over the total number of hypotheses tested, that is, they control the family-wise error [(FWE)](https://en.wikipedia.org/wiki/Family-wise_error_rate) rate (see [Westfall and Young, 1993](@ref "References")).

While other multiple comparisons correction procedures controlling the FWE exists, such as the well-known [Bonferroni](https://en.wikipedia.org/wiki/Bonferroni_correction) or [Holm-Bonferroni](https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method), they assume independence of the hypotheses. Instead, permutation tests do not.
Actually, they naturally adapt to *any degree and form of correlation among the hypotheses*, thus they result more powerful when the hypotheses are not independent (see [power](@ref "Power")).

Thanks to these characteristics, permutation tests conctitutes the ideal framework for large-scale, possibly correlated, statistical hypothesis testing. This is the case, for instance of *genomics* (gene expression level) and *neuroimaging* (brain voxel activation), where up to hundreds of thousands of hypotheses, often largely correlated, must be tested simultaneusly.

---

## Overview

*PermutationTests.jl* implements several *univariate* permutation tests and for all of them the corresponding *multiple comparisons* permutation tests based on *max-statistic*, with or without the *step down* procedure
([Holmes et *al.*, 1996; Westfall and Young, 1993](@ref "References")). 

For multiple comparisons tests, only tha case when all hypotheses are homogeneous is considered (same test statistic, same number of observations divided in the same number of groups/measurements).

Here is the list of available tests:


| Univariate and Multiple Comparisons Permutation Tests | 
|:----------|
| Pearson product-moment correlation | 
| Trend correlation (fit of any kind of regression) |
| Point bi-serial correlation* |
| Student's t for independent samples | 
| 1-way ANOVA for independent samples | 
| Χ² for 2xK contingency tables* |
| Fisher exact test* (2x2 contingency tables) | 
| Student's t for repeated-measures | 
| 1-way ANOVA for repeated-measures | 
| Cochran Q*|
| McNemar*|
| One-sample Student's t  | 
| Sign*|
|                * for dicothomous data |

When the number of permutations is high, *PermutationTests.jl* computes approximate (Monte Carlo) p-values. All test functions switches automatically from systematic to Monte Carlo permutations, but the user can force them to use either one permutation listing procedure.

For all kinds of implemented univariate tests, *PermutationTests.jl* always employs the most computationally effective equivalent test-statistic and only the required (non-redundant) systematic permutations, yielding very efficient permutation tests (see [Benchmarks](@ref)).


