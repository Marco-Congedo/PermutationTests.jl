# Main module

Here below are the types declared in the package that you may need to use.

---
#### TestDirection

```julia
abstract type TestDirection end 
```
This is the abstract type for the three possible *test directionalities*. Instances of this type are julia [singletons](https://docs.julialang.org/en/v1/manual/types/#man-singleton-types) and are used as arguments of test functions.

| Type | Test Directionality |  
|:----------|:----------|
| `Both` | bi-directional test (default for all tests) |
| `Right` | right-tailed test |
| `Left` | left-tailed test |

---
#### Assignment

```julia
abstract type Assignment end
```
This is the abstract type for the two possible *test designs*. Instances of this type are julia [singletons](https://docs.julialang.org/en/v1/manual/types/#man-singleton-types). The test design is determined automatically by the test functions and is one of the field of the structure returned by the test functions.

| Type | Test Design |  
|:----------|:----------|
| `Balanced` | equal number of subjects in all groups/tratments/blocks... |
| `Unbalanced`  | unequal number of subjects in the groups/tratments/blocks... |

---
#### TestResult

```julia
abstract type TestResult end
```
This is the abstract type for the structures returned by all test functions.

---
#### UniTest

All functions carrying out [univariate tests](@ref "Univariate tests") return the following structure of type [TestResult](@ref):

```julia
struct UniTest <: TestResult
    p::Float64
    stat::statistic where statistic <: Statistic
    obsstat::Float64 
    minp::Float64 
    nperm::Int64 
    testtype::Symbol 
    direction::testDirection where testDirection <: TestDirection
    design::assignment where assignment <: Assignment
end    

```

**Fields**:

 - `.p`: the p-value of the test. 
 - `.stat`: the actual test-statistic that has been used as a singleton, see [Statistic](@ref). 
 - `.obsstat`: observed statistic, the value of the test-statistic `.stat` for the input data. 
 - `.nperm`: the number of permutations used to obtain the p-value. 
 - `.minp`: the minimum attainable p-value, given by 1 / `nperm`. 
 - `.testtype`: either `:exact` or `:approximate`, depending on whether all possible or a random number of permutations, respectively, have been used. 
 - `.direction`: the test-directionality set by the user as a singleton, see [TestDirection](@ref). 
 - `.design`: the test design, as established by the test function as a singleton, see [Assignment](@ref). 

---

#### MultcompTest

All functions carrying out [multiple comparisons tests](@ref "Multiple comparisons tests") return the following structure of type [TestResult](@ref):

```julia
struct MultcompTest <: TestResult
    p::Vector{Float64} 
    stat::statistic where statistic <: Statistic
    obsstat::Vector{Float64} 
    minp::Float64 
    nperm::Int64 
    testtype::Symbol 
    direction::testDirection where testDirection <: TestDirection
    design::assignment where assignment <: Assignment
    nulldistr::Vector{Float64}
    rejections::Vector{Vector{Int64}} 
    stepdown::Bool 
    fwe::Float64
end    
```

**Fields**:

The fields `.stat`, `.minp`, `.nperm`,`.testtype`, `.direction` and `.design` are the same as per the [UniTest](@ref) test result structure. 

For the others, special to multiple comparisons tests, say the test concerned ``M`` multiple-comparisons hypotheses, used ``P`` permutations and the step down procedure performed ``S`` steps. ``S`` will be equal to one if the step down procedure is not used or if it stops after the first step:

 - `.p`: the vectors of ``M`` p-values of the test. 
 - `.obsstat`: observed statistics, the values of the ``M`` test-statistics `.stat` for the input data. 
 - `.nulldistr`: if ``S`` is equal to 1 this is the vector holding the null-distribution of the test-statistic computed by data permutation, otherwise it is the vector holding the null-distribution obtained at the last step of the step down procedure.  
 - `.rejections`: a vector holding ``S`` vectors, each one holding the indeces of the hypotheses that have been rejected at the corresponding step of the step down procedure. 
 - `.stepdown`: `true` if the step down procedure has been used, `false` otherwise. 
 - `.fwe`: the family-wise error rate set by the user. This is used only if `.stepdown` is `true`. By default it is fixed to `0.05` for all multiple comparisons tests.

---
#### Statistic

```julia
abstract type Statistic end
```
This is the abstract type for all *test statistics*. This type is useful if you want to [create your own test](@ref "Create your own test"). If you don't, you can skip the rest of this page. 

All test-statistic are julia [singletons](https://docs.julialang.org/en/v1/manual/types/#man-singleton-types), that is, immutable composite type (struct) with no fields. For each typical test-statistic used for paramtric tests, there exist an associated *equivalent test-statistic* that gives the same p-value with a particular univariate permutation test. This happens because the test-statistic can usually be decomposed in several terms, some of which are invariant to data permutation for a given test directionality (see [TestDirection](@ref)) and/or test design (see [Assignment](@ref)). 

Univariate tests in *Permutations.jl* make use of this equivalence to always yield the fastest test possible. In the following table, we list the usual test-statistic used for parametric tests in **bold** and the equivalent test-statistic used for univariate permutation tests just below, along with the test directionality and/or test design for which they apply. 


| Test-statistic | Type | When it applies|  
|:----------|:----------|:----------|
| **Pearson correlation coefficient r** | `PearsonR` |  always |
| Cross Product | `CrossProd` |  directional test* |
| Covariance| `Covariance` |  bi-directional test |
| **1-way ANOVA for independent samples F** | `AnovaF_IS` |  always |
| Sum of group totals squared | `SumGroupTotalsSq_IS` |  balanced design |
| Sum of group totals squared/N| `SumGroupTotalsSqN_IS` |  always |
| **Student t for independent samples** | `StudentT_IS` |  always |
| Sum of group totals squared | `Group1Total_IS` |  directional |
| Sum of group totals squared/N| `SumGroupTotalsSqN_IS` |  bi-directional test |
| **1-way ANOVA for repeated measures F** | `AnovaF_RM` |  always |
| Sum of treatement totals squared | `SumTreatTotalsSq_RM` |  always |
| **One-sample Stedent t**| `StudentT_1S` |  always |
| Sum of observations| `Sum` |  always |
| |  |  * always if the data is standardized |

The equivalent statistic that has been employed for a given univariate test is returned as one of the field
of the [test result](@ref "UniTest").

To know a-priori what equivalent statistic will be used for a test, see [`eqStat`](@ref).

Note that equivalent statistics are not used for multiple comparisons tests, as they are not dimensionless
and their use could unduly favor some hypotheses over others.

Test statistics are computed internally by the test functions of the package.
If you need to compute them for other purposes, see the page on [test statistics](@ref "Test statistics").

---

#### Statistic groups

All implemented test-statistics are grouped in four union types:

 - `BivStatistic` : all bivariate test-statistics
 - `IndSampStatistic` : all indepedent samples test-statistics
 - `RepMeasStatistic` : all repeated measures test-statistics
 - `OneSampStatistic` : all one-sample test-statistics, actually a special case of `RepMeasStatistic`.

Regardless the actual test-statstic used in the permutation test, the union type it belongs to indicates 
the permutation scheme it is subjected to, 
which in turn determines the number of possible data permutations.
The test statistics composing the four groups are listed in the table here below:

| `BivStatistic`| `IndSampStatistic`    | `RepMeasStatistic`    | `OneSampStatistic`|  
|:--------------|:----------------------|:----------------------|:------------------|
|`PearsonR`     |`AnovaF_IS`            |`AnovaF_RM`            |`StudentT_1S`      |
|`CrossProd`    |`SumGroupTotalsSq_IS`  |`SumTreatTotalsSq_RM`  |`Sum`              |
|`Covariance`   |`SumGroupTotalsSqN_IS` |                       |                   |
|               |`StudentT_IS`          |                       |                   |
|               |`Group1Total_IS`       |                       |                   |
|               |`SumGroupTotalsSqN_IS` |                       |                   |


```julia

# Singletons
StudentT_1S() isa IndSampStatistic # false

StudentT_1S() isa OneSampStatistic # true

StudentT_1S() isa StudentT_1S # true

# These are types instead and are never used as such
StudentT_1S <: IndSampStatistic # false

StudentT_1S <: OneSampStatistic # true
```

