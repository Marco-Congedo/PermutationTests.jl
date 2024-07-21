# Multiple comparisons tests

With a multiple comparisons test ``M`` null hypotheses are tested simultaneously. 
By deriving the null distribution with the max-statistic approach, the following permutation tests
control the family-wise error [(FWE)](https://en.wikipedia.org/wiki/Family-wise_error_rate) rate, that is, the probability to commit one or more [Type I error](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors) at the nominal level:

---

#### Multiple comparisons test functions

| test |  use function  |  alias  |
|:-----|:----------------|:---------------|
| Pearson correlation |  [`correlationMcTest`](@ref)| `rMcTest`|
| Trend correlation | [`trendMcTest`](@ref)| |
| Point bi-serial correlation | [`pointBiSerialMcTest`](@ref) |  |
| Student's t for independent samples | [`studentMcTestIS`](@ref) | `tMcTestIS` |
| 1-way ANOVA for independent samples | [`anovaMcTestIS`](@ref) | `fMcTestIS` |
| Chi-squared | [`chiSquaredMcTest`](@ref) | `Χ²McTest` |
| Fisher exact| [`fisherExactMcTest`](@ref)| |
| Student's t for repeated measures | [`studentMcTestRM`](@ref) | `tMcTestRM` |
| 1-way ANOVA for repeated measures  | [`anovaMcTestRM`](@ref) | `fMcTestRM` |
| Cochran Q | [`cochranqMcTest`](@ref) | `qMcTest` |
| McNemar| [`mcNemarMcTest`](@ref)| |
| One-sample Student's t | [`studentMcTest1S`](@ref) | `tMcTest1S` |
| Sign test | [`signMcTest`](@ref) |  |

You may also find useful the tests we have created as examples of how to create new tests:

| Test | 
|:----------|
| [Post-hoc tests for 1-way repeated-measures ANOVA](@ref "Example 2: Post-hoc tests for 1-way repeated-measures ANOVA") |
| [Chatterjee correlation](@ref "Example 6: multiple comparisons Chatterjee correlation") |
| [Distance correlation](@ref "Example 8: multiple comparisons distance correlation") |


---

For creating other tests, see [Create your own test](@ref).

For univariate tests, see [Univariate tests](@ref).

---

#### Common kwargs for multiple comparisons tests
The following optional keyword arguments are common to all multiple comparisons test functions:

 - `direction`: an instance of [TestDirection](@ref), either `Right()`, `Left()` or `Both()`. The default is `Both()`. 
 - `nperm`: an integer providing the number of random permutations to be used for an approximate test. It defaults to `20000`. 
 - `switch2rand`: an integer setting the upper limit of permutations to be listed exhaustively. It defaults to `1e8`. If the number of possible permutations exceeds `switch2rand`, the approximate test with `nperm` random permutations will be performed, otherwise an exact test with all possible permutations will be performed. In order to force an approximate test, set `switch2rand` to a small integer such as `1`. In order to know the number of possible permutations, see [`nrPerms`](@ref). 
 - `seed`: an integer. It applies only to approximate tests. Set to `0` to use a random seed for generating random permutations. Any natural number results instead in a reproducible test. It defaults to `1234`. 
 - `verbose`: a boolean. Print some information in the REPL while running the test. Set to false if you run benchmarks. The default is `true`.
 - `stepdown` : a boolean. If true (default) the step-down procedure is used. This is at least as powerfu as the standard procedure and still controls the FWE.
 - `fwe` : a real number in (0, 1). This is used by the step-down prcedure to control the family-wise error (FWE) rate at this level. By default it is 0.05.
 - `threaded` : a boolean. If true (default) some computations are multi-threaded if the product of the number of hypotheses, observations, and permutations exceed 500 millions.

---

#### Multiple comparisons tests API

---
## Correlation test
```@docs
correlationMcTest
correlationMcTest!
```

---
## Trend test
```@docs
trendMcTest
trendMcTest!
```

---
## Point bi-serial correlation test
```@docs
pointBiSerialMcTest
```

---
## Student's t-test for independent samples
```@docs
studentMcTestIS
```

---
## 1-way ANOVA for independent samples
```@docs
anovaMcTestIS
```


---
## Chi-squared test
```@docs
chiSquaredMcTest
```

---
## Fisher exact test
```@docs
fisherExactMcTest
```

---
## Student's t-test for repeated measures
```@docs
studentMcTestRM
studentMcTestRM!
```

---
## 1-way ANOVA for repeated measures
```@docs
anovaMcTestRM
```

---
## Cochran Q test
```@docs
cochranqMcTest
```

---
## McNemar test
```@docs
mcNemarMcTest
```

---
## One-sample Student's t-test
```@docs
studentMcTest1S
studentMcTest1S!
```

---
## Sign test
```@docs
signMcTest
```
