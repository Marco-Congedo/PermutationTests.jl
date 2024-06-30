# Univariate tests

They are used when a single null hypothesis is to be tested. The following tests are available.

---

#### Univariate test functions

| Test |  Function  |  Alias  |
|:-----|:----------------|:---------------|
| Pearson correlation |  [`correlationTest`](@ref)| `rTest`|
| Trend correlation | [`trendTest`](@ref)| |
| Point bi-serial correlation | [`pointBiSerialTest`](@ref) |  |
| Student's t for independent samples | [`studentTestIS`](@ref) | `tTestIS` |
| 1-way ANOVA for independent samples | [`anovaTestIS`](@ref) | `fTestIS` |
| Chi-squared | [`chiSquaredTest`](@ref) | `Χ²Test` |
| Fisher exact| [`fisherExactTest`](@ref)| |
| Student's t for repeated measures | [`studentTestRM`](@ref) | `tTestRM` |
| 1-way ANOVA for repeated measures  | [`anovaTestRM`](@ref) | `fTestRM` |
| Cochran Q | [`cochranqTest`](@ref) | `qTest` |
| McNemar| [`mcNemarTest`](@ref)| |
| One-sample Student's t | [`studentTest1S`](@ref) | `tTest1S` |
| Sign | [`signTest`](@ref) |  |

---


For creating other tests, see [Create your own test](@ref).

For multiple comparisons tests, see [Multiple comparisons tests](@ref).

---

#### Common kwargs for univariate tests
The following optional keyword arguments are common to all univariate test functions:

 - `direction`: an instance of [TestDirection](@ref), either `Right()`, `Left()` or `Both()`. The default is `Both()`. 
 - `equivalent`: a boolean. If `true` (default), the fastest equivalent statistic will be used. See [Statistic](@ref). 
 - `nperm`: an integer providing the number of random permutations to be used for an approximate test. It defaults to `20000`. 
 - `switch2rand`: an integer setting the upper limit of permutations to be listed exhaustively. It defaults to `1e8`. If the number of possible permutations exceeds `switch2rand`, the approximate test with `nperm` random permutations will be performed, otherwise an exact test with all possible permutations will be performed. In order to force an approximate test, set `switch2rand` to a small integer such as `1`. In order to know in advance the number of possible permutations, see [`nrPerms`](@ref). 
 - `seed`: an integer. It applies only to approximate tests. Set to `0` to use a random seed for generating random permutations. Any natural number results instead in a reproducible test. It defaults to `1234`. 
 - `verbose`: a boolean. Print some information in the REPL while running the test. Set to false if you run benchmarks. The default is `true`.

---

#### Univariate tests API

---
## Correlation test
```@docs
correlationTest
correlationTest!
```

---
## Trend test
```@docs
trendTest
trendTest!
```

---
## Point bi-serial correlation test
```@docs
pointBiSerialTest
```

---
## Student's t-test for independent samples
```@docs
studentTestIS
```

---
## 1-way ANOVA for independent samples
```@docs
anovaTestIS
```

---
## Chi-squared test
```@docs
chiSquaredTest
```

---
## Fisher exact test
```@docs
fisherExactTest
```

---
## Student's t-test for repeated measures
```@docs
studentTestRM
studentTestRM!
```

---
## 1-way ANOVA for repeated measures
```@docs
anovaTestRM
```

---
## Cochran Q test
```@docs
cochranqTest
```

---
## McNemar test
```@docs
mcNemarTest
```

---
## One-sample Student's t-test
```@docs
studentTest1S
studentTest1S!
```

---
## Sign test
```@docs
signTest
```
