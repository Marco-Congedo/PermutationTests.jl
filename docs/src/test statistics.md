# Test statistics

The test statistics used by *PermutationTests.jl* are computed by the several methods of the 
`statistic` function listed here below.
You do not need to use them to carry out permutation tests, however they are exported 
because they may turn useful. For general usage of this package you can skip this page.

All `statistic` methods take a singleton of the [Statistic](@ref) type and other arguments to yield a specific test-statistic.

All test statistics are grouped in four [Statistic groups](@ref).

---

```@docs
statistic
```

