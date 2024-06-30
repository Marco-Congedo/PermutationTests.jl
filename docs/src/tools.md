# Tools

Collection of utilities. At times, you will be referred to this page by the documentation,
but for general usage of this package you can skip this page.

---

## ns

The argument `ns` is often employed in the functions below, but also by some functions carrying out permutation tests.
This argument is special to each [group of statistics](@ref "Statistic groups"). In the table here below
``N`` denotes the number of observations and ``K`` the number of groups or measurements:

| Statistics group  | `ns` argument |  
|:------------------|:--------------|
|`BivStatistic`     | an integer indicating the number of bivariate observations ``N`` |
|`IndSampStatistic` | a vector of integers, holding the group numerosity ``[N_1,..., N_K]`` |
|`RepMeasStatistic` | a tuple with form ``ns=(n=N, k=K)`` |
|`OneSampStatistic` | an integer indicating the number of observations ``N`` |

---


```@docs
assignment
eqStat
allPerms
nrPerms
genPerms
table2vec
```

