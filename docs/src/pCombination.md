# p-values Combination

Combination functions for p-values.

Given ``N`` independent p-values, a **combination** of these p-values tests the omnibus null hypothesis ``H_0``
that they are drawn from a uniform ``[0, 1]`` distribution ([Nicholas et al., 2017](@ref "References")). The combination is a p-value itself,
returing the probability to observe a set of p-values as small as those observed.

!!! tip "Usage"
    This is useful whenever there are several independent observational units, each with its own distributon functions.
    For example, we may want to perform a statistical test for each subject separately, still we wish an evidence for an effect
    at the group level.

Each combination function is sensitive to a specific configuration of the observed p-values under the alternative hypothesis.
In particular:
- The **Liptak** function is sensitive to small departures from ``H_0`` of many p-values.
- The **Fisher** function is sensitive to large departures from ``H_0`` of a few of the smallest p-values.
- The **Pearson** function is sensitive to large departures from ``H_0`` of a few of the largest p-values.
- The **Tippett** combination is sensitive to large departures from ``H_0`` of the smallest p-value only.

This can be appreciated plotting the above combination functions in the case of ``N=2``. In the figure here below,
we plot the combination of all possible values in ``[0, 1]`` of two p-values (x and y axes) according to the **Liptak**, **Fisher**,
**Pearson** and **Tippett** function:

![Figure 1_pCombination](assets/Fig1_pCombination.png)

!!! tip "Meta-combinations"
    Several p-value combinations can be combined again (meta-combination) in order to define special rejection regions.

Each combination function has its own distribution under ``H_0``:
- **Liptak** is distributed as a *standard normal distribution*.
- **Fisher** and **Pearson** are distributed as a *χ2* (*chi²*) with ``2N`` degrees of freedom.
- **Tippett** is distributed as a *Beta* with parameter *α=1* and *β=N*.

!!! tip "Stouffer combination"
    The Liptak combination of p-values derived from z-statistics is equivalent to the Stouffer's combination of the z-statistics.

The p-value combination functions are abstract types:
```julia
abstract type LiptakComb <: PcombFunc end
abstract type FisherComb <: PcombFunc end
abstract type PearsonComb <: PcombFunc end
abstract type TippettComb <: PcombFunc end
```

They are all children of abstract type
```julia
abstract type PcombFunc end
```

```@docs
combine
```