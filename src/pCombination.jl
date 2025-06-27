# v 0.2 June 2025
# Part of the PermutationTests.jl package.
# Copyright Marco Congedo, CNRS, University Grenoble Alpes.

#=
see: Nicholas A. Heard and Patrick Rubin-Delanchy(2017) 
Choosing Between Methods of Combining p-values 
https://arxiv.org/pdf/1707.06897
=#

import Distributions: Normal, Beta, Chisq, cquantile, cdf, ccdf

# all p-value combinations are right-sided given a vector of p-values.
# They all return the probability to observe a set of p-values as small
# as those observed.

# Liptak. Equivalent to Stouffer computed on z-scores. Distributed as a Normal(0, 1).
# Fisher. Multiplicative combination. Distributed as a chi²(2). Sensitive to the smallest p-values
# Pearson. Multiplicative combination. Distributed as a chi²(2). Sensitive to the largest p-values
# Tippett (min). distributed as a Beta(1, n). Sensitive only to the smallest p-value


"""
```julia
combine(::F, p::AbstractVector{T}) where {F<: PcombFunc, T<:Real}
```

Combination of the p-values given as a vector `p` of real numbers according to combination function `F`.
    
**Examples**
```julia
using PermutationTests

N=100
p=[rand() for n=1:N] # random p-values under H0

pf = combine(FisherComb, p)

pl = combine(LiptakComb, p)

# Liptak Meta-Combination of pf and pl
pm = combine(LiptakComb, [pf, pl])
```
"""
combine(::Type{LiptakComb}, 𝐩::AbstractVector{T}) where T<:Real =
        ccdf(Normal(), sum(cquantile.(Normal(), 𝐩))/sqrt(length(𝐩)))

combine(::Type{FisherComb}, 𝐩::AbstractVector{T}) where T<:Real =
        ccdf(Chisq(2*length(𝐩)), -2.0*sum(log, 𝐩))

combine(::Type{PearsonComb}, 𝐩::AbstractVector{T}) where T<:Real =
        cdf(Chisq(2*length(𝐩)), -2.0*sum(log, (1.0.-𝐩)))

combine(::Type{TippettComb}, 𝐩::AbstractVector{T}) where T<:Real =
        cdf(Beta(1, length(𝐩)), minimum(𝐩))
              
