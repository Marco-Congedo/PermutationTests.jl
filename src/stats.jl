#=

Stats.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home

#############################
# Low-Level Test Statistics #
#############################

Run `test_stats()` in unit _test.jl to test this unit

========================================
  EXPORTED:

SCALARS functions of vectors
---------------------------------------------------------------
∑(𝐱)                  sum (sum_i (x_i))
∑of²(𝐱)               sum of squares (sum_i(x_i²))
∑∑of²(𝐱)              sum and sum of squares in one pass
μ(𝐱)                  arithmetic average (1/N ∑(𝐱))
dispersion(𝐱)         sum of squares of deviations from the mean
σ²(𝐱)                 variance (1/N sum_i(x_i-μ(𝐱))²)
σ(𝐱)                  standard deviation sqrt (1/N sum_i(x_i-μ(𝐱))²)
∑(𝐱, 𝐲)               sum of 𝐱 + 𝐲 (sum_i (x_i y_i))
Π(𝐱)                  product (product_i(x_i))
∑ofΠ(𝐱, 𝐲)            sum of products (sum_i (x_i*y_i))


statistic(𝐱, 𝐲, [statistic]; standardized=false, centered=true, means=(), sds=(), k, ns, ∑Y²kn, ∑y², ∑S²k) 
                      compute a statistic of independent variable 𝐱 and dependent variable 𝐲.
    [statistic] can be any statistic: CrossProd, Covariance, PearsonR, .....
    Depending on the statistic, some kwargs may be passed to speed up computations


VECTOR functions of vectors
---------------------------------------------------------------
μ0(𝐱)           recenter (zero mean: 𝐱-μ(𝐱)𝟏)
σ1(𝐱)           equalize (unit standard deviation: 𝐱/σ(𝐱))
μ0σ1(𝐱)         standardize (zero mean and unit standard deviation: (𝐱-μ(𝐱)𝟏)/σ(𝐱))
_∑y²(𝐲)
_∑Y²kn_∑y²_∑S²k

OTHERS
---------------------------------------------------------------
_∑Y²kn
_∑S²k


  NON-EXPORTED UTILITIES:
_getmeans
_anova_IS_fixed_params
_getsds
_cond∑_IS
_cond∑_RM
_anova_IS_fixed_params
=============================
=#


############################################
# Low-level statistical computations
############################################

# sum
"""
```julia
function (x::UniData)
```

Alias for julia `sum(x)`. Sum of the elements in `x`.
"""
∑(𝐱::UniData) = sum(𝐱)


# sum of squares
"""
```julia
function ∑of²(x::UniData)
```

Sum of squares of the elements in `x`
"""
∑of²(𝐱::UniData) = sum(abs2, 𝐱)


# sum and sum of squares in one pass
"""
```julia
function ∑∑of²(x::UniData)
```

Compute the sum and sum of squares of the elements in `x` in one pass and return them as a 2-tuple.

*Examples*
```julia
using PermutationTests
x=randn(10)
s, s² = ∑∑of²(x)
```
"""
function ∑∑of²(𝐱::UniData)
    s, s² = 0., 0.
    @simd for x ∈ 𝐱
        @inbounds s += x
        @inbounds s² += abs2(x)
    end
    s, s²
end


# mean
"""
```julia
function μ(x::UniData) 
```

Alias for julia `mean(x)`. Arithmetic mean of the elements in `x`.
"""
μ(𝐱::UniData) = mean(𝐱)


# create a centered copy of 𝐱, giving or not the mean to go faster 
"""
```julia
function μ0(x::UniData; 
    mean::Realo=nothing) 
```

Return `x` centered (zero mean).

Center around `mean` if it is provided. `mean` must be a real number. 
"""
μ0(𝐱::UniData; mean::Realo=nothing) =  mean === nothing ? 𝐱.-μ(𝐱) : 𝐱.-mean  


# sum of squares of deviations from the mean, with the mean given
"""
```julia
function dispersion(x::UniData, mean::Real)
```
Sum of squared deviations of the elements in `x` around `mean`.

`mean` must be provided and must be a real number.

*Examples*
```julia
using PermutationTests
x=randn(10)
d=dispersion(x, μ(x))
```
"""
function dispersion(𝐱::UniData, mean::Real)
    d=0.
    @simd for x in 𝐱
        @inbounds d += abs2(x-mean)
    end
    d
end


"""
```julia
function dispersion(x::UniData; 
    centered::Bool=false, mean::Realo=nothing)
```

If `centered` is true, return the sum of the squares of the elements in `x`,

otherwise, 

if `mean` is provided return the sum of squared deviations of the elements in `x` from `mean`,

otherwise,

return the sum of squared deviations of the elements in `x` from their mean.

*Examples*
```julia
using PermutationTests
x=randn(10)
x0=μ0(x)
μx=μ(x)
m=μ(x)

# in decreasing order of efficiency
d1=dispersion(x0)
d2=dispersion(x, mean=μx)
d3=dispersion(x)
d1==d2==d3 # -> true
```
"""
dispersion(𝐱::UniData; centered::Bool=false, mean::Realo=nothing) = 
    centered ? ∑of²(𝐱) : (mean===nothing ? dispersion(𝐱, μ(𝐱)) : dispersion(𝐱, mean)) 


"""
```julia
function σ²(x::UniData; 
    corrected::Bool=false, 
    centered::Bool=false, 
    mean::Realo=nothing) 
```

Population variance of the ``N`` elements in x (default), or its unbiased sample estimator if `corrected` is true.

If `centered` is true return the sum of the squared elements in `x` divided by ``N``, or ``N-1`` if 
`corrected` is true,

otherwise,

call julia standard function `var` with arguments `corrected` and `mean`.
"""
σ²(𝐱::UniData; corrected::Bool=false, centered::Bool=false, mean::Realo=nothing) = 
        centered ? ∑of²(𝐱)/(corrected ? (length(𝐱)-1) : length(𝐱)) : var(𝐱; corrected, mean) 


"""
```julia
function σ(x::UniData; 
    corrected::Bool=false, 
    centered::Bool=false, 
    mean::Realo=nothing) 
```

Population standard deviation of the elements in `x` or its unbiased sample estimator.

Call [`σ²`](@ref) with the same arguments and return its square root.
"""
σ(𝐱::UniData; corrected::Bool=false, centered::Bool=false, mean::Realo=nothing) = 
    sqrt(σ²(𝐱; corrected, centered, mean))   


"""
```julia
σ1(x::UniData; 
    corrected::Bool=false, 
    centered::Bool=false, 
    mean::Realo=nothing, 
    sd::Realo=nothing)
```

Return `x` equalized (unit standard deviation).

If `sd` is provided, return `x` equalized with respect to it,

otherwise call [`σ`](@ref) with optional arguments `corrected`, `centered` and `mean` and equalize `x`.

*Examples*
```julia
using PermutationTests
x=randn(3)
μx=μ(x) # mean
σx=σ(x) # sd
xμ0=μ0(x) # centered data

# in decreasing order of efficiency
e1=σ1(xμ0; sd=σx, centered=true)
e2=σ1(xμ0; centered=true)
e3=σ1(x; mean=μx, sd=σx)
e4=σ1(x; mean=μx)

println(σ(e4)≈1.0 ? "OK" : "Error") # -> "OK"
```
"""
σ1(𝐱::UniData; corrected::Bool=false, centered::Bool=false, mean::Realo=nothing, sd::Realo=nothing) =
    𝐱./(sd===nothing ? σ(𝐱; corrected, centered, mean) : sd)

    
"""
```julia
μ0σ1(x::UniData; 
    corrected::Bool=false, 
    centered::Bool=false, 
    mean::Realo=nothing, 
    sd::Realo=nothing)
```

Return `x` standardized (zero mean and unit standard deviation).

If `centered` is true, equalize `x` by calling [`σ1`](@ref) with optional arguents `corrected`, `centered`
and `sd`,

otherwise standardize `x` after computing the mean, if it is not provided as `mean`, and the standard deviation,
if it is not provided, calling `σ`[@ref] with optional arguments `corrected` and `centered`. 

*Examples*
```julia
using PermutationTests
x=randn(3)
μx=μ(x) # mean
σx=σ(x) # sd
xμ0=μ0(x) # centered data

# in decreasing order of efficiency
z1=μ0σ1(xμ0; sd=σx, centered=true)
z2=μ0σ1(xμ0; centered=true)
z3=μ0σ1(x; mean=μx, sd=σx)
z4=μ0σ1(x; mean=μx)
z5=μ0σ1(x; sd=σx)
z6=μ0σ1(x)

println(μ(z6)≈0. ? "OK" : "Error") # -> "Ok"
println(σ(z6)≈1. ? "OK" : "Error") # -> 'OK"
```
"""
μ0σ1(𝐱::UniData; corrected::Bool=false, centered::Bool=false, mean::Realo=nothing, sd::Realo=nothing) =
    if centered 
        return σ1(𝐱; corrected, centered=true, sd=(sd===nothing ? σ(𝐱; corrected, centered=true) : sd))
    else
        m = mean === nothing ? μ(𝐱) : mean
        s = sd === nothing ? σ(𝐱; corrected, centered=false, mean=m) : sd
        return [(x-m)/s for x in 𝐱] # can use Folds.map for this
    end
   
"""
```julia
function Π(x::UniData)
```

Product of the elements in `x`. Alias of julia function `prod`.
"""
Π(𝐱::UniData) = prod(𝐱)


"""
```julia
function ∑ofΠ(x::Union{UniData, Tuple}, y::UniData) 
```

Inner product of `x` and `y`.
"""
function ∑ofΠ(𝐱::Union{UniData, Tuple}, 𝐲::UniData) 
    p = 0.
    @simd for i in eachindex(𝐲) # here eachindex(𝐱, 𝐲) gives an error
        @inbounds p += 𝐱[i]*𝐲[i]
    end
    p
end



############################################
# Test Statistics
############################################


# --------------------------------------------------------------------
# Pearson Product Moment R and equivalent Cross-Product and Covariance
# --------------------------------------------------------------------

_getmeans(𝐱::UniData, 𝐲::UniData; means::Tuple=()) =
    isempty(means) ? (μ(𝐱), μ(𝐲)) : (means[1], means[2])


_getsds(𝐱::UniData, 𝐲::UniData; 
        centered::Bool=false, means::Tuple=(), sds::Tuple=()) =
    if isempty(sds) # standard deviations
        if centered
            return σ(𝐱, centered=true), σ(𝐲, centered=true) 
        else
            return isempty(means) ? (σ(𝐱), σ(𝐲)) : (σ(𝐱, mean=means[1]), σ(𝐲, mean=means[2]))
        end
    else
        return sds[1], sds[2]
    end 

"""
```julia
function statistic(x::UniData, y::UniData, stat::PearsonR; 
                standardized::Bool=false, 
                centered::Bool=false, 
                means::Tuple=(), 
                sds::Tuple=(), kwargs...)
```
Pearson product-moment correlation *r* statistic of input data vector `x` and `y`. 

If the means and/or standard deviations of `x` and `y` are passed as tuple `means` and `sds`, respectively, they are not computed.

If `centered` is true, the two vectors are assumed to have zero mean.

If `standardized` is true both x and y are assumed standardized, 
thus only the cross-product needs to be computed. This is by far the most efficient way 
if this function is to be called repeatedly on the same data input vectors and for computing
p-values by data permutations as the cross-product is en equivalent test statistic for the Pearson correlation 
if the data is standardized.

*Examples*
```julia
using PermutationTests
x, y = randn(10), randn(10);
c1 = statistic(x, y, PearsonR())

zx=μ0σ1(x);
zy=μ0σ1(y);
c2 = statistic(zx, zy, PearsonR(); standardized=true)
```

see [`μ0σ1`](@ref)

"""
statistic(𝐱::UniData, 𝐲::UniData, stat::PearsonR; 
        standardized::Bool=false, centered::Bool=false, means::Tuple=(), sds::Tuple=(), kwargs...) =
    standardized ? ∑ofΠ(𝐱, 𝐲) / length(𝐱) :
    statistic(𝐱, 𝐲, Covariance(); centered, means, kwargs...) / prod(_getsds(𝐱, 𝐲; centered, means, sds))


"""
```julia
function statistic(x::UniData, y::UniData, stat::Covariance; 
                centered::Bool=false, 
                means::Tuple=(), 
                kwargs...)
```
Covariance of data input vectors x and y.

Equivalent to the Pearson correlation *r* test statistic for bi-directional correlation tests, 
see [Statistic](@ref).

For the optional keyword arguments, see the method [`statistic(x, y, stat::PearsonR)`](@ref).

*Examples*
```julia
using PermutationTests
x, y = randn(10), randn(10)
c1 = statistic(x, y, Covariance())

c2 = statistic(μ0(x), μ0(y), Covariance(); centered=true)

μx=μ(x)
μy=μ(y)
c3 = statistic(x, y, Covariance(); means=(μx, μy))
```

see [`μ0`](@ref), [`μ`](@ref) 

"""
statistic(𝐱::UniData, 𝐲::UniData, stat::Covariance; 
        centered::Bool=false, means::Tuple=(), kwargs...) =
    centered ? ∑ofΠ(𝐱, 𝐲)/length(𝐱) : (∑ofΠ(𝐱, 𝐲)/length(𝐱))-prod(_getmeans(𝐱, 𝐲; means))


"""
```julia
function statistic(x::UniData, y::UniData, stat::CrossProd; 
                kwargs...)
```

Cross-product (inner product) of data input vectors x and y.

Equivalent to the Pearson correlation *r* test statistic for directional correlation tests and always 
if the data is standardized, see [Statistic](@ref).

*Examples*
```julia
using PermutationTests
x, y = randn(10), randn(10)
c = statistic(x, y, CrossProd())
```
"""
statistic(𝐱::UniData, 𝐲::UniData, stat::CrossProd; kwargs...) = ∑ofΠ(𝐱, 𝐲) 

# --------------------------------------------------------------------
# ANOVA for Independent Samples F and equivalent sum of group totals, 
# sum of square group total and mean, Student T (all used also for Chi², 
# Fisher Exact Test statistic and point bi-serial correlation tests)
# --------------------------------------------------------------------

# Return 2 tuple (k, ns), where k is the number of groups and ns is the numerosity vector
# If askwargs=true is passed, return instead the named tuple (k=k, ns=ns)
function _anova_IS_fixed_params(𝐱::IntVec; askwargs=false) 
    k=length(unique(𝐱)) 
    ns=[count(x->x==i, 𝐱) for i=1:k]
    return askwargs ? (k=k, ns=ns) : (k, ns)
end

"""
```julia
function statistic(x::IntVec, y::UniData, stat::AnovaF_IS; 
                k::Into=nothing, 
                ns::Union{Vector{Int}, Nothing}=nothing, 
                kwargs...) 
```

*F* statistic of the 1-way ANOVA for independent samples, see [Edgington (1995)](@ref "References"), p. 60. 

The data is given as an unique vector `y`, holding the observations for each group in the natural order,
thus for ``K`` groups `y` holds ``N=N_1+ \\ldots +N_K`` elements, where ``N_k`` is the number of observations 
for the ``k^{th}`` group.

`x` is the [`membership(::IndSampStatistic)`](@ref) vector.

`k` is the number of groups (independent samples). If omitted it will be computed.

`ns` is the group numerosity vector with form `ns=[N1,..., NK]`. If omitted it will be computed.

*Examples*
```julia
using PermutationTests
g1=[1.0, 2.0, 3.0, 4.0] # observations for group 1 
g2=[1.1, 2.8, 3.2, 4.4, 5.3] # observations for group 2 
ns=[length(g1), length(g2)] # N1, N2
y=vcat(g1, g2)
x=membership(AnovaF_IS(), ns)
F=statistic(x, y, AnovaF_IS())
# or, to avoid finding k and ns from y
F2=statistic(x, y, AnovaF_IS(); k=2, ns=ns)

# The square of the t test statistic for independent samples for a bi-directional test
# is equal to the above F test statistics

t=statistic(x, y, StudentT_IS(); ns=ns)
println(t^2≈F ? "OK" : "error")

```
"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::AnovaF_IS; 
            k::Into=nothing, ns::Union{Vector{Int}, Nothing}=nothing, kwargs...) 
    k===nothing && (k=length(unique(𝐱))) 
    ns===nothing && (ns=[count(x->x==i, 𝐱) for i=1:k])
    n=length(𝐱)
    ∑𝐲s = Vector{Float64}(undef, k)
    SSTv = 0.
    @inbounds for i=1:k
        e=𝐲[𝐱.==i] # it is better to allocate e here and perform sum and sum of squares sequentially on e
        s, s² = ∑∑of²(e)
        ∑𝐲s[i] = s
        SSTv += s²
    end
    ms  = sum(∑𝐲s)^2/n
    SSB = sum(∑𝐲s[i]^2/ns[i] for i=1:k) - ms # sum of squares between groups Edgington, page 60
    SST = SSTv - ms # Total sum of squares
    SSW = SST - SSB # sum of squares within groups
    if SSW ≈ 0. 
        return Inf #throw(ErrorException(📌*"Function statistic(x, y, AnovaF_IS()): The sum of squares within is equal to zero"))
    else
        return (SSB*(n-k))/(SSW/(k-1)) # F statistic : (SSB / (k-1)) / (SSW / (n-k)) 
    end
    # NB some dicothomous tables yield a sum of squares within SSW = 0, thus F=Inf
end

# Sum the elements of 𝐲 if the corresponding elements of 𝐱 are equal to k
# For Ind Samp ANOVA-like statistics the data is arranged as s1g1, s2g1,..., sn1g1, s2g1, s2g2,...,sn2g2,... ... snNgK.
# This function then return the total for a group indexed by k. 
# 𝐱 must be of the form [repeat(1, N1); repeat(2, N2),..., repeat(2, NK)] or any valid permutation, for example 
# for N1=3, N2=2 and K=2 the observed statistic is given by 𝐱=[1 1 1 2 2] and a valid permutation is [2 1 2 1 1]
# NOT SAFE
function _cond∑_IS(𝐱::IntVec, 𝐲::UniData, k::Int)
    s=0.
    @simd for i ∈ eachindex(𝐱) # before it was i=1:length(𝐱)
        @inbounds 𝐱[i]==k && (s += 𝐲[i])
    end
    s
end

"""
```julia
function statistic(x::IntVec, y::UniData, stat::SumGroupTotalsSq_IS; 
                k::Into=nothing, 
                kwargs...)
```
Sum of squared group totals. 

Equivalent to the 1-way ANOVA for independent samples *F* test statistic in balanced designs, see [Statistic](@ref). 

For optional keyword argument `k`, see the [`statistic`](@ref) method for `AnovaF_IS` here above.

"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::SumGroupTotalsSq_IS; 
        k::Into=nothing, kwargs...)
    k===nothing && (k=length(unique(𝐱)))
    sum(abs2, (_cond∑_IS(𝐱, 𝐲, i) for i=1:k)) 
end

#= Robust version of SumGroupTotalsAbs_IS. It does not yield the same p-value as the F statistic
function statistic(𝐱::IntVec, 𝐲::UniData, stat::SumGroupTotalsAbs_IS; 
        k::Into=nothing, kwargs...)
    k===nothing && (k=length(unique(𝐱)))
    sum(abs, (_cond∑_IS(𝐱, 𝐲, i) for i=1:k)) 
end
=#

# Equivalent to StudentT_IS for bi-diretional tests. 
"""
```julia
function statistic(x::IntVec, y::UniData, stat::SumGroupTotalsSqN_IS; 
                k::Into=nothing, 
                ns::Union{Vector{Int}, Nothing}=nothing, 
                kwargs...)
```
Sum of squared group totals divided each by the group numerosity. 

Equivalent to the 1-way ANOVA for independent samples *F* test statistic in general, see [Statistic](@ref).

For optional keyword arguments `k` and `ns` and for examples, 
see the [`statistic`](@ref) method for `AnovaF_IS` here above.

"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::SumGroupTotalsSqN_IS; 
        k::Into=nothing, ns::Union{Vector{Int}, Nothing}=nothing, kwargs...)
    k===nothing && (k=length(unique(𝐱))) 
    ns===nothing && (ns=[count(x->x==i, 𝐱) for i=1:k])
    sum(((_cond∑_IS(𝐱, 𝐲, i) for i=1:k) .|> abs2)./ns) 
end


"""
```julia
function statistic(x::IntVec, y::UniData, stat::StudentT_IS; 
                k::Into=nothing, 
                ns::Union{Vector{Int}, Nothing}=nothing, 
                kwargs...)
```

Student's *t* statistic for independent samples, see [Edgington (1995)](@ref "References"), p. 3.

For all arguments except `stat` and examples see the [`statistic`](@ref) method for `AnovaF_IS` here above.


"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::StudentT_IS; 
        k::Into=nothing, ns::Union{Vector{Int}, Nothing}=nothing, kwargs...)
    k===nothing && (k=length(unique(𝐱))) 
    ns===nothing && (ns=[count(x->x==i, 𝐱) for i=1:k])
    length(ns)≠2 && throw(ArgumentError(📌*"Function statistic (StudentT_IS): `ns` must contain 2 integer"))
    𝐠1, 𝐠2 = 𝐲[𝐱.==1], 𝐲[𝐱.==2] # data for the two groups
    μ1, μ2 = μ(𝐠1), μ(𝐠2) # mean for the two groups
    pσ² = (dispersion(𝐠1, μ1) + dispersion(𝐠2, μ2)) / (ns[1] + ns[2] - 2) # pooled variance
    return (μ1 - μ2) / sqrt(pσ²/ns[1] +pσ²/ns[2])  
end

"""
```julia
function statistic(x::IntVec, y::UniData, stat::Group1Total_IS; 
                kwargs...)            
```

Sum of observations for group 1.

Equivalent to the Students's *t* test statistic for independent samples for diretional tests, see [Statistic](@ref).

For arguments `x` and `y` and for examples, see the [`statistic`](@ref) method for `AnovaF_IS` here above.

"""
statistic(𝐱::IntVec, 𝐲::UniData, stat::Group1Total_IS; kwargs...) = _cond∑_IS(𝐱, 𝐲, 1) 


# --------------------------------------------------------------------
# Repeated Measures Statistics
# --------------------------------------------------------------------

# Sum the elements of 𝐲 if mod1(x, k)==k, where x are the elements of 𝐱 corresponding to the elements of 𝐲.
# For Rep Meas ANOVA-like statistics the data is arranged as s1t1, s1t2, ..., s1tk, s2t1, s2t2,...,s2tk,... ... sntk.
# This function then return the total for each treatement k. 
# 𝐱 must be of the form [1,..., nk]o r any valid permutation, for example for n=2, k=3 
# the observed statistic is given by 𝐱=[1 2 3 4 5 6] and a valid permutation is [2 1 3 6 5 4].
# Note that the permutations are restricted within successive k elements, i.e., inside the bars: | 2 1 3 | 6 5 4 | ... |.
function _cond∑_RM(𝐱::IntVec, 𝐲::UniData, k::Int, ns::@NamedTuple{n::Int, k::Int})
    s=0.
    @simd for i ∈ k:ns.k:ns.n*ns.k
        @inbounds s += 𝐲[𝐱[i]]
    end
    s
end


# ∑ of al observations, squared and dividen by k*n. For computing 1-way Rep Meas ANOVA F statistic
_∑Y²kn(𝐲, ns) = abs2(sum(𝐲))/(ns.k*ns.n)

# Pre-compute some data for the [`statistic`](@ref) methods computing the `StudentT_1S` test-statistic
# and the `AnovaF_RM` test-statistic.
# Return the sum of all squared observations in `y`. 

_∑y²(𝐲) = ∑of²(𝐲) 

# sum of subject totals; squared and divided by n. For computing 1-way Rep Meas ANOVA F statistic
_∑S²k(𝐲, ns) = sum(abs2(sum(view(𝐲[i:i+ns.k-1], :)))/ns.k for i=1:ns.k:ns.n*ns.k)  

# The three above in one pass as a vector
# Pre-compute some data for the [`statistic`](@ref) computing the `AnovaF_RM` test-statistic.
# Return a vector of three elements [_∑Y²kn, `_∑y²`, `_∑S²k`].
_∑Y²kn_∑y²_∑S²k(𝐲, ns) = [abs2(sum(𝐲))/(ns.k*ns.n), ∑of²(𝐲), sum(abs2(sum(view(𝐲[i:i+ns.k-1], :)))/ns.k for i=1:ns.k:ns.n*ns.k)]


"""
```julia
function statistic(x::IntVec, y::UniData, stat::AnovaF_RM; 
                ns::@NamedTuple{n::Int, k::Int}, 
                ∑Y²kn::Realo=nothing, 
                ∑y²::Realo=nothing, 
                ∑S²k::Realo=nothing, 
                kwargs...) 
```
*F* statistic of 1-way ANOVA for repeated measures, see [Edgington (1995)](@ref "References"), p. 102. 

The data is given as an unique vector `y` concatenaning the ``N`` observations for the ``K`` measures 
(treatments, time, ...) in the natural order, that is, the ``K`` treatments for observation 1, ..., 
the ``K`` tratments for observation ``N``. Thus, `y` holds ``N \\cdot K`` elements. 

`x` is the [`membership(::RepMeasStatistic)`](@ref) vector.

`ns` is a julia [named tuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) 
with form `(n=N, k=K)` (see examples below).

`∑Y²kn`, `∑y²` and `∑S²k` can be optionally provided to speed up computations since these quantities are
invariant by data permutations. The exported function `_∑Y²kn_∑y²_∑S²k` can be used for this purpose, 
see the examples below. 

*Examples*
```julia
using PermutationTests
# K=2 (measurements), N=4 (observations)
o1=[1.0, 2.0]; # first observation 
o2=[2.0, 2.8]; # second observation 
o3=[3.0, 2.6]; # third observation
o4=[4.0, 4.1]; # fourth observation
ns=(n=4, k=2); # four obs. and two measurements
y=vcat(o1, o2, o3, o4);
x=membership(AnovaF_RM(), ns);
F=statistic(x, y, AnovaF_RM(); ns=ns)

# pre-compute some data
pcd=_∑Y²kn_∑y²_∑S²k(y, ns);
F2=statistic(x, y, AnovaF_RM(); ns=ns, ∑Y²kn=pcd[1], ∑y²=pcd[2], ∑S²k=pcd[3])

# The t test statistic for repeated measures is the same as the one-sample 
# t test statistic on the difference of the two measurements. 
# The square of those statistics for a bi-directional test are equal to 
# the above F test statistics. 

x=membership(StudentT_1S(), ns.n)
t=statistic(x, y[1:2:ns.n*2-1].-y[2:2:ns.n*2], StudentT_1S())
println(t^2≈F ? "OK" : "error")

```
"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::AnovaF_RM; 
        ns::@NamedTuple{n::Int, k::Int}, ∑Y²kn::Realo=nothing, ∑y²::Realo=nothing, ∑S²k::Realo=nothing, kwargs...) 

    # quantities that are invariant by permutation
    ∑Y²kn   === nothing && (∑Y²kn = _∑Y²kn(𝐲, ns)) 
    ∑y²     === nothing && (∑y² = _∑y²(𝐲))  
    ∑S²k    === nothing && (∑S²k = _∑S²k(𝐲, ns))

    ∑T²n = sum(abs2(_cond∑_RM(𝐱, 𝐲, i, ns))/ns.n for i=1:ns.k)
    SSB = ∑T²n-∑Y²kn
    SSe = ∑y²-∑T²n-∑S²k+∑Y²kn
    return (SSB/SSe)*(ns.n-1) # F statistic : (SSB / (k-1)) / (SSE / (n-1)(k-1)) 
end


"""
```julia
function statistic(x::IntVec, y::UniData, stat::SumTreatTotalsSq_RM; 
                ns::@NamedTuple{n::Int, k::Int}, 
                kwargs...)
```
Sum of squared treatment totals.

Equivalent to the 1-way ANOVA for repeated measures *F* test statistic in general, see [Statistic](@ref).

For optional keyword argument `ns` see the [`statistic`](@ref) method for `AnovaF_RM` here above.

"""
function statistic(𝐱::IntVec, 𝐲::UniData, stat::SumTreatTotalsSq_RM; 
        ns::@NamedTuple{n::Int, k::Int}, kwargs...)
    sum(abs2, (_cond∑_RM(𝐱, 𝐲, i, ns) for i=1:ns.k)) 
end

# --------------------------------------------------------------------
# One-sample Statistics
# --------------------------------------------------------------------
# NB: these statistics behave differently for exact and approximate tests.
# For exact tests 𝐱 is a tuple. For approximate tests it is the useual UniData type, but it is ignored

# 1 sample t-test. See https://en.wikipedia.org/wiki/Student's_t-test 

# StudentT_1S; H0: μ=0. Optionally provide ∑y² to go faster.
function _studentT_1S(𝐲::UniData; ∑y²::Realo=nothing)
    n = length(𝐲)
    m = μ(𝐲)
    if ∑y² === nothing
        s = σ(𝐲; mean=m, corrected=true)
        return (m * sqrt(n))/s
    else
        ∑y²n = ∑y²/n
        m² = abs2(m)
        diff = ∑y²n - m²
        if abs(diff) < 1e-1 # try to avoid catastrophic cancellation
            s = σ(𝐲; mean=m, corrected=true)
            return (m * sqrt(n))/s
        else
            s = sqrt(diff*(n/(n-1))) # ∑y² is invariant by permutations 
            return (m * sqrt(n))/s
        end
    end
end

# this version taking a tuple as argument is needed for systematic permutations 
# as the iterator for permutations yields tuples
"""
```julia
function statistic(x::Tuple, y::UniData, stat::StudentT_1S; 
                ∑y²::Realo=nothing, 
                kwargs...) 
```
Student's one-sample *t* statistic.

`y` is the input data.

`x` is a tuple holding as many 1.0 as elements in `y`.

`∑y²` can be optionally provided to speed up computations, since this quantity is
invariant by data permutations. The exported function `_∑y²` can be used for this purpose, 
see the examples below. 

*Examples*
```julia
using PermutationsTest
y=randn(6);
x=(1., 1., 1., 1., 1., 1.);
t=statistic(x, y, StudentT_1S()) 

pcd=_∑y²(y)
t2=statistic(x, y, StudentT_1S(); ∑y²=pcd) 
println(t≈t2 ? "OK" : "Error")
```

"""
statistic(𝐱::Tuple, 𝐲::UniData, stat::StudentT_1S; ∑y²::Realo=nothing, kwargs...) =
    _studentT_1S(𝐱 .* 𝐲; ∑y²) 

"""
```julia
function statistic(x::UniData, y::UniData, stat::StudentT_1S; 
                ∑y²::Realo=nothing, 
                kwargs...)
```

Student's one-sample *t* statistic.

`y` is the input data.

`x` is the [`membership(::OneSampStatistic)`](@ref) vector.

`∑y²` can be optionally provided to speed up computations since this quantity is
invariant by data permutations. The exported function `_∑y²` can be used for this purpose, 
see the examples below. 

*Examples*
```julia
using PermutationsTest
y=randn(6);
x=membership(StudentT_1S(), length(y));
t=statistic(x, y, StudentT_1S()) 

pcd=_∑y²(y)
t2=statistic(x, y, StudentT_1S(); ∑y²=pcd) 
println(t≈t2 ? "OK" : "Error")
```
"""    
statistic(𝐱::UniData, 𝐲::UniData, stat::StudentT_1S; ∑y²::Realo=nothing, kwargs...) =
    _studentT_1S(𝐲; ∑y²)

# Sum (used for exact tests)
"""
```julia
function statistic(x::Tuple, y::UniData, stat::Sum; 
                kwargs...) 
```
"""
statistic(𝐱::Tuple, 𝐲::UniData, stat::Sum; kwargs...) = 
    ∑ofΠ(𝐱, 𝐲) # safer to use (𝐱 ⋅ 𝐲)

# Sum (used for approximate tests)
"""
```julia
function statistic(x::UniData, y::UniData, stat::Sum; 
                kwargs...)  
```

Sum of the elements in input data vector `y`.

Equivalent to the one-sample *t* test statistic in general, see [Statistic](@ref).

`x` is the [`membership(::OneSampStatistic)`](@ref) vector.


"""
statistic(𝐱::UniData, 𝐲::UniData, stat::Sum; kwargs...) = 
    sum(𝐲)

# maximum of `obsStats` considering only those elements which corrisponding element in `accepted` is true.
# used by _permTest! in unit multcompTest.jl
function _condMax(obsStats::UniData, accepted::BitArray)
    m = 0
    @simd for i ∈ eachindex(obsStats, accepted)
        @inbounds accepted[i] && obsStats[i]>m && (m = obsStats[i])
    end
    return m
end

####################################################################################


