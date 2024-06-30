#=
tools.jl Unit of the PermutationTests.jl Package

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home


#######################################
# General Tools for Permutation tests #
######################################

========================================
  EXPORTED:

assignment      Return Balanced() or Unbalanced() for a given test test
eqStat          equivalent statistics
allPerms        number of permutations for a given test
nrPerms         number of non-redundant permutations for a given test
membership      create the membership vector for Independent Samples ANOVA tests and equivalent
genPerms        generate systematic permutations
table2vec       Convert a table to vectors to be used as input to test functions.
flip            flip the sign for a Number and negate for a Bool

  UTILITIES
_fstat
_direction
_check_ns - used by _permTest!
_test_params - used by _permTest!
_prepare_permtest! - used by _permTest!

=============================
=#


# ----- #
"""
```julia
flip(x::Bool)

flip(x::Union{R, I}) where {R<:Real, I<:Int}
```

Invert the sign of a real number or integer and negate a boolean. 

This function is be needed only to call [`_permTest!`](@ref) or [`_permMcTest!`](@ref) if you 
[create your own test](@ref "Create your own test"). 
"""
flip(x::Bool) = !x

flip(x::Union{R, I}) where {R<:Real, I<:Int} = -x
 
# ----- #


# ----- #
# Function to be applied to the permutation statistic before comparing it to the observed one so as to obtain 
# a test with direction Left, Right or Both. See [TestDirection](@ref) and _permTest!.
_fstat(direction::Left) = flip

_fstat(direction::Right) = identity

_fstat(direction::Both) = abs
# ----- #


# inverse function of _fstat
# ----- #
function _direction(fstat::Function) 
     if fstat === abs 
        return Both()
     elseif fstat === identity
        return Right()
     elseif fstat === flip
        return Left()
     else throw(ArgumentError(📌*" Function _direction: invalid input"))
     end
end
# ----- #


# ----- #
"""
```julia
function assignment(stat::Union{BivStatistic, OneSampStatistic}, ns::Int) 

function assignment(stat::IndSampStatistic, ns::Vector{Int}) 

function assignment(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int})
```

Analyse the [ns](@ref) argument for test statistic given as singleton `stat` and return a 
singleton of type [Assignment](@ref), `Balanced()` if the design is balanced 
(equal subjects in all groups/measurements), `Unbalanced()` otherwise. 

Only for test-statistics of the `IndSampStatistic` [group](@ref "Statistic groups") the result may be `Unbalanced()`; 
for all the others the result is always `Balanced()`. 

The complete list of test statistics is [here](@ref "Statistics").

*Examples*
```julia
using PermutationTests
assignment(PearsonR(), 12) # -> Balanced()
assignment(AnovaF_IS(), [5, 5, 6]) # -> Unbalanced()
assignment(AnovaF_IS(), [6, 6, 6]) # -> Balanced()
assignment(StudentT_IS(), [5, 6]) # -> Unbalanced()
assignment(StudentT_IS(), [5, 5]) # -> Balanced()
assignment(AnovaF_RM(), (n=10, k=3)) # -> Balanced()
assignment(StudentT_1S(), 8) # -> Balanced()
```

"""
assignment(stat::Union{BivStatistic, OneSampStatistic}, ns::Int) = Balanced()

assignment(stat::IndSampStatistic, ns::Vector{Int}) = length(unique(ns))==1 ? Balanced() : Unbalanced()

assignment(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int}) = Balanced()
# ----- #


# ----- #
"""
```julia
function eqStat(stat, direction::TestDir, design::Assign) 
    where {TestDir<: TestDirection, Assign <: Assignment}
```

Return the most efficient statistic equivalent to `stat`, for the given 
singleton `direction` of type [TestDirection](@ref) and singleton `design` of type [Assignment](@ref).

`stat` is a singleton of one of the five typical test statistics, that is, `PearsonR()`, `AnovaF_IS()`, `StudentT_IS`, `AnovaF_RM()` or 
`StudentT_1S()`, see [Statistic](@ref).

To avoid errors, use the result of [`assignment`](@ref) as argument `design`. 

Do not use `StudentT_RM()` as `stat`; *PermutationTests.jl* never use this test statistic. instead,
t-tests for repeated measures are carried out as one-sample t-tests on the difference of the two measurements. 

*Examples*
```julia
using PermutationTests
a=assignment
ns=12
eqStat(PearsonR(), Both(), a(PearsonR(), ns))       # -> Covariance()
eqStat(PearsonR(), Right(), a(PearsonR(), ns))      # -> CrossProd()
ns=[6, 6, 6]
eqStat(AnovaF_IS(), Both(), a(AnovaF_IS(), ns))     # -> SumGroupTotalsSq_IS()
ns=[6, 7, 6]
eqStat(AnovaF_IS(), Both(), a(AnovaF_IS(), ns))     # -> SumGroupTotalsSqN_IS()
ns=[6, 7]
eqStat(StudentT_IS(), Both(), a(StudentT_IS(), ns)) # -> SumGroupTotalsSqN_IS()
eqStat(StudentT_IS(), Left(), a(StudentT_IS(), ns)) # -> Group1Total_IS()
ns=(n=10, k=3)
eqStat(AnovaF_RM(), Both(), a(AnovaF_RM(), ns))     # -> SumTreatTotalsSq_RM()
eqStat(AnovaF_RM(), Right(), a(AnovaF_RM(), ns))    # -> SumTreatTotalsSq_RM()
ns=7
eqStat(StudentT_1S(), Both(), a(StudentT_1S(), ns)) # -> Sum()
eqStat(StudentT_1S(), Left(), a(StudentT_1S(), ns)) # -> Sum()
```
"""
eqStat(stat::PearsonR, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    direction isa Both ? Covariance() : CrossProd()

eqStat(stat::AnovaF_IS, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    design isa Balanced ? SumGroupTotalsSq_IS() : SumGroupTotalsSqN_IS()

eqStat(stat::StudentT_IS, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    direction isa Both ? SumGroupTotalsSqN_IS() : Group1Total_IS()
   
eqStat(stat::AnovaF_RM, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    SumTreatTotalsSq_RM()

eqStat(stat::StudentT_RM, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    StudentT_RM() # No eq stat is developed in this case as the test should be done as a one-sample test on the difference

eqStat(stat::StudentT_1S, direction::TestDir, design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
    Sum()
# ----- #


# ----- #

"""
```julia
function allPerms(stat::BivStatistic, ns::Int) 

function allPerms(stat::IndSampStatistic, ns::Vector{Int})

function allPerms(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int})

function allPerms(stat::OneSampStatistic, ns::Int)

```
Total number of systemetic permutations (exact test) for a test statistic given as singleton `stat` 
of a given [group](@ref "Statistic groups").

For the `ns` argument see [ns](@ref).

For an exact tests with ``N`` observations and ``K`` groups/measuremets, 
the total number of systematic permutations is

 - for `stat` a `BivStatistic` :``\\quad N!``  
 - for `stat` a `IndSampStatistic` : ``\\quad \\frac{N!}{N_1! \\cdot \\ldots \\cdot N_K!}``
 - for `stat` a `RepMeasStatistic` : ``\\quad K!^N`` 
 - for `stat` a `OneSampStatistic` : ``\\quad 2^N``

For the number of non-redundant permutations, which is the actual number of permutations that are listed
for exact tests, see [`nrPerms`](@ref).

*Examples*
```julia
# Total number of possible permutations for a 1-way ANOVA for indepedent samples 
# test with a total of 18 observations in three balanced groups.
allPerms(AnovaF_IS(), [6, 6, 6]) # -> 17153136

# Total number of possible permutations for a correlation test with 12 observations.
allPerms(PearsonR(), 12) # -> 479001600
```
"""
allPerms(stat::BivStatistic, ns::Int) = factorial(big(ns))

allPerms(stat::IndSampStatistic, ns::Vector{Int})  = multinomial(ns...)

allPerms(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int})  = factorial(big(ns.k))^ns.n

allPerms(stat::OneSampStatistic, ns::Int)  = exp2(BigInt(ns))
# ----- #

# ----- #
"""
```julia
function nrPerms(stat::Union{BivStatistic, OneSampStatistic}, ns::Int, total, 
                direction::TestDir, design::Assign=Balanced()) 

function nrPerms(stat::IndSampStatistic, ns::Vector{Int}, total, 
                direction::TestDir, design::Assign)

function nrPerms(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int}, total, 
                direction::TestDir, design::Assign=Balanced())
                
    where {TestDir<: TestDirection, Assign <: Assignment}

```
Number of non-redundant systemetic permutations. This is the actual number of permutations 
that are listed for an exact test.

Depending on the test direction and design, there may exists redundant permutations, *i.e.*, 
permutations that yields the same test statistic. 
These permutations are therefore eliminated, yielding a faster test.

`stat` is the test statistic used by the test given as a singleton. 
It belongs to one [group](@ref "Statistic groups") 
of test statistics.

For the `ns` argument see [ns](@ref).

`total` is the total number of systematic permutations. In order to avoid errors, it should be given 
as the resut of the [`allPerms`](@ref) function.

`direction` is a singleton of type [TestDirection](@ref).

`design` is a singleton of type [Assignment](@ref). In order to avoid errors, 
use the result of [`assignment`](@ref) as argument `design`. 

With `stat` a `IndSampStatistic` and a balanced design, the non-redundant permutations are ``\\frac{total}{K!}``, 
with ``K`` the number of groups.

With `stat` a `RepMeasStatistic` and a bi-directional tests, the non-redundant permutations are ``\\frac{total}{K!}``,
with ``K`` the number of measures.

For `stat` a `BivStatistic` or a `OneSampStatistic` return `total`, that is, there is no possible redundancy.

Note that the elimination of redundant permutations is rarely implemented in software packages for
permutation tests.

*Examples*
```julia
using PermutationTests
# Number of possible permutations for a 1-way ANOVA for indepedent samples 
# test with a total of 18 observations in three balanced groups.
ns=[6, 6, 6]
total=allPerms(AnovaF_IS(), ns) # -> 17153136
# Number of non-redundant permutations that will be actually listed
# to perform the test.
design=assignment(AnovaF_IS(), ns)
nrPerms(AnovaF_IS(), ns, total, Both(), design) # -> 2858856
```
"""
nrPerms(stat::Union{BivStatistic, OneSampStatistic}, ns::Int, total, direction::TestDir, 
            design::Assign=Balanced()) where {TestDir<: TestDirection, Assign <: Assignment} = total

nrPerms(stat::IndSampStatistic, ns::Vector{Int}, total, direction::TestDir, 
            design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} =
                direction isa Both && design isa Balanced ? total ÷ (factorial(big(length(ns)))) : total

nrPerms(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int}, total, direction::TestDir, 
            design::Assign=Balanced()) where {TestDir<: TestDirection, Assign <: Assignment} =
                direction isa Both ? total ÷ (factorial(big(ns.k))) : total

# ----- #


# ----- #

"""
```julia
function membership(stat::IndSampStatistic, ns::Vector{Int})
```
Create the appropriate argument `x` to be used by functions [`_permTest!`](@ref) and [`_permMcTest!`](@ref) when you 
[create your own test](@ref "Create your own test") using the permutation scheme of test statistics belonging 
to the `IndSampStatistic` [group](@ref "Statistic groups").

For `stat` a `IndSampStatistic`, `ns` is a group numerosity vector, *i.e.*, a vector of positive integers 
`[N1,...,NK]`, where ``K`` is the number of groups and ``N_k`` is the number of observations for the ``k^{th}`` group (see [ns](@ref)).  

Return the group membership vector `[repeat (1, N1);...; repeat(K, Nk)]`

If `rev=reverse` is passed as keyword argument, return instead
group membership vector `[repeat (K, N1);...; repeat(1, Nk)]`.
This is used to run t-tests for independent samples using the `PearsonR()`
[Statistic](@ref), see for example [`studentTestIS`](@ref).

*Examples*
```julia
using PermutationsTests
membership(AnovaF_IS(), [3, 4])
# return the vector [1, 1, 1, 2, 2, 2, 2]
```
"""
membership(stat::IndSampStatistic, ns::Vector{Int}; rev=identity) = 
    vcat([ones(Int, rev(ns)[i])*i for i ∈ rev(eachindex(ns))]...) # it was for i=1:length(ns)

"""
```julia
function membership(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int})
```
Create the appropriate argument `x` to be used by functions [`_permTest!`](@ref) and [`_permMcTest!`](@ref) when you 
[create your own test](@ref "Create your own test") using the permutation scheme of test statistics belonging 
to the `RepMeasStatistic` [group](@ref "Statistic groups").

For `stat` a `RepMeasStatistic`, `ns` is a named tuple, such as `(n=N, k=K)`, where ``N`` is the number of observations (e.g., *subjects*)
and ``K`` the number of measurements (or *treatments*, *times*, ect.), see [ns](@ref).

Return `collect(1:N*K)`.

*Examples*
```julia
using PermutationsTests
membership(AnovaF_RM(), (n=2, k=4))
# return the vector [1, 2, 3, 4, 5, 6, 7, 8]
```
"""
membership(stat::RepMeasStatistic, ns::@NamedTuple{n::Int, k::Int}) = collect(1:ns.n*ns.k) #repeat(1:ns.k, ns.n)

"""
```julia
function membership(stat::OneSampStatistic, ns::Int)
```
Create the appropriate argument `x` to be used by functions [`_permTest!`](@ref) and [`_permMcTest!`](@ref) when you 
[create your own test](@ref "Create your own test") using the permutation scheme of test statistics belonging 
to the `OneSampStatistic` [group](@ref "Statistic groups").

For `stat` a `OneSampStatistic`, `ns` is the number of observations (e.,g., *subjects*) given as an integer,
see [ns](@ref).

Return `ones(Int, ns)`.

*Examples*
```julia
using PermutationsTests
membership(Sum(), 5)
# return the vector [1, 1, 1, 1, 1]
```
"""
membership(stat::OneSampStatistic, ns::Int) = ones(Int, ns)
# ----- #


# ----- #
"""
```julia
function genPerms(stat::BivStatistic, x::UniData, 
            ns::Int, direction::TestDir, design::Assign) 
    
function genPerms(stat::IndSampStatistic, x::UniData, 
            ns::Vector{Int}, direction::TestDir, design::Assign)
                
function genPerms(stat::RepMeasStatistic, x::UniData, 
            ns::@NamedTuple{n::Int, k::Int}, direction::TestDir, design::Assign)
            
function genPerms(stat::OneSampStatistic, x::UniData, 
            ns::Int, direction::TestDir, design::Assign)

    where {TestDir<: TestDirection, Assign <: Assignment}
```

Generate a *lazy* iterator over all possible (systematic) permutations according to the permutation scheme 
to be used for test statistic `stat`, which is given as a singleton of type [Statistic](@ref).

Note that data permutation in *PermutationsTests.jl* are always obtained by lazy iterators, that is, permutations are never
listed physically. This is the main reason why the package is fast.

You do not need these functions for general usage of the package, however you need to know them 
if you wish to [create your own test](@ref "Create your own test"). 

There exists a permutation scheme for each [group](@ref "Statistic groups") the test statistics `stat` belong to.

`x` is the [`membership`](@ref) vector.

For the `ns` argument see [ns](@ref).

`direction` is a singleton of type [TestDirection](@ref).

`design` is a singleton of type [Assignment](@ref). In order to avoid errors, 
use the result of [`assignment`](@ref) as argument `design`. 

With `stat` a `BivStatistic` (*e.g.*, `PearsonR()`), the iterations unfold the ``N!`` possible reorderings of the 
``N`` elements in vector `x` and argument `ns` is ignored. In this case `x` is either a trend (e.g., [1,...,N] 
for a linear trend) or a variable to be permuted for correlation/trend tests.

With `stat` a `IndSampStatistic` (*e.g.*, `AnovaF_IS()`), the iterations unfold the 
``\\frac{N!}{N_1 \\cdot \\ldots \\cdot N_K}`` possible arrangements of ``N`` elements in ``K`` groups. 
In this case `x` is ignored and `ns` is a vector holding the K group 
numerosities (i.e., [N1,...,Nk]). For this scheme, if the design is balanced, *i.e.*, all elements of ns are equal, 
some permutations are redundant, see [`nrPerms`](@ref).

With `stat` a `RepMeasStatistic` (*e.g.*, `AnovaF_RM()`), the iterations unfold the ``(K!)^N`` reordering of ``K`` measures 
(e.g., treatments) in all ``N`` observation (e.g., subjects). In this case `x` is ignored and `ns` is a 
[named tuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) with form `(n=N, k=K)`.
For this scheme, if the test is bi-directional some permutations are redundant, see [`nrPerms`](@ref).

With `stat` a `OneSampStatistic` (*e.g.*, `StudentT_1S()`), the iteartions unfold the ``2^N`` flip-sign patterns. 
In this case `x` is ignored and `ns=N`, where ``N`` is the number of observations.

!!! note "nota bene"
    Only for `stat` belonging to the `OneSampStatistic` group, 
    the iterator generates tuples and not arrays (see examples below).

!!! tip "keep in mind"
    Regardless the permutation scheme, the first generated iteration always correspons to the "permutation"
    of the data as it has been observed (that is, the only permutation that actually does not permute the data).


*Examples*
```julia
using PermutationsTests

x=membership(PearsonR(), 3) # -> [1, 2, 3]
# The observations are always listed in the natural order 1, 2,...
Pr = genPerms(PearsonR(), x, 0, Both(), Balanced()) 
collect(Pr) # physically list the iterations
# yields the 6 = 3! elements [1, 2, 3]...[3, 2, 1]
# Here the integers represent the permuted position of the obervations in x

PfIS = genPerms(AnovaF_IS(), Int[], [2, 3], Both(), Unbalanced()) 
collect(PfIS)
# yields the 10 = 5!/2!3! elements [1, 1, 2, 2, 2]...[2, 2, 2, 1, 1]
# Here the integers represent the permuted group to which the corresponding 
# obervations in the data vector `y` belongs. see `_permTest!`.
# If the design is balanced and the test is bi-directional, 
# only 1/k! of the permutations are retained, thus 
PfIS_ = genPerms(AnovaF_IS(), Int[], [2, 2], Both(), Balanced()) 
collect(PfIS_) 
# yields the 3 = (4!/2!2!)/2! elements [1, 1, 2, 2], [1, 2, 2, 1] and [2, 1, 2, 1]

PfRM = genPerms(AnovaF_RM(), Int[], (n=2, k=3), Right(), Balanced()) 
collect(PfRM)
# yields the 36 = (k!)^n elements [1, 2, 3, 4, 5, 6], [1, 3, 2, 4, 5, 6], 
# ..., [3, 2, 1, 6, 5, 4]
# Here the integers represent the treatement for each subject to which the 
# corresponding obervation in the data vector `y` belongs. see `_permTest!`.
# If the test is bi-directional, 
# only 1/k! of the permutations are retained, thus 
PfRM_ = genPerms(AnovaF_RM(), Int[], (n=2, k=3), Both(), Balanced()) 
collect(PfRM_)
# yields instead the 6 = (k!)^n / elements [1, 2, 3, 4, 5, 6], [1, 3, 2, 4, 5, 6], 
# ..., [3, 2, 1, 4, 5, 6]

Pt1S = genPerms(StudentT_1S(), Int[], 3, Both(), Balanced())
collect(Pt1S)
# yields the 8=2^3 elements (1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1), 
#..., (-1, -1, -1) 
# Here the integers represent the sign to be applied to each observation 
# at each data permutation. 
# Note that only for OneSampStatistic, the iterator generates tuples and not arrays.
```
"""
genPerms(stat::BivStatistic, 𝐱::UniData, ns::Int, direction::TestDir, 
            design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} = permutations(𝐱)

function genPerms(stat::IndSampStatistic, 𝐱::UniData, ns::Vector{Int}, direction::TestDir, 
                design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} 
    n, k, msp = sum(ns), length(ns), multiset_permutations
    # for balanced design in IS ANOVA-like tests, only the permutations respecting condition i % factorial(k) == 1 are needed
    memb = vcat([ones(Int, ns[i])*i for i ∈ eachindex(ns)]...)
    return design isa Balanced && direction isa Both ? (p for (i, p) in enumerate(msp(memb, n)) if i % factorial(k) == 1) : msp(memb, n)
end

# the concatenation of iterators' product of (permutations(1:k), permutations(1k+1:1k+k), ..., permutations((n-1)k+1:(n-1)*k+k))
function genPerms(stat::RepMeasStatistic, 𝐱::UniData, ns::@NamedTuple{n::Int, k::Int}, direction::TestDir, 
                design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} 
    # for bi-directional tests RM ANOVA-like tests, only the first non_redundant = all ÷ factorial(k) permutations are needed
    if !(direction isa Both) 
        return (vcat(p...) for p in Iterators.product((permutations((i*ns.k)+1:(i*ns.k)+ns.k ) for i in 0:ns.n-1)...)) # yields : [1, 2, 3, 4, 5, 6], [2, 1, 3, 4, 5, 6],...   
    else
        needed = allPerms(stat, ns) ÷ factorial(ns.k) # take only the first `needed` permutations :
        return Iterators.take((vcat(p...) for p in Iterators.product((permutations((i*ns.k)+1:(i*ns.k)+ns.k ) for i in 0:ns.n-1)...)), needed) 
    end
end

genPerms(stat::OneSampStatistic, 𝐱::UniData, ns::Int, direction::TestDir, 
            design::Assign) where {TestDir<: TestDirection, Assign <: Assignment} = 
    Iterators.product(((1., -1.) for i=1:ns)...)
# Iterators.map(SVector, Iterators.product(((1., -1.) for i=1:ns)...)) !!!!!! xxx Arrays
# ----- #


# ----- #
# check the `ns` argument, which is given as input for several test functions
function _check_ns(𝐲, ns, stat::Union{BivStatistic, OneSampStatistic}) 
    ns isa Int || throw(ArgumentError(📌*" Function _check_ns: for bivariate and OneSampStatistic statistics `ns` must be an integer"))
    ns≠length(𝐲) && throw(ArgumentError(📌*" Function _check_ns: for bivariate and OneSampStatistic statistics `ns` must be an integer equal to the length of 𝐲"))
end

function _check_ns(𝐲, ns, stat::IndSampStatistic) 
    ns isa IntVec || throw(ArgumentError(📌*" Function _check_ns: for Independent Samples statistics `ns` must be a vector of integer"))
    sum(ns)≠length(𝐲) && throw(ArgumentError(📌*" Function _check_ns: for Independent Samples statistics the elements in `ns` must sum up to the length of vector argument 𝐲"))
end

function _check_ns(𝐲, ns, stat::RepMeasStatistic) 
    ns isa NamedTuple || throw(ArgumentError(📌*" Function _check_ns: for Repeated Measure statistics `ns` must be a NamedTuple"))
    ns.n*ns.k≠length(𝐲) && throw(ArgumentError(📌*" Function _check_ns: for Repeated Measure statistics the `k`(# of treatments) and `n`(# of subjects) fields in `ns` must sum up to the length of vector argument 𝐲"))
end
# ----- #

# ----- #
# Return the test type (:exact or :approximate) and the number of permutations, given the number of observations,
# the group numerosity vector n, the test statistic `stat`, the default number of permutations 
# for an approximate test `nperm` and the upper limit for the number of permutations to allow an exact test.
# If observations > 30 or the number of systematic permutations exceeds `switch2rand` the the approximate 
# test with nperm permutations is chosen, otherwise the exact test is chosen.
function _test_params(stat::Stat, ns, fstat::Function, observations, nperm, switch2rand) where Stat<:Statistic
    design = assignment(stat, ns) # Balanced() or Unbalanced()
    direction = _direction(fstat)
    if observations > 26 # avoid to compute large factorials
        testtype = :approximate
        nonRed = nperm
    else
        totalperm = allPerms(stat, ns) # Total systematic permutations.
        nonRed = BigInt(nrPerms(stat, ns, totalperm, direction, design)) # non-redundant systematic permutations
        if nonRed>switch2rand
            testtype = :approximate
            nonRed = nperm
        else
            nperm=Int(nonRed) # (totalperm)  # Typecast to Int as allPerms uses BigInt. NB thi change argument nperm
            testtype = :exact
        end
    end

    return testtype, nperm, direction, design
end
# ----- #


# ----- #
# make checks and prepare the ensuing _perm_test! function
function _prepare_permtest!(𝐱, 𝐲::UniData, ns, stat, fstat, nperm, switch2rand, seed, standardized, centered)#, means, sds)

    # check ns argument 
    _check_ns(𝐲, ns, stat)
    
    # determine whether the permutations are to be systematic or random (testtype), their number (npern) ,
    # the test direction (direction) and whether the design is balanced or unbalanced (design)
    testtype, nperm, direction, design  = _test_params(stat, ns, fstat, length(𝐲), nperm, switch2rand)

    # for OneSampStatistic and approximate tests 𝐱 is not used. For all other cases 𝐱 and 𝐲 must have equal length)
    !((stat isa OneSampStatistic) & (testtype == :approximate)) && length(𝐱) ≠ length(𝐲) && throw(ArgumentError(📌*"Function _prepare_permtest!: the length of vector argument 𝐱 and 𝐲 is not equal"))
  
    rng = nothing
    testtype == :approximate && (rng = MersenneTwister(seed==0 ? rand(UInt32) : seed))
    # use Random.seed!(1234) to reset the seed of the random number generator;

    getkwargs(𝐱, 𝐲, ns, stat::CrossProd)    = (k=0, ns=0)
    getkwargs(𝐱, 𝐲, ns, stat::Covariance)   = (k=0, ns=0, centered=centered)
    getkwargs(𝐱, 𝐲, ns, stat::PearsonR)     = (k=0, ns=0, standardized=standardized, centered=centered)#, means=means, sds=sds)
    getkwargs(𝐱, 𝐲, ns, stat::IndSampStatistic) = (_anova_IS_fixed_params(𝐱; askwargs=true)..., ft=0) # ft is dummy, but otherwise does not compile
    getkwargs(𝐱, 𝐲, ns, stat::Union{RepMeasStatistic, OneSampStatistic}) = (k=0, ns=ns)
    kwargs = getkwargs(𝐱, 𝐲, ns, stat)
    
    return testtype, nperm, direction, design, kwargs, rng
end


# ----- #

"""
```julia
function table2vec(table::Matrix{I}, stat::IndSampStatistic) 

function table2vec(tables::AbstractVector{Matrix{I}}, stat::IndSampStatistic) 

function table2vec(table::Matrix{I}, stat::RepMeasStatistic) 

function table2vec(tables::AbstractVector{Matrix{I}}, stat::RepMeasStatistic) 

    where I<:Int

```
Format tables of dicothomous data so as to create data input for test functions. 

Return the 2-tuple holding the formatted table as a vector (for univariate test) or
a vector of vectors (for multiple comparisons tests) and the appropriate argument [ns](@ref).

You do not need to use these functions in general, as they are called internally by the test functions.
However they may be useful if you [create your own test](@ref "Create your own test") for dicothomous data.

`stat` is a singleton of the [Statistic](@ref) type, that is, a test-statistic belonging to one
of the [group of test statistics](@ref "Statistic groups").

For the usage of these functions, see the documentation of the test functions using them:

 **Univariate test functions**

 - [`chiSquaredTest`](@ref)
 - [`fisherExactTest`](@ref)
 - [`cochranqTest`](@ref)
 - [`mcNemarTest`](@ref)

 **Multiple comparisons test functions**
 
 - [`chiSquaredMcTest`](@ref)
 - [`fisherExactMcTest`](@ref)
 - [`cochranqMcTest`](@ref)
 - [`mcNemarMcTest`](@ref)
"""
function table2vec(table::Matrix{I}, stat::IndSampStatistic) where I<:Int
    size(table, 1) > 2 && throw(ArgumentError(📌*" Function table2vec: the contingency table may have at most two rows (but any K≥2 columns)"))
    return Float64.(vcat([fill(1-(r-1), table[r, c]) for c in axes(table, 2) for r in axes(table, 1)]...)), vec(sum(table, dims=1))
end

function table2vec(tables::AbstractVector{Matrix{I}}, stat::IndSampStatistic) where I<:Int
    size(tables[1], 1) > 2 && throw(ArgumentError(📌*" Function table2vec: the contingency tables may have at most two rows (but any K≥2 columns)"))
    return [Float64.(vcat([fill(1-(r-1), t[r, c]) for c in axes(t, 2) for r in axes(t, 1)]...)) for t ∈ tables], vec(sum(tables[1], dims=1))
end


function table2vec(table::Matrix{I}, stat::RepMeasStatistic) where I<:Int
    N, K = size(table)
    K < 2 && throw(ArgumentError(📌*" Function table2vec: the table must have at least two columns (treatements)"))
    N < K && throw(ArgumentError(📌*" Function table2vec: the number of rows of the table (subjects) must be graeter than the number of columns (treatments)"))
    return Float64.(vcat(table'...)), (n=N, k=K)
end


function table2vec(tables::AbstractVector{Matrix{I}}, stat::RepMeasStatistic) where I<:Int
    N, K = size(tables[1])
    K < 2 && throw(ArgumentError(📌*" Function table2vec: the tables must have at least two columns (treatements)"))
    N < K && throw(ArgumentError(📌*" Function table2vec: the number of rows of the tables (subjects) must be graeter than the number of columns (treatments)"))
    return [Float64.(vcat(t'...)) for t ∈ tables], (n=N, k=K)
end

# ----- #
