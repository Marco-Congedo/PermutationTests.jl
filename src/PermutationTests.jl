#=
Main Module of the PermutationTests.jl Package
v0.1.0

MIT License
Copyright (c) 2024,
Marco Congedo, CNRS, Grenoble, France:
https://sites.google.com/site/marcocongedo/home
=#

module PermutationTests

using Base.Threads
using Statistics: mean, var
using Random: rand, shuffle, shuffle!, MersenneTwister
using Combinatorics: permutations, multiset_permutations, multinomial
using Folds: Folds, Folds.maximum

# import

export
# from this module
Statistic,
BivStatistic,
    PearsonR, 
        CrossProd, Covariance,
IndSampStatistic,
    AnovaF_IS, 
        SumGroupTotalsSq_IS, SumGroupTotalsSqN_IS, StudentT_IS, Group1Total_IS,
RepMeasStatistic,
    AnovaF_RM, 
        SumTreatTotalsSq_RM, StudentT_RM,
OneSampStatistic,
    StudentT_1S, 
        Sum,
TestDirection,
Right,
Left,
Both,
Assignment,
Balanced,
Unbalanced,
TestResult,
UniTest,
MultcompTest,

# from stats.jl
∑,
∑of²,
∑∑of²,
μ,
dispersion,
σ²,
σ,
∑,
Π,
∑ofΠ,
statistic,
_∑y²,
_∑Y²kn_∑y²_∑S²k,
μ0,
σ1,
μ0σ1,
# from tools.jl
flip,
assignment,
eqStat,
allPerms,
nrPerms,
membership,
genPerms,
table2vec,
# from uniTests.jl
_permTest!,
# from uniTests_API.jl
correlationTest!, rTest!,
correlationTest,  rTest,
trendTest!, trendTest,
pointBiSerialTest,
anovaTestIS, fTestIS,
studentTestIS, tTestIS, 
chiSquaredTest, Χ²Test, 
fisherExactTest,
anovaTestRM, fTestRM,
cochranqTest, qTest,
mcNemarTest,
studentTestRM, tTestRM,
studentTestRM!, tTestRM!,
studentTest1S!, tTest1S!,
studentTest1S, tTest1S,
signTest!, signTest,
# from multcompTests.jl
_observedStats,
_permutedStat,
_permMcTest!,
# from multcompTests_API.jl
correlationMcTest!, rMcTest!,
correlationMcTest,  rMcTest,
trendMcTest!, trendMcTest,
pointBiSerialMcTest,
anovaMcTestIS, fMcTestIS,
studentMcTestIS, tMcTestIS, 
chiSquaredMcTest, Χ²McTest, 
fisherExactMcTest,
anovaMcTestRM, fMcTestRM,
cochranqMcTest, qMcTest,
mcNemarMcTest,
studentMcTestRM, tMcTestRM,
studentMcTestRM!, tMcTestRM!,
studentMcTest1S!, tMcTest1S!,
studentMcTest1S, tMcTest1S,
signMcTest!, signMcTest


# Consts
const titleFont     = "\x1b[38;5;208m" # orange
const diceFont      = "\x1b[38;5;210m" # orange
const separatorFont = "\x1b[38;5;223m" # orange clair
const defaultFont   = "\x1b[0m"
const greyFont      = "\x1b[90m"
const 📌            = titleFont*"PermutationTests.jl "*defaultFont
const dice = ("⚀", "⚁", "⚂", "⚃", "⚄", "⚅")

# Types
# Types for input data: All Real, Integers and Boolean (Bits)
DataType = Union{R, I, Bool} where {R<:Real, I<:Int}
UniData = Union{AbstractVector{R}, AbstractVector{I}, AbstractVector{Bool}} where {R<:Real, I<:Int}
UniDataVec = Union{AbstractVector{Vector{R}}, AbstractVector{Vector{I}}, AbstractVector{Vector{Bool}}} where {R<:Real, I<:Int}
UniDataVec² = Union{AbstractVector{Vector{Vector{R}}}, AbstractVector{Vector{Vector{I}}}, AbstractVector{Vector{Vector{Bool}}}} where {R<:Real, I<:Int}

# Useful types
IntVec = AbstractVector{I} where I<: Int
Into = Union{I, Nothing} where I <: Int
Realo = Union{R, Nothing} where R <: Real
nsType=Union{Int, Vector{Int}, @NamedTuple{n::Int, k::Int}}

# Abstract Type which children are all statistics
abstract type Statistic end 

# singletons (statistics) and Unions of singletons (statistics and their equivalent statistics) using the same permutation scheme :
# _________________________________________________________________________________________________________________________________
struct PearsonR             <: Statistic end # PRINCIPAL: Pearson product moment correlation coefficient (H0: R≠0)
struct CrossProd            <: Statistic end # EQUIVALENT to PearsonR for the case of directional test
struct Covariance           <: Statistic end # EQUIVALENT to PearsonR for the case of bi-directional test
BivStatistic = Union{Covariance, PearsonR, CrossProd} ### Group of Bivariate Statistics
# __________________________________________
struct AnovaF_IS            <: Statistic end # PRINCIPAL: 1-way ANOVA for independent samples (H0: µ1 = µ2 = µ3 = … = µk)
                                             # With dichotomous data yields 1-saided and bi-directional Χ² test and Fisher Exact Tests for >2 groups
struct SumGroupTotalsSq_IS  <: Statistic end # EQUIVALENT to AnovaF_IS for balanced designs (groups of equal numerosity)
struct SumGroupTotalsSqN_IS <: Statistic end # EQUIVALENT to AnovaF_IS in general and to StudentT_IS in the case of bi-directional tests
struct StudentT_IS          <: Statistic end # PRINCIPAL: Student T for independent samples (H0: μ1≠μ2)
                                             # With dichotomous data yields 1-saided and bi-directional Χ² test and Fisher Exact Tests for 2 groups
                                             # For ranked data gives the same p-value as the Mann-Whitney U test
struct Group1Total_IS       <: Statistic end # EQUIVALENT to StudentT_IS for the case of directional tests μ(G1) > μ(G2)
IndSampStatistic = Union{AnovaF_IS, SumGroupTotalsSq_IS, SumGroupTotalsSqN_IS, StudentT_IS, Group1Total_IS} ### Group of Ind Samp Statistics
# __________________________________________
struct AnovaF_RM            <: Statistic end # PRINCIPAL: 1-way repeated-measure ANOVA (H0: µ1 = µ2 = µ3 = … = µk) 
                                             # With dichotomous data yields 1-saided and bi-directional Cochran Q test (K>2) and to the McNemar and sign test (K=2) 
                                             # For ranked data gives the same p-value as the Friedman test
struct SumTreatTotalsSq_RM  <: Statistic end # EQUIVALENT to AnovaF_RM in all cases
struct StudentT_RM          <: Statistic end # PRINCIPAL: Student T for paired samples (repeated-measure) (H0: μ1≠μ2) 
                                             # EQUIVALENT data gives the same p-value as the Wilcoxon signed rank test
RepMeasStatistic = Union{AnovaF_RM, SumTreatTotalsSq_RM, StudentT_RM} ### Group of Repeated-Measure Statistics
# __________________________________________
struct StudentT_1S          <: Statistic end # PRINCIPAL: Student T for one sample (H0: μ≠0)
struct Sum                  <: Statistic end # EQUIVALENT to Student T for one sample (H0: μ≠0)
OneSampStatistic = Union{StudentT_1S, Sum}   ### Group of One-sample statistics
# _________________________________________________________________________________________________________________________________

AllStatistics=Union{BivStatistic, IndSampStatistic, RepMeasStatistic, OneSampStatistic}

# Parameter-free statistics, those that do not need precomputed data to be computed faster 
ParFreeStatistics=Union{BivStatistic, IndSampStatistic, SumTreatTotalsSq_RM, StudentT_RM, Sum}

abstract type TestDirection end # Right-, Left-Directional or Bi-Directional
# singletons :
struct Right        <: TestDirection end # 
struct Left         <: TestDirection end # 
struct Both         <: TestDirection end # 

abstract type Assignment end # Balanced or Unbalanced
struct Balanced     <: Assignment end # 
struct Unbalanced   <: Assignment end # 

abstract type TestResult end

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

struct MultcompTest <: TestResult
    p::Vector{Float64} # different from UniTest
    stat::statistic where statistic <: Statistic
    obsstat::Vector{Float64} # different from UniTest
    minp::Float64 
    nperm::Int64 
    testtype::Symbol 
    direction::testDirection where testDirection <: TestDirection
    design::assignment where assignment <: Assignment
    nulldistr::Vector{Float64} # in addition to UniTest
    rejections::Vector{Vector{Int64}} # in addition to UniTest
    stepdown::Bool # in addition to UniTest
    fwe::Float64 # in addition to UniTest
end    


include("stats.jl")
include("tools.jl")
include("uniTests.jl")
include("uniTests_API.jl")
include("multcompTests.jl")
include("multcompTests_API.jl")


# Override Base.Show for UniTest and MultcompTest
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, t::Union{UniTest, MultcompTest})
    println(io, titleFont, t isa UniTest ? ("\n"*rand(dice)*" Univariate Permutation test") : 
                                           ("\n"*reduce(*, [rand(dice)*" " for i=1:3])*"... Multiple Comparison Permutation test"))
    if t isa UniTest
        println(io, defaultFont, "p-value = ", t.p < 0.0001 ? "<0.001" : round(t.p; digits=3)) 
    else
        steps=max(1, length(t.rejections)) # number of steps
        stepsstr=steps==1 ? " step" : " steps"
        nrej=length(vcat(t.rejections...)) # number of rejections
        nhyp=length(t.p) # number of hypotheses
        proprej=round((nrej/nhyp)*100, digits=2) # proportion of rejections
        println(io, defaultFont, "Rejected ", nrej, " out of ", nhyp, " hypotheses ", "(", proprej, " %) in ", steps, stepsstr, " with FWE=", round(t.fwe, digits=8))
    end
    print(io, greyFont, " ")
    println(io, reduce(*, [rand(dice)*"      " for i=1:10]))

    print(io, separatorFont,".p ", greyFont, t isa UniTest ? "(p-value) " : "(p-values)")
    print(" ") 
    print(io, separatorFont,".stat ", greyFont, "(test statistic)")
    print("   ")
    println(io, separatorFont,".obsstat ", greyFont, t isa UniTest ? "  (observed statistic)" : " (observed statistics)")
    print(io, separatorFont,".minp ", greyFont, " (minimum attainable p-value)")
    print("   ")
    println(io, separatorFont,".nperm ", greyFont, "(number of permutations)")
    print(io, separatorFont,".testtype ", greyFont, " (exact or approximated)")  
    print("  ")
    println(io, separatorFont,".direction ", greyFont, "(Both, Left or Right)")
    println(io, separatorFont,".design ", greyFont, " (Balanced or Unbalanced)")
    if t isa MultcompTest
        print(io, separatorFont,".nulldistr ", greyFont, "  (null distribution)")
        print("  ")
        println(io, separatorFont,".rejections ", greyFont, "(for each step down)")
        print(io, separatorFont,".stepdown ", greyFont, "   (true if stepdown)")
        print("  ")
        println(io, separatorFont,".fwe ", greyFont, "(family-wise err. stepdown)")
    end
end # show
  
println("\n⭐ "," Welcome to the ", 📌, "package", " ⭐\n")
@info " "
println(" Your Machine `", separatorFont, gethostname(), defaultFont, "` (",Sys.MACHINE, ")")
println(" runs on kernel ",Sys.KERNEL, " with word size ", Sys.WORD_SIZE,".")
println(" CPU  Threads: ", separatorFont, Sys.CPU_THREADS, defaultFont)
println(" Base.Threads: ", separatorFont, Threads.nthreads(), defaultFont)

end # module
