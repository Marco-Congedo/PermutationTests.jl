#   Unit "_test.jl", additional resource to the PemutationTests.jl package for the julia language
#
#   MIT License
#   Copyright (c) 2024,
#   Marco Congedo, CNRS, Grenoble, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit tests the validity of all univariate and multiple comaprison tests implemented
#   in the PermutationsTests.jl package. 

#   The general methodology is to create data under H0 and perform the test.
#   The procedure is repeated many times in order to estimate the rejection threshold
#   at many type I error rates (univariate tests) and the average proportion of rejections 
#   for a given family-wise error rate (multiple comparison tests).

#   REQUIRED PACKAGES : PermutationsTests, Plots, Distributions

using Dates:now, toms, Minute # julia built-in
using Random:MersenneTwister # julia built-in
using Plots, Plots.Measures
using Distributions: Binomial
using PermutationTests
plotly()

# Grid of histograms
function Hplot(P, N, time, filename; noplot=false, layout=(2, 3))
    plots=(Plots.stephist(P[:, i], xlims=(0, 1), label=false, xlabel="p-values", ylabel="Frequency", title="N=$(N[i])" , 
            titlefont=10, bin=40, bottom_margin = 10mm, left_margin = 10mm) for i=1:size(P, 2))
    p=Plots.plot(plots..., size=(1000, 600), layout = layout)
    #png(filename)
    noplot || return p
end

# Grid of p-p plots
function ppplot(x, P, N, time, filename; noplot=false, layout=(2, 2))
    plots=(Plots.plot([x, sort(P[:, i])], x, xlims=(0, 1), ylims=(0, 1), label=false, #label=["Expected" "Observed"], 
    xlabel="Expected p-value", ylabel="Observed p-value", 
        title="N=$(N[i])" , titlefont=10, bottom_margin = 10mm, left_margin = 10mm) for i=1:size(P, 2))
    p=Plots.plot(plots..., size=(1000, 600), layout = layout)
    #png(filename)
    noplot || return p
end


## UNIVARIATE TESTS

# set the monte carlo simulations
Nrep=1000
x=[i/Nrep for i=1:Nrep] # x-data for p-p plots
s2r=Int(1e5) # switch 2 random if more than 1e5 permutations


# Correlation test
N=[6, 12, 24, 48] # number of observations
begin
    ⌚ = now()
    P = [rTest(randn(n), randn(n); switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni r test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni r test_L"); layout=(2, 2))


# 1-way ANOVA for indepedent samples
N=[9, 18, 36, 72] # number of observations
begin
    ⌚ = now()
    P = [fTestIS(randn(n), [n÷3, n÷3, n-((n÷3)*2)]; switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni fIS test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni fIS test_L"); layout=(2, 2))


# t-test for indepedent samples
N=[9, 18, 36, 72] # number of observations
begin
    ⌚ = now()
    P = [tTestIS(randn(n), [n÷2, n-(n÷2)]; switch2rand=s2r, verbose=false, asPearson=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni tIS test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni tIS test_L"); layout=(2, 2))


# Chi-Square test
K=[2, 2, 4, 4] # number of columns in table
N=[10, 20, 10, 20]
rng = MersenneTwister(1234)
begin
    ⌚ = now()
    P = [Χ²Test(rand(rng, Binomial(n, 0.5), 2, k), switch2rand=Int(1E6), verbose=false).p for r=1:Nrep, (n, k) ∈ zip(N, K)]
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni Chi test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni Chi test_L"); layout=(2, 2))


# 1-way ANOVA for repeated measures
N=[9, 18, 36, 72] # number of observations
begin
    ⌚ = now()
    P = [fTestRM(randn(n), (n=n÷3, k=3); switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni fRM test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni fRM test_L"); layout=(2, 2))


# Cochran Q test
N=[9, 18, 36, 72] # number of observations
k=3
begin
    ⌚ = now()
    P = [qTest(Int.(rand(Bool, n÷k, k)); switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni q test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni q test_L"); layout=(2, 2))


# one sample t-test
N=[9, 18, 36, 72] # number of observations
begin
    ⌚ = now()
    P = [tTest1S(randn(n); switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2))
Hplot(P, N, time, joinpath(@__DIR__, "uni t1S test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni t1S test_L"); layout=(2, 2))


# Sign Test
N=[9, 18, 36, 72] # number of observations
begin
    ⌚ = now()
    P = [signTest(rand(Bool, n); switch2rand=s2r, verbose=false).p for r=1:Nrep, n ∈ N]   
    time = now()-⌚ 
end
times=string(round(toms(time) / toms(Minute(1)); digits=2)) 
Hplot(P, N, time, joinpath(@__DIR__, "uni sign test_H"); layout=(2, 2))
ppplot(x, P, N, time, joinpath(@__DIR__, "uni sign test_L"); layout=(2, 2))


# one-sided tests ...


############################################################################################"""""


## MULTPLE COMPARISON TESTS

# set up monte carlo simutaions

FWE=0.05 # Family-wise error
Nrep=100 # number of repetitions
s2r=Int(1e5) # switch 2 random if more than 1e5 permutations
M=[2, 5, 10, 100] # number of hypoteses
R=Matrix{Float64}(undef, length(M), 6)


# Correlation test
n=6
nrPerms(PearsonR(), n, allPerms(PearsonR(), n), Both())
for i=1:length(M)
    rejected = sum(count(≤(FWE), rMcTest(randn(n), [randn(n) for j=1:M[i]]; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 1] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 1].>FWE)>0 ? (@warn "the r test does not control the FWE" R[:, 1]) : 
                        (@info "the r test controls the FWE" R[:, 1])
# result: R = [ 0.03, 0.008, 0.002, 0.0007]


# ANOVA for independent samples
n=8
ns=[n÷3, n÷3, n-((n÷3)*2)]
nrPerms(AnovaF_IS(), ns, allPerms(AnovaF_IS(), ns), Both(), Unbalanced())
for i=1:length(M)
    rejected = sum(count(≤(FWE), fMcTestIS([randn(n) for j=1:M[i]], ns; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 2] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 2].>FWE)>0 ? (@warn "the ANOVA IS test does not control the FWE" R[:, 2]) : 
                        (@info "the ANOVA IS test controls the FWE" R[:, 2])
# result: R = [ 0.025, 0.008, 0.008, 0.0005]

# t-test for independent samples
n=12
ns=[n÷2, n-(n÷2)]
nrPerms(StudentT_IS(), ns, allPerms(StudentT_IS(), ns), Both(), Balanced())
for i=1:length(M)
    rejected = sum(count(≤(FWE), tMcTestIS([randn(n) for j=1:M[i]], ns; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 3] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 3].>FWE)>0 ? (@warn "the t-test IS does not control the FWE" R[:, 3]) : 
                        (@info "the t-test IS controls the FWE" R[:, 3])
# result: R = [ 0.03, 0.014, 0.003, 0.0001]

# Chi-Square
n=3 # expected frequency in each cell
k=3 # number of columns in tables
rng = MersenneTwister(1234)

# generate tables with fixed column sum equal to twice the expected frequency `n`
# the first element of each column is generate randomly and the second is such that the col sum is 2n
function gettable(n, k, rng)
    table=zeros(Int64, 2, k)
    for i=1:size(table, 2)
        table[1, i]=rand(rng, Binomial(n, 0.5))
        table[2, i]=(2*n)-table[1, i]
    end
    return table
end

for i=1:length(M)
    rejected = sum(count(≤(FWE), Χ²McTest([gettable(n, k, rng) for j=1:M[i]]; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 4] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 4].>FWE)>0 ? (@warn "the Chi-Suqared test does not control the FWE" R[:, 4]) : 
                        (@info "the Chi-Suqared test controls the FWE" R[:, 4])

# ANOVA for repeated measures
n=12
ns=(n=4, k=3)
nrPerms(AnovaF_RM(), ns, allPerms(AnovaF_RM(), ns), Both())
for i=1:length(M)
    rejected = sum(count(≤(FWE), fMcTestRM([randn(n) for j=1:M[i]], ns; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 5] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 5].>FWE)>0 ? (@warn "the ANOVA RM test does not control the FWE" R[:, 5]) : 
                        (@info "the ANOVA RM test controls the FWE" R[:, 5])
# result: R = [ 0.005, 0.01, 0.008, 0.0003]


# one-sample t-test
n=9
nrPerms(StudentT_1S(), n, allPerms(StudentT_1S(), n), Both())
for i=1:length(M)
    rejected = sum(count(≤(FWE), tMcTest1S([randn(n) for j=1:M[i]]; 
            switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, stepdown=true).p) for r=1:Nrep)
    println("done for m=", M[i])
    R[i, 6] = (rejected/(Nrep*M[i])) # proportion of rejected hypotheses
end
sum(R[:, 6].>FWE)>0 ? (@warn "the ANOVA RM test does not control the FWE" R[:, 6]) : 
                        (@info "the ANOVA RM test controls the FWE" R[:, 6])
# result: R = [  0.015, 0.008, 0.008, 0.0007]

kwargs=(minorgrid=true, label=["r" "F(is)" "t(is)" "Χ²" "F(rm)" "t"], xlabel="Number of hypotheses", 
        leg=:topright, lw=2, xticks = ([1:1length(M);], string.(M)),
        tickfontsize=11, labelfontsize=14, legendfontsize=10, right_margin = 10mm, size=(900, 600))


plot(R; ylabel="Average proportion of rejections", title="FWE control for α=0.05", ylims=(0, 0.05), kwargs...)


plot(R; ylabel="average proportion of rejections")
