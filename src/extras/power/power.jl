#   Unit "_power.jl", additional resource to the PemutationTests.jl package for the julia language
#
#   MIT License
#   Copyright (c) 2024,
#   Marco Congedo, CNRS, Grenoble, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit compares the power of all univariate and multiple comparison tests implemented
#   in the PermutationsTests.jl package to parametric equivalent tests. 

#   The general methodology is to create data under under H1 with some affect size (under H0 and H1 for multiple
#   comparison tests) and perform the tests.
#   The procedure is repeated many times in order to estimate the probability to reject the
#   hypothesis for the univariate tests and the proportion of rejected hypotheses among the false hypotheses 
#   for multiple comparison tests. It also includes some more in depth analysis of the behavior of the tests.

#   REQUIRED PACKAGES : PermutationsTests, Plots, StatsPlots, LaTeXStrings, Distributions, LinearAlgebra, PrettyTables

using Plots, Plots.Measures, StatsPlots, LaTeXStrings
using Distributions
using LinearAlgebra: I, Hermitian
using PrettyTables
using PermutationTests 

plotly() # select Plots.jl backend


##** UTILITIES FOR Correlation

# Create M multivariate Normal Standard variables of N elements with unit variance and correlation matrix such that
# 1) if M=2, the first variable has correlation ρ with the second
# 2) if M>2, the first variable is given correlation ρ with all the others 
# and variables 2 to M are given correlation ρ2 in between each other.
# The data is returned to be compatible as input of rTest and rMcTest,
# that is, as two vectors 𝐱, 𝐲 if M==2 and first vector x, vectors 2 to M as a vector of vector 𝐘 if M>2.
# Example: x, y = gen_correlated_data(10, 2, sigr(10, 0.05), 0) # M=2
# Example: x, y = gen_correlated_data(10, 20, sigr(10, 0.05), sigr(10, 0.1)) # M>2
function gen_correlated_data(N, M, ρ, ρ2)

    R=Matrix(1.0I, M, M) # deired correlation matrix
    for i=2:size(R, 2) # correlation between the first and all other variables. it suffice to write the upper off diagonal elements
        R[1, i]=ρ
    end
    if M>2
        for j=2:size(R, 1), i=j+1:size(R, 2) # correlation among variables 2:M. it suffice to write the upper off diagonal elements
            R[j, i]=ρ2
        end
    end
    D=MvNormal(fill(0.0, M), Matrix(Hermitian(R)))
    
    Y=Matrix(rand(D, N)) # generate random sample

    if M == 2
        return Y[1, :], Y[2, :]
    else
        return Y[1, :], [Y[i, :] for i=2:M]
    end
end

# Return the correlation of two variables of N elements must have so as to be significantly correlated 
# at the α level using the unidirectional t-test for the correlation coefficient: t = (r*(sqrt(n-2)))/sqrt(1-r²)
# Solving the above equation for r given t, yields r=sqrt(x/(x+1)), where x=t²/(n-2)
function sigr(N, α)
    x=abs2(cquantile(TDist(N-2), α))/(N-2) 
    return sqrt(x/(x+1))
end

# use to check the above equation : r2t(sigr(N, α), N) must give the t(N-2 df)-value yielding α under the curve on its right
r2t(r, N)=(r*(sqrt(N-2)))/sqrt(1-abs2(r))

# Fisher z-transformation r -> Normal(0, 1)
#Fisherz(r, N) = atanh(r)/sqrt((1/N)+6/(2*N^2))

# asymptotyc correlation test based on t-test with N-2 df
rTestAsy(𝐱, 𝐲, N) = ccdf(TDist(N-2), r2t(cor(𝐱, 𝐲), N))

# asymptotyc correlation test based on the Fisher transformation of correlation coefficient
#rTestAsy(𝐱, 𝐲, N) = ccdf(Normal(), Fisherz(cor(𝐱, 𝐲), N))

##**

##-- UTILITIES FOR T-TESTS INDEPENDENT Samples

# return M pairs of independent samples with N÷2 elements each.
# 1) if M=1 the mean of the second sample is shifted bu μ, i.e., if μ1 is the mean of the distribution, the second sample has mean μ1+μ.
# 2) if M>1 the mean of the second sample of the first nH1 pairs only shifted bu μ. By defalt nH1=M, i.e., all samples are shifted.
# Furthermore, if M>1, the expected correlation between the M variables within the first and second group is ρ. 
# The data is returned to be compatible as input of tTestIS and tMcTestIS,
# that is, as a vector 𝐲 if M==1 and a vector of M vectors 𝐲 if M>1, where 𝐲 concatenates the two samples.
# Example: Y = gen_shifted_data(10, 1, sigt(10, 0, 0.05), 0) # M=1
# Example: Y = gen_shifted_data(10, 3, sigt(10, 0, 0.05), sigr(10, 0.1)) # M>1, all under H1
# Example: Y = gen_shifted_data(10, 20, sigt(10, 0, 0.05), sigr(10, 0.1); nH1=5) # M>1, 5 under H1 and 15 under H0
function gen_shifted_data(N::Int, M::Int, μ::Float64, ρ; nH1=:all)
    Nd2=N÷2
    if M==1
        D=Normal()
        return vcat(rand(D, Nd2), rand(D, Nd2).+μ)
    else
        R=Matrix(1.0I, M, M) # desired correlation matrix for ALL variables
        for j=1:size(R, 1), i=j+1:size(R, 2) # correlation among all variables 2:M. it suffice to write the upper off diagonal elements
            R[j, i]=ρ
        end

        D=MvNormal(fill(0.0, M), Matrix(Hermitian(R)))
        μ1 = nH1===:all ? repeat([μ], M) : [repeat([μ], nH1); zeros(Float64, M-nH1)]
        X=rand(D, Nd2)
        Y=rand(D, Nd2)
        return [vcat(X[i, :], Y[i, :].+μ1[i]) for i=1:M]
    end
end

# Return the mean a sample of N÷2 elements must have so as to be significantly smaller at the α level 
# as compared to μ1, the mean of a sample of N÷2 elements to be compared to using the 
# unidirectional t-test for independent samples.
# The UNBIASED variance of the two samples is assumed equal to 1.
# According to a t-test with pooled unbiased standard deviation = 1, sigt(N, μ1, α)/sqrt(2/N) = t,
# where t is the value of the student t-statistic yielding a p-value equal to α. 
# NB: Generating random samples of the second group with this mean the expected power of a most uniformly powerful
# test is 0.5 using the Normal distribution for generating the data. This is so because because the expected proportion
# of samples with mean more extreme or equal to this mean is 0.5.  
sigt(N, μ1, α) = μ1 - sqrt(2/N) * cquantile(TDist(N-2), α)

# parametric right-directional t-test for uncorrelated samples
tTestISAsy(𝐱, 𝐲, N, ns) = ccdf(TDist(N-2), statistic(𝐱, 𝐲, StudentT_IS(); k=2, ns=ns))

##--


### UNIVARIATE TESTS
############################################################

# set the monte carlo simulations
Nrep=1000
#x=[i/Nrep for i=1:Nrep] # x-data for p-p plots
s2r=Int(1e5) # switch 2 random if more than 1e5 permutations

## Correlation test

N=[6, 9, 36, 100] # sample size
A=Matrix{Float64}(undef, Nrep, length(N)) # p-values of asymptotic test
P=Matrix{Float64}(undef, Nrep, length(N)) # p-values of permutation test
α=0.01 # expected correlation p-value
T1E=0.05 # type 1 error
for n=1:length(N), r=1:Nrep
    𝐱, 𝐲 = gen_correlated_data(N[n], 2, sigr(N[n], α), 0) # setting 0 instead of sigr(N[n], α) the hyp is under H0 
    P[r, n] = rTest(𝐱, 𝐲; switch2rand=s2r, direction=Right(), verbose=false).p
    A[r, n] = rTestAsy(𝐱, 𝐲, N[n]) 
end

# proportion of rejected hypotheses
propRejP=[count(≤(T1E), P[:, n])/Nrep for n=1:length(N)]
propRejA=[count(≤(T1E), A[:, n])/Nrep for n=1:length(N)]

nam = repeat("N=" .* string.(N), outer = 2)
gbkwargs=(tickfontsize=11, labelfontsize=14, legendfontsize=10, leg=:bottom,
            title="r-test, α=$α", label=["Param" "Perm"])
if Plots.backend()==Plots.PlotlyBackend()
    groupedbar(nam, hcat(propRejA, propRejP), ylabel="Power "*"(1-β)", leg=:bottom, tickfontsize=11, 
                title="Pearson correlation test, α=$(T1E), 𝔼p=$α", label=["Parametric" "Permutation"], legendfontsize=12, size=(800, 500))
#    hline!([0.5], label="𝔼(1-β)", lw=2, c=:grey) # the expected power when α = T1E is 0.5
else
    groupedbar(nam, hcat(propRejA, propRejP), label=["Param" "Perm"], ylabel="Power "*L"(1-β)", title="r-test, α=$(T1E), , 𝔼p=$α")
#    hline!([0.5], label=L"E(1-β)", lw=2, c=:grey) # the expected power when α = T1E is 0.5
end

hkwargs=(xlim=(0, 0.4), bins=100, alpha=0.7, color=[:cyan :yellow], label=["Param" "Perm"])
lkwargs=(lw=4, c=:grey, label="nominal "*L"α")
histogram(hcat(A[:, 1], P[:, 1]); title="Histogram of p-values (N=$(N[1]))", hkwargs...) # N=6
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 2], P[:, 2]); title="Histogram of p-values (N=$(N[2]))", hkwargs...) # N=9
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 3], P[:, 3]); title="Histogram of p-values (N=$(N[3]))", hkwargs...) # N=36
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 4], P[:, 4]); title="Histogram of p-values (N=$(N[3]))", hkwargs...) # N=100
vline!([T1E]; lkwargs...)

####################################################################
# Analysis of the behavior of the t-test correlation parametric test.
# This block of code is sot essential

N=[6, 9, 36, 100] # sample size
A=Matrix{Float64}(undef, Nrep, length(N)) # p-values of asymptotic test
P=Matrix{Float64}(undef, Nrep, length(N)) # p-values of permutation test
pTheo=[0.9998*(1-rand())+0.0001 for r=1:Nrep] # theoretical p-values ∈ (0.0001, 0.9999]
for n=1:length(N), r=1:Nrep
    𝐱, 𝐲 = gen_correlated_data(N[n], 2, sigr(N[n], pTheo[r]), 0) 
    P[r, n] = rTest(𝐱, 𝐲; switch2rand=s2r, direction=Right(), verbose=false).p
    A[r, n] = rTestAsy(𝐱, 𝐲, N[n]) 
end

# p observed - p theoretical
scatter(pTheo, A.-pTheo; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p observed - p theoretical", ms=2, ma=0.61, smooth=true, lw=3)
scatter(pTheo, P.-pTheo; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p observed - p theoretical", ms=2, ma=0.61, smooth=true, lw=3)

# the two below are equivalent: p obs Asy - p the - p obs Perm - p the  = p obs Asy - p obs Perm
D=A-P
scatter(pTheo, (A.-pTheo)-(P.-pTheo); label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p-value difference (Param-Perm)", ms=2, ma=0.61, lw=3)
scatter(pTheo, D; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p-value difference (Param-Perm)", ms=2, ma=0.61, lw=3) # smooth=true

# mean of the differences A-B for each value of N
meandiff=[mean(D[:, n]) for n=1:length(N)]
# one-sample t-test on D to test H0: A-B = 0 for each value of N
pdiff = [tTest1S(D[:, n]).p for n=1:length(N)]

# Another way to do it: difference of z-tranformed p-values obtained by the asymptotic and permutation test
Dz=[cquantile(Normal(), A[r, n])-cquantile(Normal(), P[r, n]) for r=1:Nrep, n=1:length(N)]
scatter(pTheo, Dz; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="z(p-value) difference (Param-Perm)", ms=2, ma=0.61, smooth=true, lw=3)

# mean of the differences z(A)-z(B) for each value of N
meandiff=[mean(Dz[:, n]) for n=1:length(N)]
# one-sample t-test on Dz to test H0: z(A)-z(B) = 0 for each value of N
pdiff = [tTest1S(Dz[:, n]).p for n=1:length(N)]

####################################################################

## t-test for indepedent samples

N=[12, 24, 48, 200] # total sample size, must be even
A=Matrix{Float64}(undef, Nrep, length(N)) # p-values of parametric test
P=Matrix{Float64}(undef, Nrep, length(N)) # p-values of permutation test
α=0.01 # expected t-test p-value
T1E=0.05 # type 1 error
for n=1:length(N) 
    ns=[N[n]÷2, N[n]÷2]
    𝐱=membership(StudentT_IS(), ns)
    for r=1:Nrep
        𝐲 = gen_shifted_data(N[n], 1, sigt(N[n], 0., α), 0) 
        P[r, n] = tTestIS(𝐲, [N[n]÷2, N[n]÷2]; switch2rand=s2r, direction=Right(), verbose=false, asPearson=false).p
        A[r, n] = tTestISAsy(𝐱, 𝐲, N[n], ns)
    end
end

# proportion of rejected hypotheses
propRejP=[count(≤(T1E), P[:, n])/Nrep for n=1:length(N)]
propRejA=[count(≤(T1E), A[:, n])/Nrep for n=1:length(N)]
nam = repeat("N=" .* string.(N), outer = 2)
gbkwargs=(tickfontsize=11, labelfontsize=14, legendfontsize=12, leg=:bottom,
            title="r-test, α=$(T1E), 𝔼p=$α", label=["Param" "Perm"])
if Plots.backend()==Plots.PlotlyBackend()
    groupedbar(nam, hcat(propRejA, propRejP), ylabel="Power "*"(1-β)", leg=:bottom, tickfontsize=11, 
                title="Student's t test Independent Samples, , α=$(T1E), 𝔼p=$α", label=["Parametric" "Permutation"], legendfontsize=12, size=(800, 500))
#    hline!([0.5], label="𝔼(1-β)", lw=2, c=:grey) # the expected power when α = T1E is 0.5
else
    groupedbar(nam, hcat(propRejA, propRejP), label=["Param" "Perm"], ylabel="Power "*L"(1-β)", title="r-test, , α=$(T1E), 𝔼p=$α")
#    hline!([0.5], label=L"E(1-β)", lw=2, c=:grey) # the expected power when α = T1E is 0.5
end


hkwargs=(xlim=(0, 0.4), bins=100, alpha=0.7, color=[:cyan :yellow], label=["Param" "Perm"])
lkwargs=(lw=4, c=:grey, label="nominal "*L"α")
histogram(hcat(A[:, 1], P[:, 1]); title="Histogram of p-values (N=$(N[1]))", hkwargs...) # N=12
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 2], P[:, 2]); title="Histogram of p-values (N=$(N[2]))", hkwargs...) # N=24
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 3], P[:, 3]); title="Histogram of p-values (N=$(N[3]))", hkwargs...) # N=48
vline!([T1E]; lkwargs...)
histogram(hcat(A[:, 4], P[:, 4]); title="Histogram of p-values (N=$(N[3]))", hkwargs...) # N=100
vline!([T1E]; lkwargs...)

###############################################################################################
# Analysis of the behavior of the t-test.
# This block of code is sot essential

N=[12, 24, 48, 200] # total sample size, must be even
A=Matrix{Float64}(undef, Nrep, length(N)) # p-values of parametric test
P=Matrix{Float64}(undef, Nrep, length(N)) # p-values of permutation test
pTheo=[0.9998*(1-rand())+0.0001 for r=1:Nrep] # theoretical p-values ∈ (0.0001, 0.9999]
for n=1:length(N) 
    ns=[N[n]÷2, N[n]÷2]
    𝐱=membership(StudentT_IS(), ns)
    for r=1:Nrep 
        𝐲 = gen_shifted_data(N[n], 1, sigt(N[n], 0., pTheo[r]), 0.) 
        P[r, n] = tTestIS(𝐲, [N[n]÷2, N[n]÷2]; switch2rand=s2r, direction=Right(), verbose=false).p
        A[r, n] = tTestISAsy(𝐱, 𝐲, N[n], ns) 
    end
end

# p observed - p theoretical
scatter(pTheo, A.-pTheo; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p obseved-theoretical", ms=2, ma=0.61, smooth=true, lw=3)
scatter(pTheo, P.-pTheo; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p observed-theoretical", ms=2, ma=0.61, smooth=true, lw=3)

# the two below are equivalent: ((p obs Para) - (p the)) - ((p obs Perm) - (p the))  = (p obs Para) - (p obs Perm)
D=A-P
scatter(pTheo, (A.-pTheo)-(P.-pTheo); label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p-value difference (Para-Perm)", ms=2, ma=0.61, lw=3)
scatter(pTheo, D; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="p-value difference (Para-Perm)", ms=2, ma=0.61, lw=3) #smooth=true

# mean of the differences A-B for each value of N
meandiff=[mean(D[:, n]) for n=1:length(N)]
# one-sample t-test on D to test H0: A-B = 0 for each value of N
pdiff = [tTest1S(D[:, n]).p for n=1:length(N)]

# Another way to do it: difference of z-tranformed p-values obtained by the parametric and permutation test
Dz=[cquantile(Normal(), A[r, n])-cquantile(Normal(), P[r, n]) for r=1:Nrep, n=1:length(N)]
scatter(pTheo, Dz; label="N=".*string.(N'), xlims=(0, 1), xlabel="theoretical p", ylabel="z(p-value) difference (Para-Perm)", ms=2, ma=0.61, smooth=true, lw=3)

# mean of the differences z(A)-z(B) for each value of N
meandiff=[mean(Dz[:, n]) for n=1:length(N)]
# one-sample t-test on Dz to test H0: z(A)-z(B) = 0 for each value of N
pdiff = [tTest1S(Dz[:, n]).p for n=1:length(N)]
###############################################################################################



############################################################
## MULTIPLE COMPARISONS TESTS
############################################################

# set up monte carlo simutaions

Nrep=100 # number of repetitions
s2r=Int(1e5) # switch 2 random if more than 1e5 permutations

## Correlation test
n=24 
α=0.01 # expected p-value of the correlation between x and all y variables
αρ=0.1 # # expected p-value of the correlation between all y variables
FWE=0.05 # Family-wise error
M=[5, 20, 100] # number of hypotheses
P=Matrix{Float64}(undef, Nrep, length(M))
A=Matrix{Float64}(undef, Nrep, length(M))
nrPerms(PearsonR(), n, allPerms(PearsonR(), n), Right())
for i=1:length(M)
    for r=1:Nrep
        𝐱, 𝐲 = gen_correlated_data(n, M[i], sigr(n, α), αρ≈1 ? 0.0 : sigr(n, αρ)) 
        P[r, i] = count(≤(FWE), rMcTest(𝐱, 𝐲; 
                switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, direction=Right(), stepdown=true).p)
        A[r, i] = count(≤(FWE/M[i]), [rTestAsy(𝐱, y, n) for y ∈ 𝐲]) # bonferroni corrected
    end
    println("done for m=", M[i])
end
# average number of rejections / number of hypotheses under H1 (power) 
P_=(sum(P, dims=1)./Nrep)./M'
A_=(sum(A, dims=1)./Nrep)./M'

# make a table
labels=["Permutation" "Parametric"]
println(" Average rejections / hypotheses under H1 (Power)")
print(" correlation test, ", n," observations")
pretty_table(hcat(Int.(M), P_', A_'); formatters = ft_printf("%5.2f", 2:3), header = ( hcat("M", labels),), header_crayon = crayon"yellow bold")

# difference between the rejection of the parametric - permutation test
D=P-A
scatter(Base.OneTo(Nrep), D, label="M=".*string.(M'), xlabel="Simulations", ylabel="rejections Perm - rejections Param", 
                ms=3, ma=0.61, tickfontsize=11, labelfontsize=13, legendfontsize=11)


## t-test for independent samples

M=[5, 20, 100] # number of hypotheses
ns=[8, 8] # must be balanced, see gen_shifted_data
n=sum(ns) # n1+n2
𝐱=membership(StudentT_IS(), ns)
α=0.01 # expected p-value of the correlation between x and all y variables
αρ=0.1 # # expected p-value of the correlation between all y variables
nH1=[M[i]÷2 for i=1:length(M)]
FWE=0.05 # Family-wise error

P=Matrix{Float64}(undef, Nrep, length(M))
A=Matrix{Float64}(undef, Nrep, length(M))
nrPerms(StudentT_IS(), ns, allPerms(StudentT_IS(), ns), Right(), Balanced())
for i=1:length(M)
    for r=1:Nrep
        𝐲 = gen_shifted_data(n, M[i], sigt(n, 0., α), αρ≈1 ? 0.0 : sigr(n, αρ); nH1=nH1[i]) 
        P[r, i] = count(≤(FWE), tMcTestIS(𝐲, ns; 
                switch2rand=s2r, verbose=false, threaded=true, fwe=FWE, direction=Right(), stepdown=true, asPearson=false).p)
        A[r, i] = count(≤(FWE/M[i]), [tTestISAsy(𝐱, y, n, ns) for y ∈ 𝐲]) # bonferroni corrected 
    end
    println("done for m=", M[i])
end
# average rejections / number of hypotheses under H1 (power) - make a table
P_=(sum(P, dims=1)./Nrep)./nH1' 
A_=(sum(A, dims=1)./Nrep)./nH1'

# make a table
labels=["Permutation" "Parametric"]
println(" Average rejections / hypotheses under H1 (Power)")
print(" Student's t test for Ind. Samples, ", ns," subjects per group")
pretty_table(hcat(Int.(M), P_', A_'); formatters = ft_printf("%5.2f", 2:3), header = ( hcat("M", labels),), header_crayon = crayon"yellow bold")


# difference between the rejection of the parametric - permutation test
D=P-A
scatter(Base.OneTo(Nrep), D, label="M=".*string.(M'), xlabel="Simulations", ylabel="rejections Perm - rejections Param", 
                ms=3, ma=0.61, tickfontsize=11, labelfontsize=13, legendfontsize=11)

###