#   Unit "benchmarks.jl", additional resource to the PemutationTests.jl package for the julia language
#
#   MIT License
#   Copyright (c) 2024,
#   Marco Congedo, CNRS, Grenoble, France:
#   https://sites.google.com/site/marcocongedo/home

#   CONTENTS :
#   This unit benchmark the speed of execution of the univarate and multiple comparison tests
#   implemented in PermutationTests.jl. It creates tables and plots of the minimum execution time in milliseconds,
#   memory usage in Megabytes and memory allocation count.
#   NB: for univariate tests we use @threads for speeding up the monte carlo simulations, whereas for
#   multiple comparisons tests we don't as the main testing function is multi-threaded, thus we leave the logic 
#   processing unit free to be used by the main testing function.    

#   REQUIRED PACKAGES : PermutationsTests, BenchmarkTools, PrettyTables, Distributions, Plots
#   NB: if with Plots.jl you use the plotly() back-end, as in the code, 
#       you may need to install packages 'PlotlyBase' and 'PlotlyKaleido'

using Random:MersenneTwister # julia built-in
using PrettyTables, Plots, Plots.Measures
using BenchmarkTools, Base.Threads
using Distributions:Binomial
using PermutationTests
plotly() # select Plots.jl backend

# time will be given in milliseconds
function addResults!(np, m, T, G, M, A, i, j; timeInSec=false)
    T[i, j] = timeInSec ? m.time/(1e9) : m.time/(1e6) # milliseconds, seconds if timeInSec=true 
    G[i, j] = m.gctime/1e6 #milliseconds
    M[i, j] = round(Int, m.memory/1e3) # Megabyte
    A[i, j] = m.allocs
end


rng = MersenneTwister(1234)

# Create a table
function gettable(n, k, rng)
    table=zeros(Int64, 2, k)
    nk=n÷k
    for i=1:size(table, 2)
        table[1, i]=rand(rng, Binomial((n/2)÷k, 0.5))
        table[2, i]=table[1, i]<nk ? (nk-table[1, i]) : 0
    end
    # correct the table is the sum of the columns sum > N
    for i=1:size(table, 2)
        cs=sum(table[:, i])
        while cs>nk
            if table[1, i]>0
                table[1, i]-=1
            end
            cs=sum(table[:, i])
            if cs>nk
                if table[2, i]>0
                    table[2, i]-=1
                end
            end
        end
    end
    # correct the table if k is even and n is odd
    s=sum(table)
    while s<n
        table[1, 1]+=1
        s=sum(table)
    end
    return table
end


## UNIVARIATE TESTS
####################################################################
N=[9, 90, 900, 9000, 90000]
# All tests are run with `numPerm` random permutations
numPerm=10000
nTests=11
nn=length(N)

T=Matrix{Float64}(undef, nn, nTests)    # time
G=Matrix{Float64}(undef, nn, nTests)    # garbage collector time
M=Matrix{Int64}(undef, nn, nTests)      # Memory
A=Matrix{Int64}(undef, nn, nTests)      # Allocations

# Correlation test
@threads for i = 1:length(N)
    n = N[i]
    x = randn(n)
    y = randn(n)
    t = @benchmark rTest($x, $y; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 1)
    println("Correlation Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# 1-way ANOVA for indepedent samples
@threads for i = 1:length(N)
    n = N[i]
    x = randn(n) 
    ns = [n÷3, n÷3, n-((n÷3)*2)]
    t = @benchmark fTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 2)
    println("1-way ANOVA IS Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# t-test for indepedent samples
@threads for i = 1:length(N)
    n = N[i]
    x = randn(n) 
    ns = [n÷2, n-(n÷2)]
    t = @benchmark tTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false, asPearson=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 3)
    println("Student t IS Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# t-test for indepedent samples by Pearson Correlation
@threads for i = 1:length(N)
    n = N[i]
    x = randn(n) 
    ns = [n÷2, n-(n÷2)]
    t = @benchmark tTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false, asPearson=true)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 4)
    println("Student t IS Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# Chi-Square
@threads for i = 1:length(N)
    n = N[i]
    k = 3
    x = gettable(n, k, rng)
    t = @benchmark Χ²Test($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 5)
    println("Chi-Square Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# Fisher Exact Test
@threads for i = 1:length(N)
    n = N[i]
    k = 2
    x = gettable(n, k, rng)
    t = @benchmark fisherExactTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 6)
    println("Fisher exact Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# 1-way ANOVA for repeated samples
@threads for i = 1:length(N)
    n = N[i]
    k = 3
    x = randn(n*k)
    ns = (n=n, k=k)
    t = @benchmark fTestRM($x, $ns; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 7)
    println("1-way ANOVA RM Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# Cochran Q test
@threads for i = 1:length(N)
    n = N[i]
    k = 3
    x = Int.(rand(Bool, n, k))
    t = @benchmark qTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 8)
    println("Cochran Q Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# McNemar test
@threads for i = 1:length(N)
    n = N[i]
    k = 2
    x = Int.(rand(Bool, n, k))
    t = @benchmark mcNemarTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 9)
    println("McNemar Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# one sample t-test
@threads for i = 1:length(N)
    n = N[i]
    x = randn(n)
    t = @benchmark tTest1S($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 10)
    println("One-Sample Student-t Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# Sign Test
@threads for i = 1:length(N)
    n = N[i]
    x = rand(Bool, n)
    t = @benchmark signTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 11)
    println("Sign Test. Done for N = ", N[i], ", Minimum Time(ms): ", minimum(t).time/1e6)
end

# Tables
labels=["Pearson r" "ANOVA F IS" "St. t IS" "St. t IS(r)" "χ²" "Fisher Exact" "ANOVA F RM" "Cochran Q" "McNemar" "St. t 1S" "Sign"]
print(" Minumum Execution Time (ms), ", numPerm," Random Permutations")
pretty_table(hcat(Int.(N), T); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("N", labels),), header_crayon = crayon"yellow bold")
print(" Memory Allocations (MB), ", numPerm," Random Permutations")
pretty_table(hcat(Int.(N), M); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("N", labels),), header_crayon = crayon"yellow bold")
print(" Memory Allocation count, ", numPerm," Random Permutations")
pretty_table(hcat(Int.(N), A); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("N", labels),), header_crayon = crayon"yellow bold")

# Plots
labels=["Pearson r" "ANOVA F IS" "Student t IS" "Student t IS as r" "χ²" "Fisher Exact" "ANOVA F RM" "Cochran Q" "McNemar" "One-sample Student t" "Sign"]
kwargs=(yscale=:log10, minorgrid=true, label=labels, xlabel="Number of observations", leg=:topleft, 
        left_margin = 25mm, bottom_margin = 5mm, right_margin = 5mm, lw=2, xticks = ([1:1length(N);], string.(N)),
        size=(900, 600), tickfontsize=11, labelfontsize=14, legendfontsize=10,)

# temp: remove once the tests are done also for the Chi-Square test
plot(T; ylabel="Minimum execution time (ms)", title="$numPerm permutations", ylims=(1, 1e5), kwargs...)

plot(M; ylabel="Memory allocation (MB)", title="$numPerm permutations", ylims=(1e2, 1e5), kwargs...)

plot(M; ylabel="allocation count", title="$numPerm permutations", ylims=(1e5, 1e7*3), kwargs...)


## Multiple Comparison TESTS
####################################################################
M_=[10, 100, 1000, 10000]
N=12
numPerm=10000
nTests=11
nm=length(M_)

T=Matrix{Float64}(undef, nm, nTests)    # time
G=Matrix{Float64}(undef, nm, nTests)    # garbage collector time
M=Matrix{Int64}(undef, nm, nTests)      # Memory
A=Matrix{Int64}(undef, nm, nTests)      # Allocations


# Correlation test
x = randn(N)
for i=1:length(M_)
    y = [randn(N) for j=1:M_[i]]
    t = @benchmark rMcTest($x, $y; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 1; timeInSec=true)
    println("Correlation Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# 1-way ANOVA for indepedent samples 
for i=1:length(M_)
    x = [randn(N) for j=1:M_[i]]
    ns=[N÷3, N÷3, N-((N÷3)*2)]
    t = @benchmark fMcTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 2; timeInSec=true)
    println("1-way ANOVA IS Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# t-test for indepedent samples
for i=1:length(M_)
    x = [randn(N) for j=1:M_[i]]
    ns=[N÷2, N-(N÷2)]
    t = @benchmark tMcTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false, asPearson=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 3; timeInSec=true)
    println("Student t IS Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# t-test for indepedent samples as PearsonR
for i=1:length(M_)
    x = [randn(N) for j=1:M_[i]]
    ns=[N÷2, N-(N÷2)]
    t = @benchmark tMcTestIS($x, $ns; switch2rand=1, nperm=numPerm, verbose=false, asPearson=true)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 4; timeInSec=true)
    println("Student t IS Test as Pearson. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# Chi-Square Test
for i=1:length(M_)
    k = 3
    x = [gettable(N, k, rng) for j=1:M_[i]]
    t = @benchmark Χ²McTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 5; timeInSec=true)
    println("Χ² Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# Fisher Exact Test
for i=1:length(M_)
    k = 2
    x = [gettable(N, k, rng) for j=1:M_[i]]
    t = @benchmark fisherExactMcTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 6; timeInSec=true)
    println("Fisher Exact Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# 1-way ANOVA for repeated measures
for i=1:length(M_)
    x = [randn(N) for j=1:M_[i]]
    ns=(n=N÷3, k=3)
    t = @benchmark fMcTestRM($x, $ns; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 7; timeInSec=true)
    println("1-way ANOVA RM Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# Cochran Q Test
for i=1:length(M_)
    k = 3
    x = [Int.(rand(Bool, N÷k, k)) for j=1:M_[i]]
    t = @benchmark qMcTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 8; timeInSec=true)
    println("Cochran Q Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# McNemar Test
for i=1:length(M_)
    k = 3
    x = [Int.(rand(Bool, N÷k, k)) for j=1:M_[i]]
    t = @benchmark qMcTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 9; timeInSec=true)
    println("McNemar Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# One-sample t-test
for i=1:length(M_)
    x = [randn(N) for j=1:M_[i]]
    t = @benchmark tMcTest1S($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 10; timeInSec=true)
    println("One-Sample Student-t Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# Sign Test
for i=1:length(M_)
    x = [rand(Bool, N)  for j=1:M_[i]]
    t = @benchmark signMcTest($x; switch2rand=1, nperm=numPerm, verbose=false)
    addResults!(numPerm, minimum(t), T, G, M, A, i, 11; timeInSec=true)
    println("Sign Test. Done for M = ", M_[i], ", Minimum Time(s): ", minimum(t).time/1e9)
end

# Tables
labels=["Pearson r" "ANOVA F IS" "St. t IS" "St. t IS(r)" "χ²" "Fisher Exact" "ANOVA F RM" "Cochran Q" "McNemar" "St. t 1S" "Sign"]
print(" Minumum Execution Time (s), ", numPerm," Random Permutations, N=", N)
pretty_table(hcat(Int.(M_), T); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("M", labels),), header_crayon = crayon"yellow bold")
print(" Memory Allocations (MB), ", numPerm," Random Permutations, N=", N)
pretty_table(hcat(Int.(M_), M); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("M", labels),), header_crayon = crayon"yellow bold")
print(" Memory Allocation count, ", numPerm," Random Permutations, N=", N)
pretty_table(hcat(Int.(M_), A); formatters = ft_printf("%5.2f", 2:4), header = ( hcat("M", labels),), header_crayon = crayon"yellow bold")

# Plots
labels=["Pearson r" "ANOVA F IS" "Student t IS" "Student t IS as r" "χ²" "Fisher Exact" "ANOVA F RM" "Cochran Q" "McNemar" "One-sample Student t" "Sign"]
kwargs=(yscale=:log10, minorgrid=true, label=labels, xlabel="Number of hypotheses", leg=:topleft, 
        left_margin = 35mm, bottom_margin = 5mm, right_margin = 5mm, lw=2, xticks = ([1:1length(M_);], string.(M_)),
        size=(900, 600), tickfontsize=11, labelfontsize=14, legendfontsize=10,)

# temp: remove once the tests are done also for the Chi-Square test
plot(T; ylabel="Minimum Execution time (s)", title="$numPerm permutations", ylims=(0.1, 1e2), kwargs...)

plot(M; ylabel="Memory allocation (MB)", title="$numPerm permutations", ylims=(1e4*3, 1e8*2), kwargs...)

plot(M; ylabel="allocation count", title="$numPerm permutations", ylims=(1e4*3, 1e9), kwargs...)

