#   Unit "runtests.jl" of the PemutationTests.jl package
#
#   MIT License
#   Copyright (c) 2024,
#   Marco Congedo, CNRS, Grenoble, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENT:
#   test_stats()
#   test_unitests()
#   test_unitests_API()
#   test_multicompTests()
#   test_multicompTests_API()
#   The test of all the above is at the end of the script 

#   REQUIRED PACKAGES : PermutationsTests, Test, Statistics

using Random: randn # julia built-in
using Test: @test, @testset
using Statistics: mean, var, std, cov, cor
using PermutationTests


# test all functions in stats.jl
function test_stats(; tol=0)

    println("============ Basic Statistics functions =================")
    println()
    
    n=10
    𝐱=abs.(randn(n))
    𝐲=abs.(randn(n))
    m=mean(𝐱)
    v=var(𝐱; corrected=false)
    sd=std(𝐱; corrected=false)
    𝐱c=𝐱.-m

    println("∑(𝐱), ∑of²(𝐱), ∑∑of²(𝐱), Π(𝐱), ∑ofΠ(𝐱, 𝐲), μ(𝐱)")
    @testset "\nSums and Products                 " begin
        @test ∑(𝐱) ≈ sum(𝐱);
        @test  ∑of²(𝐱) ≈ sum(𝐱.^2);
        @test  sum(∑∑of²(𝐱)) ≈ sum(𝐱) + sum(𝐱.^2);
        @test Π(𝐱) ≈ prod(𝐱);
        @test ∑ofΠ(𝐱, 𝐲) ≈ sum(𝐱.*𝐲);
        @test μ(𝐱) ≈ mean(𝐱);
    end

    println("dispersion(𝐱), σ²(𝐱), σ(𝐱)")
    @testset "\nMeasures of dispersion            " begin
        @test isapprox(dispersion(𝐱; mean=m), v*n, atol=tol);
        @test isapprox(dispersion(𝐱c; centered=true), var(𝐱c; corrected=false)*n, atol=tol);
        @test isapprox(dispersion(𝐱), v*n, atol=tol); 
        @test isapprox(σ(𝐱), sd; atol=tol);
        @test isapprox(σ(𝐱; mean=m), sd, atol=tol);
        @test isapprox(σ(𝐱c; centered=true), std(𝐱c; corrected=false), atol=tol); 
    end
    
    println("μ0(𝐱), σ1(𝐱), μ0σ1(𝐱)")
    @testset "\nData Transformations              " begin
        @test isapprox(mean(μ0(𝐱)), 0, atol=1e-1); 
        @test isapprox(mean(μ0(𝐱; mean=m)), 0, atol=1e-1); 
        @test isapprox(σ(σ1(𝐱)), 1, atol=tol); 
        @test isapprox(σ(σ1(𝐱; mean=m, sd=sd)), 1, atol=tol);
        a=μ0(𝐱)
        @test isapprox(σ(σ1(a; centered=true, sd=sd)), 1, atol=tol); 
        @test isapprox(σ(σ1(a; centered=true)), 1, atol=tol); 
        @test isapprox(σ(σ1(𝐱; mean=m)), 1, atol=tol); 
        @test isapprox(σ(σ1(𝐱; sd=sd)), 1, atol=tol); 

        a=μ0σ1(𝐱)
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1));
        a=μ0σ1(𝐱; sd=sd)
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1));
        a=μ0σ1(𝐱; mean=m)
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1));
        a=μ0σ1(𝐱; mean=m, sd=sd)
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1));
        b=μ0(𝐱)
        a=μ0σ1(b; centered=true)
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1));
        a=μ0σ1(b; centered=true, sd=σ(b; centered=true))
        @test (isapprox(μ(a), 0, atol=1e-1) && isapprox(σ(a), 1, atol=1e-1)); 
    end
    
    
    c=cov(𝐱, 𝐲; corrected=false)
    𝐱0μ=μ0(𝐱)
    𝐲0μ=μ0(𝐲)
    c0=cov(𝐱0μ, 𝐲0μ; corrected=false)

    println("CrossProd(), Covariance(), PearsonR()")
    @testset "\nBivariate Statistics              " begin
        @test statistic(𝐱, 𝐲, CrossProd()) ≈ ∑ofΠ(𝐱, 𝐲);   
        @test statistic(𝐱0μ, 𝐲0μ, Covariance(); centered=true, means=()) ≈ c0 ;   
        @test statistic(𝐱, 𝐲, Covariance(); means=(mean(𝐱), mean(𝐲))) ≈ c;  
        @test statistic(𝐱, 𝐲, Covariance()) ≈ c;  

        r=cor(𝐱, 𝐲)
        @test statistic(𝐱, 𝐲, PearsonR()) ≈ r ;  
        @test statistic(𝐱, 𝐲, PearsonR(); means=(μ(𝐱), μ(𝐲))) ≈ r ;       
        @test statistic(𝐱, 𝐲, PearsonR(); means=(μ(𝐱), μ(𝐲)), sds=(σ(𝐱), σ(𝐲))) ≈ r ;  
        𝐱_=μ0σ1(𝐱)
        𝐲_=μ0σ1(𝐲)
        @test statistic(𝐱_, 𝐲_, PearsonR(); standardized=true) ≈ r ;  
    end

    println("AnovaF_IS(), SumGroupTotalsSq_IS(), Group1Total_IS(), StudentT_IS()")
    @testset "\nIndependent Sample Statistics     " begin
        𝐱 = [1, 1, 1, 2, 2, 2]
        𝐲 = [1.3, 2.5, 0.5, 2.5, 6.3, 1.4]
        k=length(unique(𝐱)) 
        ns=[count(x->x==i, 𝐱) for i=1:k]
        @test isapprox(statistic(𝐱, 𝐲, AnovaF_IS(); k=k, ns=ns), 1.52214, atol=1e-1); 
        @test isapprox(statistic(𝐱, 𝐲, AnovaF_IS()), 1.52214, atol=1e-1);
        𝐱 = [1, 1, 2, 2]
        𝐲 = [1, 2, 3, 4]
        @test isapprox(statistic(𝐱, 𝐲, SumGroupTotalsSq_IS(); k=2), (1+2)^2+(3+4)^2, atol=tol); 
        @test isapprox(statistic(𝐱, 𝐲, SumGroupTotalsSqN_IS(); k=2, ns=[2, 2]), ((1+2)^2)/2+((3+4)^2)/2, atol=tol); 

        𝐱=[1, 1, 1, 2, 2, 2]
        𝐲=[1., 2., 3., 1.5, 2.5, 2.9]
        @test isapprox(statistic(𝐱, 𝐲, Group1Total_IS()), 6.0, atol=tol); 
        ns=[3, 3]
        @test isapprox(statistic(𝐱, 𝐲, StudentT_IS(); k=2, ns=ns), -0.42146361521, atol=tol); 
    end

    println("\nAnovaF_RM(), SumTreatTotalsSq_RM()")
    @testset "Repeated-Measures Statistics     " begin
        # (Edgington, page 103-104)
        n, k = 3, 3
        𝐱=collect(1:n*k)
        𝐲 = [10., 11., 13., 13., 14., 17., 6., 8., 9.]
        ns=(n=n, k=k)
        @test isapprox(statistic(𝐱, 𝐲, SumTreatTotalsSq_RM(); ns), 3451, atol=tol); 


        # https://statistics.laerd.com/statistical-guides/repeated-measures-anova-statistical-guide-2.php
        n, k = 6, 3
        𝐱=collect(1:n*k)
        𝐲=[45., 50., 55., 42., 42., 45., 36., 41., 43., 39., 35., 40., 51., 55., 59., 44., 49., 56.]
        ns=(n=n, k=k)
        @test isapprox(statistic(𝐱, 𝐲, AnovaF_RM(); ns=ns), 12.53, atol=1e-1); 
    end

    println("\nStudentT_1S(), Sum()")
    @testset "One-Sample Statistics             " begin
        # https://www.machinelearningplus.com/statistics/one-sample-t-test/
        𝐲=[21.5, 24.5, 18.5, 17.2, 14.5, 23.2, 22.1, 20.5, 19.4, 18.1, 24.1, 18.5].-20
        𝐱=[1.]
        @test isapprox(statistic(𝐱, 𝐲, StudentT_1S()), 0.2006, atol=1e-1);
        sumy²=sum(abs2, 𝐲)
        @test isapprox(statistic(𝐱, 𝐲, StudentT_1S(), ∑y²=sumy²), 0.2006, atol=1e-1); 
        @test isapprox(statistic(𝐱, 𝐲, Sum(); testtype = :approximate), sum(𝐲), atol=1e-1); 
    end
    println()
    return true
end


#####################################################################################
# Test all functions in uniTests.jl

function test_unitests(; tol=1e-1)
     
    println("============ Univariate test functions =================")
    println()
    
    # Product-moment Correlation (E.S. Edington, page 208)
    println("Product-moment Correlation tests : "); 
    𝐱=[2.0, 4.0, 6.0, 8.0, 10.0]
    𝐲=[1.0, 3.0, 5.0, 8.0, 7.0]
    #_permTest!(𝐱, 𝐲, length(𝐱), CrossProd(); fstat=identity)
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), CrossProd(), CrossProd(); fstat=identity, verbose=false).p ≈ (2/120) ; print("right-directional ✔, ") 
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), CrossProd(), CrossProd(); fstat=flip, verbose=false).p ≈ (119/120) ;  print("left-directional ✔, ")
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), Covariance(), Covariance(), verbose=false).p ≈ (4/120) ;  println("bi-directional ✔ ")
    println("")

    println("Correlation - trend tests : "); 
    𝐱=[2.0, 4.0, 6.0, 8.0, 10.0]
    𝐲=[1.0, 3.0, 5.0, 8.0, 7.0]
    n=length(𝐲)
    𝐱=collect(Base.OneTo(n))
    #_permTest!(𝐱, 𝐲, length(𝐱), CrossProd(); fstat=identity)
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), CrossProd(), CrossProd(); fstat=identity, verbose=false).p ≈ (2/120) ; print("right-directional ✔, ") 
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), CrossProd(), CrossProd(); fstat=flip, verbose=false).p ≈ (119/120) ;  print("left-directional ✔, ")
    @test _permTest!(copy(𝐱), 𝐲, length(𝐱), Covariance(), Covariance(), verbose=false).p ≈ (4/120) ;  println("bi-directional ✔ ")
    _permTest!(copy(𝐱), 𝐲, length(𝐱), CrossProd(), CrossProd(); fstat=identity, switch2rand=1, nperm=5000, seed=0, verbose=false)  # run approximate test
    println("")


    # 1-Way ANOVA for independent samples (E.S. Edington, page 62)
    𝐲=[1, 2, 3, 4, 5, 7, 8, 9, 10] # NB: input data can be also integer
    ns=[2, 3, 4]
    𝐱=membership(AnovaF_IS(), ns)
    println("1-way independent samples ANOVA test : "); 
    @test _permTest!(copy(𝐱), 𝐲, ns, AnovaF_IS(), AnovaF_IS(), verbose=false).p ≈ (2/1260) ;  print("F statistic ✔, ")
    print("Equivalence of ES SumGroupTotalsSqN_IS() : "); 
    @test _permTest!(copy(𝐱), 𝐲, ns, AnovaF_IS(), AnovaF_IS(), verbose=false).p ≈ _permTest!(copy(𝐱), 𝐲, ns, SumGroupTotalsSqN_IS(), SumGroupTotalsSqN_IS(), verbose=false).p ;  println(" ✔ ")


    println("Second test:")
    # 1-Way ANOVA for independent samples (E.S. Edgington, page 70)
    𝐲=[1, 2, 3, 4, 5, 6, 7, 8, 9] # NB: input data can be also integer
    ns=[2, 3, 4]
    𝐱=membership(AnovaF_IS(), ns)
    @test _permTest!(copy(𝐱), 𝐲, ns, AnovaF_IS(), AnovaF_IS(), verbose=false).p ≈ (6/1260) ;  print("F statistic ✔, ")
    @test _permTest!(copy(𝐱), 𝐲, ns, SumGroupTotalsSqN_IS(), SumGroupTotalsSqN_IS(), verbose=false).p ≈ (6/1260) ;  println("SumGroupTotalsSqN_IS statistic ✔ ")
    println("")

    println("test the number of non-redundant permutations ")
    ns=[3, 3, 3]
    stat=AnovaF_IS()
    @test allPerms(stat, ns) / nrPerms(stat, ns, allPerms(stat, ns), Both(), Balanced()) ≈ factorial(length(ns));  print("✔, ")
    _permTest!(copy(𝐱), 𝐲, ns, SumGroupTotalsSqN_IS(), SumGroupTotalsSqN_IS(), verbose=false, switch2rand=1, nperm=5000) # run approximate test 

    # T-tests (independent samples)
    𝐲=[4., 5., 6., 7., 1., 2., 3]
    ns=[4, 3]
    𝐱=membership(StudentT_IS(), ns)
    println("t-test independent samples (bi-directional): ");
    _permTest!(copy(𝐱), 𝐲, ns, StudentT_IS(), StudentT_IS(); switch2rand=1, seed=0, verbose=false)
    @test _permTest!(copy(𝐱), 𝐲, ns, StudentT_IS(), StudentT_IS(), verbose=false).p ≈ (2/35) ;  print("StudentT_IS statistic ✔, ")
    @test _permTest!(copy(𝐱), 𝐲, ns, SumGroupTotalsSqN_IS(), SumGroupTotalsSqN_IS(), verbose=false).p ≈ (2/35) ;  println("SumGroupTotalsSqN_IS statistic ✔, ")
    @test _permTest!(copy(𝐱), 𝐲, ns, StudentT_IS(), StudentT_IS(), fstat=identity, verbose=false).p ≈ (1/35) ;  print("StudentT_IS statistic rigth-directional ✔ ")
    @test _permTest!(copy(𝐱), 𝐲, ns, Group1Total_IS(), Group1Total_IS(), fstat=identity, verbose=false).p ≈ (1/35) ;  println("SumGroupTotalsSqN_IS statistic rigth-directional ✔ ")
    _permTest!(copy(𝐱), 𝐲, ns, Group1Total_IS(), Group1Total_IS(), fstat=identity, verbose=false, switch2rand=1, nperm=5000) # run approximate test
    println("")

    # 1-Way repeated-measures ANOVA (E.S. Edgington, page 103-104)
    n, k = 3, 3
    𝐲 = [10., 11., 13., 13., 14., 17., 6., 8., 9.]
    ns=(n=n, k=k)
    𝐱=membership(AnovaF_RM(), ns)
    println("1-Way Repeated Measure ANOVA :"); 
    @test _permTest!(copy(𝐱), 𝐲, ns, AnovaF_RM(), AnovaF_RM(); verbose=false).p ≈ (1/36); print("AnovaF_RM statistic ✔, ")

    print("Equivalence of SumTreatTotalsSq_RM statistic: "); 
    @test _permTest!(copy(𝐱), 𝐲, ns, SumTreatTotalsSq_RM(), SumTreatTotalsSq_RM(), verbose=false).p ≈ _permTest!(copy(𝐱), 𝐲, ns, SumTreatTotalsSq_RM(), SumTreatTotalsSq_RM(), verbose=false).p ;  println(" ✔ ")
    _permTest!(copy(𝐱), 𝐲, ns, AnovaF_RM(), AnovaF_RM(); switch2rand=1, nperm=5000, seed=0, verbose=false) # run approximate test

    println("test the number of non-redundant permutations ")
    ns=(n=6, k=3)
    stat=AnovaF_RM()
    @test allPerms(stat, ns) / nrPerms(stat, ns, allPerms(stat, ns), Both(), Balanced()) ≈ factorial(ns.k);  print("✔, ")   
 
    # One-sample t-test
    ns = 6
    𝐲 = randn(ns).+0.7
    𝐱=membership(StudentT_1S(), ns)
    println("One-sample t-test :");
    @test _permTest!(copy(𝐱), copy(𝐲), ns, StudentT_1S(), StudentT_1S(), verbose=false).p ≈ _permTest!(copy(𝐱), copy(𝐲), ns, Sum(), Sum(), verbose=false).p ;  println(" Equivalence of Sum ✔ ")

    # just run the approximate test
    ns = 18
    𝐱=Int[]
    𝐲 = randn(ns)#.+0.7
    _permTest!(copy(𝐱), copy(𝐲), ns, StudentT_1S(), StudentT_1S(); switch2rand=1, nperm=5000, seed=0, verbose=false) # run approximate test
    
    println()
    return true
end


#####################################################################################################""
# Test all functions in uniTests_API.jl

function test_uniTests_API(; tol=1e-1)

    println("============ API for univariate test functions =================")
    println()

    println("Correlation tests : "); 
    𝐱=[2.0, 4.0, 6.0, 8.0, 10.0]
    𝐲=[1.0, 3.0, 5.0, 8.0, 7.0]
    @test rTest(𝐱, 𝐲; direction = Right(), verbose=false).p ≈ (2/120) ; print("right-dir (r) ✔, ") 
    @test rTest(𝐱, 𝐲; direction = Left(), verbose=false).p ≈ (119/120) ; print("left-dir (r) ✔, ") 
    @test rTest(𝐱, 𝐲; direction = Both(), verbose=false).p ≈ (4/120) ; println("bi-dir (r) ✔, ")
    correlationTest(𝐱, 𝐲; direction = Both(), switch2rand=1, verbose=false) 
    println("")

    println("1-way ANOVA for independent samples (method 1): "); 
    𝐲=[1., 2., 3., 4., 5., 7., 8., 9., 10.] 
    ns=[2, 3, 4]
    𝐱=membership(AnovaF_IS(), ns)
    @test fTestIS(𝐲, ns; direction = Both(), verbose=false).p ≈ (2/1260) ; print("bi-dir (F) ✔, ") 
    𝐲=[[1., 2.], [3., 4., 5.],[7., 8., 9., 10.]]
    print(" method 2 ")
    @test fTestIS(𝐲, direction = Both(), verbose=false).p ≈ (2/1260) ; print("bi-dir (F) ✔, ") 
    𝐲=[[1., 2., 3], [3.1, 4., 5.],[7., 8., 9.]]
    ns=[3, 3, 3]
    np=nrPerms(AnovaF_IS(), ns, allPerms(AnovaF_IS(), ns), Both(), Balanced()) 
    @test fTestIS(𝐲; direction = Both(), verbose=false).nperm ≈ np  ; println("non-redundant permutations ✔, ") 
    println("")

    println("t-test independent samples (method 1): ");
    𝐲=[4., 5., 6., 7., 1., 2., 3]
    ns=[4, 3]
    @test tTestIS(𝐲, ns; direction = Both(), verbose=false).p ≈ (2/35) ; print("bi-dir (T) ✔, ") 
    @test tTestIS(𝐲, ns; direction = Right(), verbose=false).p ≈ (1/35) ; println("right-dir (T) ✔, ") 
    print("method 2 ")
    𝐲=[[4., 5., 6., 7.] , [1., 2., 3]]
    @test tTestIS(𝐲; direction = Both(), verbose=false).p ≈ (2/35) ; print("bi-dir (T) ✔, ") 
    @test tTestIS(𝐲; direction = Right(), verbose=false).p ≈ (1/35) ; println("right-dir (T) ✔, ") 
    print("method 3 ")
    𝐱=[4., 5., 6., 7.]
    𝐲=[1., 2., 3]
    @test tTestIS(𝐱, 𝐲; direction = Both(), verbose=false).p ≈ (2/35) ; print("bi-dir (T) ✔, ") 
    @test tTestIS(𝐱, 𝐲; direction = Right(), verbose=false).p ≈ (1/35) ; print("right-dir (T) ✔, ") 
    #
    𝐲=[4., 5., 6., 7., 0., 1., 2., 3]
    ns=[4, 4]
    np=nrPerms(StudentT_IS(), ns, allPerms(StudentT_IS(), ns), Both(), Balanced()) 
    @test tTestIS(𝐲, ns; direction = Both(), verbose=false, asPearson=false).nperm ≈ np  ; println("non-redundant permutationns ✔, ") 
    println("")

    println("Chi-Square and Fisher Exact Test : "); 
    table=[0 2 2; 3 1 0]
    @test Χ²Test(table; direction = Both(), verbose=false).p ≈ 0.22857142857142856  ; print("bi-dir (Χ²) ✔, ") 

    table=[1 7; 6 2] # asymptotic Χ²-test: Χ²=6.3492, p=0.011743 https://stats.libretexts.org/Learning_Objects/02%3A_Interactive_Statistics/37%3A__Chi-Square_Test_For_Independence_Calculator
    @test Χ²Test(table; direction = Both(), verbose=false, asPearson=false).p ≈ 0.04055944055944056  ; print("bi-dir (Χ²) ✔, ") 
    @test Χ²Test(table; direction = Right(), verbose=false, asPearson=false).p ≈ 0.9993006993006993  ; print("right-dir (Χ²) ✔, ") 
    @test Χ²Test(table; direction = Left(), verbose=false, asPearson=false).p ≈ 0.02027972027972028 ; println("left-dir (Χ²) ✔, ") 
    
    table=[11 3; 1 9] # https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    @test Χ²Test(table; direction = Right(), verbose=false, asPearson=false).p ≈ 0.0013797280926100418 ; print("right-dir (Χ²) ✔, ")
    @test Χ²Test(table; direction = Both(), verbose=false, asPearson=false).p ≈ 0.0027594561852200836 ; println("bi-dir (Χ²) ✔, ")

    println("")

    println("1-way ANOVA for Repeated Measures : "); 
    𝐲 = [10., 11., 13., 13., 14., 17., 6., 8., 9.]
    n, k = 3, 3
    ns=(n=n, k=k)
    𝐱=membership(AnovaF_RM(), ns)
    @test fTestRM(𝐲, ns; direction = Both(), verbose=false).p ≈ (1/36); print("bi-dir (F) ✔, ")
    print(" method 2 ")
    𝐲 = [[10., 11., 13.] , [13., 14., 17.], [6., 8., 9.]]
    @test fTestRM(𝐲; direction = Both(), verbose=false).p ≈ (1/36); println("bi-dir (F) ✔, ")
    println("")

    
    println("One-sample t-test ");
    𝐲=[21.5, 24.5, 18.5, 17.2, 14.5, 23.2, 22.1, 20.5, 19.4, 18.1, 24.1, 18.5]
    @test studentTest1S!(𝐲; refmean = 20, direction = Both(), verbose=false).p ≈ 0.84765625 ; print("bi-dir (t) ✔, ")
    𝐲=[21.5, 24.5, 18.5, 17.2, 14.5, 23.2, 22.1, 20.5, 19.4, 18.1, 24.1, 18.5]
    @test studentTest1S!(𝐲; refmean = 20, direction = Right(), verbose=false).p ≈ 0.423828125 ; print("right-dir (t) ✔, ")
    𝐲=[21.5, 24.5, 18.5, 17.2, 14.5, 23.2, 22.1, 20.5, 19.4, 18.1, 24.1, 18.5]
    @test studentTest1S!(𝐲; refmean = 20, direction = Left(), verbose=false).p ≈ 0.581787109375 ; println("left-dir (t) ✔, ")
    𝐲=[21.5, 24.5, 18.5, 17.2, 14.5, 23.2, 22.1, 20.5, 19.4, 18.1, 24.1, 18.5]
    @test studentTest1S(𝐲; refmean = 20, direction = Both(), verbose=false).p ≈ 0.84765625 ; print("bi-dir (t) ✔, ")
    @test studentTest1S(𝐲; refmean = 20, direction = Right(), verbose=false).p ≈ 0.423828125 ; print("right-dir (t) ✔, ")
    @test studentTest1S(𝐲; refmean = 20, direction = Left(), verbose=false).p ≈ 0.581787109375 ; println("left-dir (t) ✔, ")
    println("")


    println("Cochran Test : "); 
    table=[1 1 0; 1 0 0; 1 1 1; 1 1 0; 1 0 1; 1 1 0]
    @test cochranqTest(table; direction = Both(), verbose=false).p ≈ 0.12345679012345678 ; print("bi-dir (Q) ✔, ") 
    
    table=[1 1; 1 0; 1 1; 1 0; 1 0; 1 1; 1 0]
    𝐲, ns = table2vec(table, AnovaF_RM())
    @test cochranqTest(table; direction = Both(), verbose=false).p ≈ 0.125 ; println("bi-dir (Q) ✔, ") 
    println("")
    

    println("Sign test :");
    𝐲 = Bool.([1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0])
    @test signTest(𝐲; direction=Both(), verbose=false).p ≈ 0.454498291015625 ; print("bi-dir (t) ✔, ") 
    @test signTest(𝐲; direction=Right(), verbose=false).p ≈ 0.8949432373046875;  print("right-dir (t) ✔ ")
    @test signTest(𝐲; direction=Left(), verbose=false).p ≈ 0.2272491455078125;  println("left-dir (t) ✔ ")
    println("")
    return true
end

######################################################################################"
# test all functions in multicompTests.jl 

function test_multicompTests(; tol=1e-1)
    println("============ Multiple comparison test functions =================")
         
    # Product-moment Correlation 
    println("")

    println("Product-moment Correlation tests : "); 
    𝐱=[2.0, 4.0, 6.0, 8.0, 10.0]
    𝐘=[[1.0, 3.0, 5.0, 8.0, 7.0], [7.0, 8.0, 5.0, 3.0, 1.0]]  # first hyp is positively and second negatively correlated to 𝐱
    #_permTest!(𝐱, 𝐲, length(𝐱), CrossProd(); fstat=identity)
    @test sum(_permMcTest!(copy(𝐱), 𝐘, length(𝐱), PearsonR(), PearsonR(); verbose=false).p)<=(0.05*length(𝐘)) ;print("bi-directional ✔, ") # both are significant for a two-tailed test
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, length(𝐱), PearsonR(), PearsonR(); fstat=identity, verbose=false).p)==[1] ; print("right-directional ✔, ") # first only significant for a right-tailed test
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, length(𝐱), PearsonR(), PearsonR(); fstat=flip, verbose=false).p)==[2] ; println("left-directional ✔, ") # second only significant for a left-tailed test  
    _permMcTest!(copy(𝐱), 𝐘, length(𝐱), PearsonR(), PearsonR(); fstat=identity, switch2rand=1, nperm=5000, verbose=false); # run the approximate test    
    println("")
    # Use this as an example: the univariate two-dir test gives p-value 0.033 for both hypothesis, hence using Bonferroni none is rejected.
    # Instead the r-max test rejectes them both at the 0.05 level at the first step.

    # 1-Way ANOVA for independent samples 
    ns=[4, 2, 6] # use even numbers here
    n=sum(ns)
    isodd(n) && throw(ArgumentError("test_multicompTests(): 1-way independent samples ANOVA test, n must be even"))
    𝐘=[collect(Base.OneTo(n)).*1.0, rand(n), rand(n)] # only first hypothesis is to be rejected
    𝐱=membership(AnovaF_IS(), ns)
    print("1-way independent samples ANOVA test : "); 
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, ns, AnovaF_IS(), AnovaF_IS(); verbose=false).p)==[1] ;  println(" ✔, ")
    _permMcTest!(copy(𝐱), 𝐘, ns, AnovaF_IS(), AnovaF_IS(); switch2rand=1, nperm=5000, seed=0, verbose=false); # run approximate test

    # T-tests (independent samples)
    ns=[4, 6] # use even numbers
    n=sum(ns)
    isodd(n) && throw(ArgumentError("test_multicompTests(): 1-way independent samples ANOVA test, n must be even"))
    𝐘=[collect(Base.OneTo(n)).*1.0, rand(n), rand(n)] # only first hypothesis is to be rejected
    𝐱=membership(StudentT_IS(), ns)
    println("t-test independent samples test : "); 
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, ns, StudentT_IS(), StudentT_IS(); verbose=false).p)==[1] ;  print(" bi-directional ✔, ")
    @test isempty(findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, ns, StudentT_IS(), StudentT_IS(); fstat=identity, verbose=false).p)) ;  print(" right-directional ✔, ")
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, ns, StudentT_IS(), StudentT_IS(); fstat=flip, verbose=false).p)==[1] ;  println(" left-directional ✔, ")
    _permMcTest!(copy(𝐱), 𝐘, ns, StudentT_IS(), StudentT_IS(); verbose=false, switch2rand=1, nperm=5000); #run approximate test 
    println("")

    # 1-Way repeated-measures ANOVA 
    n, k = 4, 3 # if change this, change also H0 and H1
    H1=[10., 11., 13., 12., 14., 17., 6., 8., 9., 7., 8., 9.]
    H0=[10., 11., 13., 17., 14., 12., 6., 8., 9., 9., 8., 7.]
    𝐘 = [H1, H0, H0] # only the first hypothesis is to be rejected
    ns=(n=n, k=k)
    𝐱=membership(AnovaF_RM(), ns)
    print("1-Way Repeated Measure ANOVA (method 1):");
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), 𝐘, ns, AnovaF_RM(), AnovaF_RM(); verbose=false).p)==[1] ; println(" ✔, ")
    _permMcTest!(copy(𝐱), 𝐘, ns, AnovaF_RM(), AnovaF_RM(); verbose=false, switch2rand=1, nperm=5000); # run approximate test
    println("")

    # One-sample t-test
    ns = 10
    H1=[1.0, 2.0, -0.2, 3.0, 1.5, 0.1, 0.2, 0.3, 0.4, 0.23]
    H0=[-0.1, 0.11, -0.2, 0.18, -0.4, 0.38, -0.16, 0.18, -0.4, 0.41]
    𝐘 = [H0, H1, H0] # only the second hypothesis is to be rejected
    𝐱=membership(StudentT_1S(), ns)
    print("One-sample t-test :");
    @test findall(p->p<=0.05, _permMcTest!(copy(𝐱), copy(𝐘), ns, StudentT_1S(), StudentT_1S(); verbose=false).p)==[2];  println("  ✔ ")
    _permMcTest!(copy(𝐱), copy(𝐘), ns, StudentT_1S(), StudentT_1S(); verbose=false, switch2rand=1, nperm=5000); # just run the approximate test       
    return true
end


function test_multcompTests_API(; tol=1e-1)

    println("============ API for multiple comparison test functions =================")
    println()

    # Product-moment Correlation 
    println("")

    println("Product-moment Correlation tests : "); 
    𝐱=[2.0, 4.0, 6.0, 8.0, 10.0]
    𝐘=[[1.0, 3.0, 5.0, 8.0, 7.0], [7.0, 8.0, 5.0, 3.0, 1.0]]  # first hyp is positively and second negatively correlated to 𝐱
    #_permTest!(𝐱, 𝐲, length(𝐱), CrossProd(); fstat=identity)
    @test sum(rMcTest(𝐱, 𝐘; verbose=false).p)<=(0.05*length(𝐘)) ;print("bi-directional ✔, ") # both are significant for a two-tailed test
    @test findall(p->p<=0.05, rMcTest(𝐱, 𝐘; direction=Right(), verbose=false).p)==[1] ; print("right-directional ✔, ") # first only significant for a right-tailed test
    @test findall(p->p<=0.05, rMcTest(𝐱, 𝐘; direction=Left(), verbose=false).p)==[2] ; println("left-directional ✔, ") # second only significant for a left-tailed test  
    println("")
    # Use this as an example: the univariate two-dir test gives p-value 0.033 for both hypothesis, hence using Bonferroni none is rejected.
    # Instead the r-max test rejectes them both at the 0.05 level at the first step.

    # 1-Way ANOVA for independent samples 
    ns=[4, 2, 6] # use even numbers here
    n=sum(ns)
    isodd(n) && throw(ArgumentError("test_multicompTests(): 1-way independent samples ANOVA test, n must be even"))
    𝐘=[collect(Base.OneTo(n)).*1.0, rand(n), rand(n)] # only first hypothesis is to be rejected
    print("1-way independent samples ANOVA test : "); 
    @test findall(p->p<=0.05, fMcTestIS(𝐘, ns; verbose=false).p)==[1] ;  println(" ✔, ")

    # T-tests (independent samples)
    ns=[4, 6] # use even numbers
    n=sum(ns)
    isodd(n) && throw(ArgumentError("test_multicompTests(): 1-way independent samples ANOVA test, n must be even"))
    𝐘=[collect(Base.OneTo(n)).*1.0, rand(n), rand(n)] # only first hypothesis is to be rejected
    println("t-test independent samples test : "); 
    @test findall(p->p<=0.05, tMcTestIS(𝐘, ns; verbose=false).p)==[1] ;  print("bi-directional ✔, ")
    @test isempty(findall(p->p<=0.05, tMcTestIS(𝐘, ns; direction=Right(), verbose=false).p)) ;  print(" right-directional ✔, ")
    @test findall(p->p<=0.05, tMcTestIS(𝐘, ns; direction=Left(), verbose=false).p)==[1] ;  println(" left-directional ✔, ")
    println("")

    # 1-Way repeated-measures ANOVA 
    n, k = 4, 3 # if change this, change also H0 and H1
    H1=[10., 11., 13., 12., 14., 17., 6., 8., 9., 7., 8., 9.]
    H0=[10., 11., 13., 17., 14., 12., 6., 8., 9., 9., 8., 7.]
    𝐘 = [H1, H0, H0] # only the first hypothesis is to be rejected
    ns=(n=n, k=k)
    print("1-Way Repeated Measure ANOVA (method 1):");
    @test findall(p->p<=0.05, fMcTestRM(𝐘, ns; verbose=false).p)==[1] ; println(" ✔, ")
    println("")

    # One-sample t-test
    ns = 10
    H1=[1.0, 2.0, -0.2, 3.0, 1.5, 0.1, 0.2, 0.3, 0.4, 0.23]
    H0=[-0.1, 0.11, -0.2, 0.18, -0.4, 0.38, -0.16, 0.18, -0.4, 0.41]
    𝐘 = [H0, H1, H0] # only the second hypothesis is to be rejected
    𝐱=membership(StudentT_1S(), ns)
    print("One-sample t-test :");
    @test findall(p->p<=0.05, tMcTest1S(𝐘; verbose=false).p)==[2];  println("  ✔ ")
    return true
end


# test all units of package PermutationTests
@testset "\nAll tests...                 " begin
    @test test_stats()==true
    @test test_unitests()==true
    @test test_uniTests_API()==true
    @test test_multicompTests()==true
    @test test_multcompTests_API()==true
end 
println("\n All tests done.")

