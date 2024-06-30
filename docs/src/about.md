# About

---
## About the code

**PermutationTests.jl** includes seven code units (.jl files):

| Unit   | Description |
|:----------|:----------|
| PermutationTests.jl | Main module, declaration of constants, types and structuress |
| stats.jl | Low-level computations for statistics and test-statistic |
| tools.jl | General tools |
| uniTests.jl | Univariate tests |
| uniTests_API.jl | API for univariate tests |
| multcompTests.jl | Multiple comparisons tests |
| multcompTests_API.jl | API for multiple comparisons tests |

In addition, units for running [benchmarks](@ref "Benchmarks"), [error control](@ref "Error Control") tests and [power](@ref "Power") analysis can be found
in the [src\extras](https://github.com/Marco-Congedo/PermutationTests.jl/tree/master/src/extras) folder.

A unit to test the main functions is available as well in the 
[test](https://github.com/Marco-Congedo/PermutationTests.jl/tree/master/test) folder.

---
## About the authors

[Marco Congedo](https://sites.google.com/site/marcocongedo), corresponding author and developer of the package, is a Research Director of [CNRS](http://www.cnrs.fr/en) (Centre National de la Recherche Scientifique), working at
[UGA](https://www.univ-grenoble-alpes.fr/english/) (University of Grenoble Alpes, France).
**Contact**: first name dot last name at gmail dot com

[Livio Finos](https://pnc.unipd.it/finos-livio/), is Full professor at the  [Department of Statistical Sciences](https://www.unipd.it/en/stat) of [Univerità di Padova, Italy](https://pnc.unipd.it/).
**Contact**: first name dot last name at unipd dot it

---
## Disclaimer

This version has been roughly tested.
Independent reviewers for both the code and the documentation are welcome.

---
## TroubleShoothing

| Problem   | Solution |
|:----------|:----------|
| [Folds.jl](https://github.com/JuliaFolds/Folds.jl) does not work properly for future versions of julia | use keyword argument `threaded=false` for multiple comparison tests |

---
## References

[R.A. Fisher](https://en.wikipedia.org/wiki/Ronald_Fisher) (1935) The Design of Experiments, Hafner.

E.S. Edgington (1995), Randomization Tests, Marcel Dekker Inc.

A.P. Holmes, R.C. Blair, J.D.G Watson, I. Ford (1996) Non-Parametric Analysis of Statistic Images From Functional Mapping Experiments. Journal of Cerebral Blood Flow and Metabolism 16:7-22

F. Pesarin (2001) Multivariate Permutation Tests with applications in Biostatistics. John Wiley & Sons.

E. J. G. Pitman (1937) Significance tests which may be applied to samples from any population, Royal Statistical Society Supplement, 4: 119-130 and 225-32 (parts I and II). 

[E. J. G. Pitman](https://en.wikipedia.org/wiki/E._J._G._Pitman) (1938) Significance tests which may be applied to samples from any population. Part III. The analysis of variance test. Biometrika. 29 (3–4): 322–335.

P.H. Westfall, S.S. Young (1993) Resampling-Based Multiple Testing: Examples and Methods for P-Value Adjustment, John Wiley & Sons.

---
## Contents

```@contents
    Pages =  [
        "index.md",
        "about.md",
        "PermutationTests.md",
        "univariate tests.md",
        "multiple comparisons tests.md",
        "package tests.md",
        "tools.md",
        "statistics.md",
        "test statistics.md",
        "create your own test.md",
        "chose a test.md"
    ]
Depth = 2
```

---
## Index

```@index
```

