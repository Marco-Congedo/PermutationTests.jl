| **Documentation**  |
|:---------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Marco-Congedo.github.io/PermutationTests.jl/dev) |

**PermutationTests.jl** is a pure-[**Julia**](https://julialang.org/) comprehensive, fast and well-documented
package for performing *univariate* and *multiple comparisons* statistical hypothesis tests based
on *permutation theory*.

---
## Installation

Execute the following command in Julia's REPL:

    ]add PermutationTests

---
## Available tests

For each *univariate* test there is its *multiple comparisons* counterpart: 
- Pearson product-moment correlation
- Trend correlation (fit of any kind of regression)
- Point bi-serial correlation*
- Student's t for independent samples
- 1-way ANOVA for independent samples
- Χ² for 2xK contingency tables*
- Fisher exact test* (2x2 contingency tables)
- Student's t for repeated-measures 
- 1-way ANOVA for repeated-measures
- Cochran Q*
- McNemar*
- One-sample Student's t  
- Sign*

* for dicothomous data 

---
## Quick start

As an example, let's run a Pearson correlation univariate test:

```
using PermutationTests
N=10 # number of observations
x, y = randn(N), randn(N) # some random Gaussian data
t = rTest(x, y)
```

The test result `t` is a structure and its fields are printed in yellow, 
for example:

![](/docs/src/assets/Result_example.png)


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
| **Documentation**  | 
|:---------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Marco-Congedo.github.io/PermutationTests.jl/dev) |
