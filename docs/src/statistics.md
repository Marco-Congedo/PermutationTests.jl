# Statistics

Collection of low-level, pure-julia, efficient functions for computing **statistics** and **data transformations**.

These functions are heavely internally by the package. You do not need them
to carry out permutation tests, however they are exported 
because they may turn useful. For general usage of this package you can skip this page altogheter.

## Scalars functions of vectors (statistics)
| Function |  descriptiion  | expression  |
|:--------|:---------------------|:---------------------|
|[`∑`](@ref)(x)           |sum| ``\sum_{i=1}^{N}x_i``|
|[`∑of²`](@ref)(x)        |sum of squares| ``\sum_{i=1}^{N}x_i^2``|
|[`∑∑of²`](@ref)(x)       |sum and sum of squares| ``\sum_{i=1}^{N}x_i`` and ``\sum_{i=1}^{N}x_i^2`` in one pass |
|[`μ`](@ref)(x)           |mean| ``\frac{1}{N} \sum_{i=1}^{N}x_i`` |
|[`dispersion`](@ref)(x)  |sum of squared deviations| ``\sum_{i=1}^{N}(x_i- \mu (x))^2`` |
|[`σ²`](@ref)(x)          |variance| ``\frac{1}{N} \sum_{i=1}^{N}(x_i- \mu (x))^2`` |
|[`σ`](@ref)(x)           |standard deviation | ``\sqrt{\sigma^2(x)}`` |
|[`∑`](@ref)(x, y)        | sum of (x + y) | ``\sum_{i=1}^{N}(x_i+y_i)`` |
|[`Π`](@ref)(c)           | product |``\prod_{i=1}^{N}x_i``|
|[`∑ofΠ`](@ref)(x, y)     | inner product of two vectors| ``x \cdot y = \sum_{i=1}^{N}(x_i \cdot y_i)``|


## Vector functions of vectors (data transformations)
| Function  |  descriptiion  |
|:--------|:---------------------|
|[`μ0`](@ref)(x)      |  recenter (zero mean): ``x_i \gets x_i- \mu (x), \forall i=1 \dots N`` |
|[`σ1`](@ref)(x)      |  equalize (unit standard deviation): ``x_i \gets \frac{x_i}{\sigma (x)}, \forall i=1 \dots N)`` |
|[`μ0σ1`](@ref)(x)    |  standardize (zero mean and unit standard deviation): ``x_i \gets \frac{x_i- \mu (x)}{\sigma (x)}, \forall i=1 \dots N)`` |


```@docs
∑
∑of²
∑∑of²
μ
dispersion
σ²
σ
Π
∑ofΠ
```

---

```@docs
μ0
σ1
μ0σ1
```