# Cost functions

``` r

library(rupturesRcpp)
```

This page lists the cost functions built into `rupturesRcpp`, their
options, and how to supply your own cost as an R function. A cost
function measures how badly a single model fits a segment; every
segmentation class and `costFactory` takes one as a `costFunc` object.

## Built-in cost functions

Create a `costFunc` object by specifying the desired (supported) cost
function and optional parameters:

``` r

costFuncObj = costFunc$new("L2")
costFuncObj$pass() #output attributes corresponding to the specified cost function.
#> $costFunc
#> [1] "L2"
```

The following table shows the list of supported cost functions
(pre-implemented ones are `PELT`-compatible; see [Segmentation
classes](https://edelweiss611428.github.io/rupturesRcpp/articles/segmentation-classes.md)).
Here, `n` is segment length.

| **Cost function** | **Description** | **Parameters/active bindings** | **Dimension** | **Time complexity** |
|----|----|----|----|----|
| `"L1"` | Sum of `L1` distances to the segment-wise median; robust to outliers. | `costFunc` | `multi` | `O(nlog(n))` |
| `"L2"` | Sum of squared `L2` distances to the segment-wise mean; faster but less robust than `L1`. | `costFunc` | `multi` | `O(1)` |
| `"SIGMA"` | Log-determinant of the empirical covariance (divided by `n`, no bias correction); models varying mean&variance. | `costFunc`, `addSmallDiag`, `epsilon` | `multi` | `O(1)` |
| `"LinearL1"` | Sum of `L1` residuals from a linear regression model, fit via Iteratively Reweighted Least Squares (IRLS); robust to outliers. | `costFunc`, `intercept`, `tol`, `maxIter` | `multi` | not `O(1)`! |
| `"LinearL2"` | Sum of squared residuals from a linear regression model with constant noise variance. | `costFunc`, `intercept` | `multi` | `O(1)` |
| `"LinearSIGMA"` | Log-determinant of the residual covariance (divided by `n`, no degrees-of-freedom correction) from a linear regression model; models varying noise covariance around a regression mean. | `costFunc`, `intercept`, `addSmallDiag`, `epsilon` | `multi` | `O(1)` |
| `"VAR"` | Sum of squared residuals from a vector autoregressive model with constant noise variance. | `costFunc`, `pVAR` | `multi` | `O(1)` |
| `"Custom"` | User-defined cost, supplied as a plain R function (see *User-defined cost functions* below). | `costFunc`, `evalFun`, `paramFun` | `multi` | depends on `evalFun` |

Every covariance in the package, in the `"SIGMA"`/`"LinearSIGMA"` costs
and in the `cov` returned by `$get_params()`/`$segments()`, is the
biased maximum-likelihood estimate: divided by the segment length `n`,
with no Bessel (`n - 1`) or degrees-of-freedom (`n - q` for `q`
regression coefficients) correction.

If active binding `costFunc` is modified by assigning to
`costFuncObj$costFunc` and the required parameters are missing, the
default parameters will be used. This does not apply to `"Custom"`:
there is no sensible default `evalFun`, so `$pass()`/`$fit()` will error
until one is set (see below).

``` r

costFuncObj$costFunc = "VAR"
costFuncObj$pass()
#> $costFunc
#> [1] "VAR"
#> 
#> $pVAR
#> [1] 1
```

## User-defined cost functions: `costFunc = "Custom"`

If none of the built-in cost functions fit, `costFunc = "Custom"` lets
you supply your own as a plain R function, with no C++ required. It
takes two active bindings:

- `evalFun` (required): called as `evalFun(segment, a, b)`, where
  `segment` is the raw matrix of rows `(a+1):b` of the fitted `tsMat`,
  and `a`/`b` are the same 0-indexed `(a,b]` bounds `$eval(a, b)` uses.
  Must return a single numeric value.
- `paramFun` (optional, default `NULL`): same calling convention, used
  by internal parameter reporting; if omitted, `"Custom"` simply reports
  no parameters.

Passing `a`/`b` through, not just `segment`, is what makes this more
than a convenience wrapper: `evalFun` can use them to align `segment`
against any other externally-captured, position-indexed data (e.g. an
exogenous series or weight vector) that the package itself is never told
about. Because each call crosses back into R, it is substantially slower
per call than the built-in costs, so prefer one of those when it fits.

As a sanity check, re-implementing `"L2"` as a `"Custom"` cost gives
identical numbers:

``` r

myL2eval = function(segment, a, b){
  segment = as.matrix(segment)
  cm = colMeans(segment)
  sum(sweep(segment, 2, cm, FUN = "-")^2)
}

customCF = costFunc$new("Custom", evalFun = myL2eval)
customCF$pass()
#> $costFunc
#> [1] "Custom"
#> 
#> $evalFun
#> function (segment, a, b) 
#> {
#>     segment = as.matrix(segment)
#>     cm = colMeans(segment)
#>     sum(sweep(segment, 2, cm, FUN = "-")^2)
#> }
#> 
#> $paramFun
#> NULL
```

``` r

set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)),
              c(rnorm(100,0), rnorm(100,5,5))) # the series from the 2-regime SIGMA worked example

customObj = PELT$new(minSize = 1L, jump = 1L, costFunc = customCF)
customObj$fit(tsMat)
customObj$eval(0, 150)
#> [1] 3943.78
```

which matches

``` r

L2Obj = PELT$new(costFunc = costFunc$new("L2"))
L2Obj$fit(tsMat)
L2Obj$eval(0, 150)
#> [1] 3943.78
```

exactly.

### Risk of data mismatch

The segmentation logic of existing modules only depends on being able to
compute the cost for an arbitrary segment `(a,b]`; it does not depend on
how the data are stored. Therefore, with a custom cost function, a
mismatch can occur if the function relies on external data that are not
part of the object passed to `$fit()`.

**Data mismatch example**

For example, a custom Poisson cost can silently use externally captured
data that do not match the data passed to `$fit()`:

``` r

set.seed(1)
counts = as.matrix(c(rpois(250, 5), rpois(250, 0)))
counts2 = as.matrix(rpois(500, 5))

poissonCost = function(segment, a, b) {
  y = as.vector(counts[(a + 1):b])
  lambda_hat = mean(y)
  if (lambda_hat <= 0) return(0)
  -2 * sum(dpois(y, lambda_hat, log = TRUE))
}

binSegObj = binSeg$new(
  minSize = 5L,
  costFunc = costFunc$new("Custom", evalFun = poissonCost)
)
binSegObj$fit(counts2) # counts2 has NO change-point by design.
binSegObj$predict(nBkps = 1)
#> [1] 250 500
```

Here, `counts` contains a change-point at 250, but `counts2` does not.
Since `poissonCost` implicitly uses `counts` rather than `counts2`, the
detected segmentation can be inconsistent with the data supplied to
`$fit()`.

``` r

binSegObj$plot()
```

![Count series with no change-point, split at t = 250 because the custom
cost reads a different
series.](cost-functions_files/figure-html/unnamed-chunk-8-1.png)

**Implicit external data example**. The actual use case is a custom cost
function that closes over data the package was never explicitly given.
For example, below, `externalSeries` is captured purely through lexical
scope (it is never passed to `$fit()`), and `evalFun` uses `a` and `b`
to align it with each candidate segment:

``` r

set.seed(1)
tsMat2 = cbind(c(rnorm(100, 0), rnorm(100, 4)))
externalSeries = as.matrix(rnorm(200)) # captured by closure, never passed to `$fit()`

externalRegCost = function(segment, a, b){
  x = externalSeries[(a+1):b, , drop = FALSE]
  sum(lm(segment ~ x)$residuals^2)
}

customObj2 = PELT$new(minSize = 2L, jump = 1L,
                       costFunc = costFunc$new("Custom", evalFun = externalRegCost))
customObj2$fit(tsMat2)
customObj2$predict(pen = 15)
#> [1] 100 200
```

This matches the built-in `"LinearL2"` cost told about `externalSeries`
directly, via `covariates`:

``` r

linObj = PELT$new(minSize = 2L, jump = 1L, costFunc = costFunc$new("LinearL2"))
linObj$fit(tsMat2, externalSeries)
linObj$predict(pen = 15)
#> [1] 100 200
```

`$describe()` reports `evalFun`/`paramFun` as `<function>`/`NULL` rather
than printing the closure itself:

``` r

customObj2$describe(printConfig = TRUE)
#> Pruned Exact Linear Time (PELT) 
#> minSize      : 2L
#> jump         : 1L
#> costFunc     : "Custom"
#> evalFun      : <function>
#> paramFun     : NULL
#> fitted       : TRUE
#> n            : 200L
#> p            : 1L
```

`"Custom"` is supported by `PELT`, `binSeg`, `Window` and `Dynp` alike.
