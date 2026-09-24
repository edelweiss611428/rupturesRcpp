# Welcome to `rupturesRcpp`

[![R-CMD-check](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/edelweiss611428/rupturesRcpp/graphs/commit-activity) [![rupturesRcpp status badge](https://edelweiss611428.r-universe.dev/rupturesRcpp/badges/version)](https://edelweiss611428.r-universe.dev/rupturesRcpp)
[![CRAN Version](https://www.r-pkg.org/badges/version/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp) 
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp) [![codecov](https://codecov.io/gh/edelweiss611428/rupturesRcpp/branch/main/graph/badge.svg)](https://app.codecov.io/gh/edelweiss611428/rupturesRcpp)

## Description

<p>The R package provides an efficient, object-oriented R6 interface for offline change point detection, implemented in C++ for high performance. This was created as part of the Google Summer of Code 2025 program (see <a href="https://github.com/edelweiss611428/rupturesRcpp/blob/gsoc-2025/README.md">edelweiss611428/rupturesRcpp at gsoc-2025</a> for the project archive).</p>

<pre>
+------------------------------------------------------------+
|                                                            |
|          Google Summer of Code 2025 Program                |
|                                                            | 
|  Project: rupturesRcpp                                     |
|  Contributor: @edelweiss611428                             |
|  Mentors: @tdhock & @deepcharles                           |
|  Organisation: The R Project for Statistical Computing     |
|                                                            |
+------------------------------------------------------------+
</pre>


## Installation

To install the newest version of the package, use the following R code: 

```r
library(devtools)
install_github("edelweiss611428/rupturesRcpp") 
```

## Getting started

To detect change-points using `rupturesRcpp` you need three main components:

- **Cost function** (`costFunc`)
- **Segmentation method** (`binSeg`, `Window`, `PELT`, `Dynp`)
- **Linear penalty threshold**

Each `component` is implemented using an R6-based object-oriented design for modularity and maintainability.

### Cost functions

Create a `costFunc` object by specifying the desired (supported) cost function and optional parameters:

```r
library("rupturesRcpp")
costFuncObj = costFunc$new("L2")
costFuncObj$pass() #output attributes corresponding to the specified cost function.
```
<pre>
$costFunc
[1] "L2"
</pre>

The following table shows the list of supported cost functions. Here, `n` is segment length.

| **Cost function** | **Description**                                                                                  | **Parameters/active bindings**           | **Dimension** | **Time complexity** |
|-------------------|--------------------------------------------------------------------------------------------------|------------------------------------------|----------------|----------------------|
| `"L1"`            | Sum of `L1` distances to the segment-wise median; robust to outliers.                            | `costFunc`                               | `multi`        | `O(nlog(n))`                 |
| `"L2"`            | Sum of squared `L2` distances to the segment-wise mean; faster but less robust than `L1`.        | `costFunc`                               | `multi`        | `O(1)`                 |
| `"SIGMA"`         | Log-determinant of empirical covariance; models varying mean&variance.                           | `costFunc`, `addSmallDiag`, `epsilon`    | `multi`        | `O(1)`                 |
| `"LinearL1"`      | Sum of `L1` residuals from a linear regression model, fit via Iteratively Reweighted Least Squares (IRLS); robust to outliers. | `costFunc`, `intercept`, `tol`, `maxIter` | `multi`        | not `O(1)`! |
| `"LinearL2"`      | Sum of squared residuals from a linear regression model with constant noise variance.            | `costFunc`, `intercept`                  | `multi`        | `O(1)`                 |
| `"LinearSIGMA"`   | Log-determinant of the residual covariance from a linear regression model; models varying noise covariance around a regression mean. | `costFunc`, `intercept`, `addSmallDiag`, `epsilon` | `multi`        | `O(1)`                 |
| `"VAR"`           | Sum of squared residuals from a vector autoregressive model with constant noise variance.        | `costFunc`, `pVAR`                       | `multi`        | `O(1)`                 |
| `"Custom"`        | User-defined cost, supplied as a plain R function -- see *User-defined cost functions* below.    | `costFunc`, `evalFun`, `paramFun`        | `multi`        | depends on `evalFun`  |

If active binding `costFunc` is modified by assigning to `costFuncObj$costFunc` and the required parameters are missing, the default parameters will be used. This does not apply to `"Custom"`: there is no sensible default `evalFun`, so `$pass()`/`$fit()` will error until one is set (see below).
```r
costFuncObj$costFunc = "VAR"
costFuncObj$pass()
```
<pre>
$costFunc
[1] "VAR"

$pVAR
[1] 1
</pre>


### Segmentation methods

After initialising a `costFunc` object, create a segmentation object such as `binSeg`, `Window`, `PELT`, or `Dynp`.

| **R6 Class**     | **Method**                | **Description**                                                                | **Parameters/active bindings**                                |
|------------------|---------------------------|--------------------------------------------------------------------------------|---------------------------------------------------------------|
| `binSeg`         | Binary Segmentation       | Recursively splits the signal at points that minimise the cost.                | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates`          |
| `Window`         | Slicing Window            | Detects change-points using local gains over sliding windows.                  | `minSize`, `jump`, `radius`, `costFunc`,`tsMat`, `covariates` |
| `PELT`           | Pruned Exact Linear Time  | Optimal segmentation with pruning for linear-time performance (`costFunc` must be `PELT`-compatible).                 | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates`          |
| `Dynp`           | Exact Dynamic Programming | Globally optimal segmentation for a specified number of change-points, via a full dynamic-programming table rather than `binSeg`'s greedy search. | `minSize`, `jump`, `nBkpsMax`, `costFunc`, `tsMat`, `covariates` |

The `covariates` argument is optional and only required for models involving both dependent and independent variables (e.g., `"LinearL2"`, `"LinearSIGMA"`, `"LinearL1"`). If not provided, the model is force-fitted using only 
an intercept term (i.e., a column of ones).

A `PELT` object, for example, can be initialised as follows:
```r
detectionObj = PELT$new(minSize = 1L, jump = 1L, costFunc = costFuncObj)
```

All segmentation objects (`binSeg`, `Window`, `PELT`, `Dynp`) implement the following methods:

- `$describe(printConfig)`: Views the (current) configurations of the object.
- `$fit(tsMat, covariates)`: Constructs a `C++` detection module corresponding to the current configurations.
- `$predict(pen, nBkps)`: Performs change-point detection given a linear penalty value, or a target number of change-points via `nBkps` (which takes precedence over `pen` when both are supplied).
- `$eval(a,b)`: Evaluates the cost of a segment (a,b].
- `$plot(d, endPts,...)`: Plots change-point segmentation in `ggplot` style.

`binSeg`, `Window`, and `Dynp` additionally implement:

- `$getHistory()`: Returns a `data.frame` of the cost after `0, 1, 2, ...` change-points -- for `binSeg`/`Window`, also which breakpoint was added at each step (see *Exact optimal segmentation via `Dynp`* below for why `Dynp`'s version omits that column).
- `$plotElbow(maxK)`: Plots `$getHistory()`'s cost trajectory against the number of change-points, for choosing `nBkps` via the "elbow method" instead of tuning `pen` directly.

`Dynp` additionally implements `$costPath()`, the raw numeric vector of exact minimal costs that `$getHistory()` wraps into a `data.frame`.

Active bindings (such as `minSize` or `tsMat`) can be modified at any time—either before or after the object is created via the `$` operator. 
For consistency, if the object has already been fitted, modifying any active bindings will automatically trigger the re-fitting process.

```r
detectionObj$minSize = 2L #Before fitting
detectionObj$fit(a_time_series_matrix) #Fitted
detectionObj$minSize = 1L #After fitting - automatically trigger `$fit()`
```

## Simulated data examples

### 2-regime SIGMA example via binary segmentation

To demonstrate the package usage, we first consider a simple 2d time series with two piecewise Gaussian regimes and varying variance.

```r
set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)),
              c(rnorm(100,0), rnorm(100,5,5)))
```
<img width="2492" height="872" alt="image" src="https://github.com/user-attachments/assets/65b5511c-070e-4b2d-872b-410679b4e395" />


As our example involves regimes with varying variance, a suitable `costFunc` option is `"SIGMA"`.  Since the segmentation objects' interfaces are similar, it is sufficient to demonstrate the usage of `binSeg` only.

```r
SIGMAObj = costFunc$new("SIGMA", addSmallDiag = TRUE, epsilon = 1e-6)
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = SIGMAObj) 
binSegObj$fit(tsMat) 
```

Once fitted, `$predict()` and `$eval()` can be used. To view the configurations of the `binSeg` object, we can use `$describe()`. 

```r
binSegObj$describe(printConfig = TRUE) 
```
<pre>
Binary Segmentation (binSeg)
minSize      : 1L
jump         : 1L
costFunc     : "SIGMA"
addSmallDiag : TRUE
epsilon      : 1e-06
fitted       : TRUE
n            : 200L
p            : 2L
</pre>

To obtain an estimated segmentation, we can use the `$predict()` method and specify a non-negative penalty value `pen`, which should be properly tuned. This returns a sorted integer vector of end-points, including the number of observations by design. 

Here, we set `pen = 100`.

```r
binSegObj$predict(pen = 100)
```
<pre>
[1] 100 200
</pre>

After running `$predict()`, the segmentation output is temporarily saved to the `binSeg` object, allowing users to use the `$plot()` method without specifying `endPts`.

```r
binSegObj$plot(d = 1:2, 
               main = "method: binSeg; costFunc: SIGMA; pen: 100")
```
<img width="2492" height="872" alt="image" src="https://github.com/user-attachments/assets/f8750edf-13d8-4363-b158-9beb744bef0b" />


### 2-regime VAR example: Modifying a `binSeg` object through its active bindings

You can also modify a `binSeg` object's fields through its active bindings. To demonstrate this, we consider a piecewise vector autoregressive example with constant noise variance.

```r
set.seed(1)
tsMat = matrix(c(filter(rnorm(100), filter = 0.9, method = "recursive"), 
                 filter(rnorm(100), filter = -0.9, method = "recursive")))
```
<img width="2492" height="872" alt="image" src="https://github.com/user-attachments/assets/1837aac4-d37d-4835-b124-abdd6b74dbc6" />


Here, the most suitable cost function is `"VAR"`. Instead of creating a new `binSeg` object, we will modify the current `binSegObj` as follows:

```r
VARObj = costFunc$new("VAR")
binSegObj$tsMat = tsMat
binSegObj$costFunc = VARObj
```
Modifying `tsMat` (or any other bindings) will automatically trigger `self$fit()` if a `tsMat` has already existed. 

```r
binSegObj$describe(printConfig = TRUE)
```

<pre>
Binary Segmentation (binSeg)
minSize      : 1L
jump         : 1L
costFunc     : "VAR"
pVAR         : 1L
fitted       : TRUE
n            : 200L
p            : 1L
</pre>

We can then perform binary segmentation with `pen = 25` and plot the segmentation results.

```r
binSegObj$predict(pen = 25)
binSegObj$plot(d = 1L, 
               main = "method: binSeg; costFunc: VAR; pen: 25")
```
<img width="2492" height="872" alt="image" src="https://github.com/user-attachments/assets/f677f835-1a99-41b3-a244-6b4e5de25f93" />


### User-defined cost functions: `costFunc = "Custom"`

If none of the built-in cost functions fit, `costFunc = "Custom"` lets you supply
your own as a plain R function -- no C++ required. It takes two active bindings:

- `evalFun` (required): called as `evalFun(segment, a, b)`, where `segment` is the
  raw matrix of rows `(a+1):b` of the fitted `tsMat`, and `a`/`b` are the same
  0-indexed `(a,b]` bounds `$eval(a, b)` uses. Must return a single numeric value.
- `paramFun` (optional, default `NULL`): same calling convention, used by internal
  parameter reporting; if omitted, `"Custom"` simply reports no parameters.

Passing `a`/`b` through -- not just `segment` -- is what makes this more than a
convenience wrapper: `evalFun` can use them to align `segment` against any other
externally-captured, position-indexed data (e.g. an exogenous series or weight
vector) that the package itself is never told about. Because each call crosses
back into R, it is substantially slower per call than the built-in costs -- prefer
one of those when it fits.

As a sanity check, re-implementing `"L2"` as a `"Custom"` cost gives identical numbers:

```r
myL2eval = function(segment, a, b){
  segment = as.matrix(segment)
  cm = colMeans(segment)
  sum(sweep(segment, 2, cm, FUN = "-")^2)
}

customCF = costFunc$new("Custom", evalFun = myL2eval)
customCF$pass()
```
<pre>
$costFunc
[1] "Custom"

$evalFun
function (segment, a, b) 
{
    segment = as.matrix(segment)
    cm = colMeans(segment)
    sum(sweep(segment, 2, cm, FUN = "-")^2)
}

$paramFun
NULL
</pre>

```r
customObj = PELT$new(minSize = 1L, jump = 1L, costFunc = customCF)
customObj$fit(tsMat) # tsMat from the 2-regime SIGMA example above
customObj$eval(0, 150)
```
<pre>
[1] 3943.78
</pre>
which matches 
```r
PELT$new(costFunc = costFunc$new("L2"))$fit(tsMat)$eval(0, 150)
```
<pre>
[1] 3943.78
</pre>
exactly.

#### Risk of data mismatch

The segmentation logic of existing modules only depends on being able to compute the cost for an arbitrary segment \((a,b]\); it does not depend on how the data are stored. Therefore, with a custom cost function, a mismatch can occur if the function relies on external data that are not part of the object passed to `$fit()`.

**Data mismatch example**

For example, a custom Poisson cost can silently use externally captured data that do not match the data passed to `$fit()`:

```r
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
```

<pre>
[1] 250 500
</pre>

Here, `counts` contains a change-point at 250, but `counts2` does not. Since `poissonCost` implicitly uses `counts` rather than `counts2`, the detected segmentation can be inconsistent with the data supplied to `$fit()`.

```{r}

binSegObj$plot()
```
<img width="2516" height="1403" alt="image" src="https://github.com/user-attachments/assets/c85c24fa-2742-4cb8-9563-4f47337ba9dd" />


**Implicit external data example**. The actual use case is a custom cost function that closes over data the package was never explicitly given. For example, below, `externalSeries` is captured purely through lexical scope—it is never passed to `$fit()`—and `evalFun` uses `a` and `b` to align it with each candidate segment:

```r
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
```
<pre>
[1] 100 200
</pre>

This matches the built-in `"LinearL2"` cost told about `externalSeries` directly, via `covariates`:

```r
linObj = PELT$new(minSize = 2L, jump = 1L, costFunc = costFunc$new("LinearL2"))
linObj$fit(tsMat2, externalSeries)
linObj$predict(pen = 15)
```
<pre>
[1] 100 200
</pre>

`$describe()` reports `evalFun`/`paramFun` as `<function>`/`NULL` rather than printing the closure itself:

```r
customObj2$describe(printConfig = TRUE)
```
```text
Pruned Exact Linear Time (PELT) 
minSize      : 2L
jump         : 1L
costFunc     : "Custom"
evalFun      : <function>
paramFun     : NULL
fitted       : TRUE
n            : 200L
p            : 1L
```

`"Custom"` is supported by `PELT`, `binSeg`, and `Window` alike.

### Elbow-method model selection: `$getHistory()`, `$plotElbow()`, and `$predict(nBkps = ...)`

`binSeg` and `Window` both build up their segmentation by adding one change-point
at a time -- `binSeg` by recursively splitting the segment that most reduces
cost, `Window` by ranking candidate local maxima by gain. `$getHistory()`
exposes that trajectory directly, so you can inspect it, or choose the number
of change-points via the "elbow method", instead of only tuning `pen`.

Continuing with the `binSegObj` (`"VAR"` cost) from the previous example:

```r
binSegObj$getHistory()
```
<pre>
   k     cost added_bkp
1  0 533.9504        NA
2  1 165.2573        99
3  2 159.6331       111
4  3 154.7607       159
5  4 149.1974       180
...
</pre>

`$plotElbow()` renders this as a `ggplot` object (cost vs. number of
change-points); look for where the marginal decrease in cost flattens out to
pick `k`.

```r
binSegObj$plotElbow()
```

Once a `k` is chosen, `$predict()` accepts it directly via `nBkps`, which
takes precedence over `pen` when both are supplied:

```r
binSegObj$predict(nBkps = 1)
```
<pre>
[1]  99 200
</pre>

`nBkps` is treated as an upper bound, not a strict requirement: for `binSeg`
it returns its own best answer among the splits it already explored, and for
`Window` the `nBkps` highest-gain local maxima it found -- neither is
guaranteed to be the *globally* optimal segmentation for that count. `Window`
in particular can only ever offer as many change-points as it found local
maxima for; if you ask for more, `$predict()` returns what's available and
reports the shortfall via a message rather than erroring.

### Exact optimal segmentation via `Dynp`

The caveat above -- `binSeg`/`Window` returning their own best answer for a given `nBkps`, not necessarily the *globally* optimal one -- is exactly what `Dynp` addresses. It finds the segmentation that is exactly optimal for a specified number of change-points, by building a full dynamic-programming table instead of searching greedily. The cost is complexity: the table takes `O(nBkpsMax * M^2)` time, where `M` is the number of `(minSize, jump)`-admissible positions, versus `binSeg`'s near-`O(n log n)` greedy search or `PELT`'s pruned search for a penalty.

```r
set.seed(1121)
signals = as.matrix(c(rnorm(100, 0, 1), rnorm(100, 5, 1)))

DynpObj = Dynp$new(minSize = 1L, jump = 1L) # nBkpsMax left NULL -> resolved at $fit()
DynpObj$fit(signals)
```
<pre>
`nBkpsMax` not set; using 20L.
</pre>

`$describe()` shows the two fields specific to `Dynp`: the user-set `nBkpsMax` (here left `NULL`) and the `resolvedNBkpsMax` actually used, capped at 20 by default.

```r
DynpObj$describe(printConfig = TRUE)
```
<pre>
Exact Dynamic Programming (Dynp) 
minSize          : 1L
jump             : 1L
nBkpsMax         : NULL (auto)
resolvedNBkpsMax : 20L
costFunc         : "L2"
fitted           : TRUE
n                : 200L
p                : 1L
</pre>

`$predict()`, `$eval()`, and `$plot()` all work exactly as with `binSeg`/`Window`/`PELT`, either via `nBkps` or via a penalty:

```r
DynpObj$predict(nBkps = 1)
```
<pre>
[1] 100 200
</pre>

```r
DynpObj$predict(pen = 100) # same split, reached via a penalty instead
```
<pre>
[1] 100 200
</pre>

```r
DynpObj$plot(main = "method: Dynp; costFunc: L2; nBkps: 1")
```

Since `$fit()` already computes the exact minimal cost for *every* change-point count from `0` to `nBkpsMax` in one pass (exposed via `$costPath()`), `$getHistory()`/`$plotElbow()` come for free here -- no extra search is needed the way `binSeg`/`Window` would need to keep growing their own nested sequence:

```r
DynpObj$getHistory()
```
<pre>
  k      cost
1 0 1511.6991
2 1  186.6857
3 2  176.1339
4 3  169.4885
5 4  163.1693
6 5  158.4455
...
</pre>

Unlike `binSeg`/`Window`, there is no `added_bkp` column: each `k`'s solution is independently exact and need not be nested inside the solution for `k+1`, so "the one breakpoint added at this step" is not generally well-defined. Use `$predict(nBkps = k)` to get the full breakpoint set for a given `k`.

## Cost evaluation without detection: `costFactory`

Sometimes you don't need a detection algorithm at all -- only fast cost evaluation and parameter estimation for segments whose boundaries you already have (e.g., cross-validating a `pen` value against known change-points, or just querying a segment's fitted parameters). `costFactory` wraps the same `C++` cost modules used internally by `PELT`/`binSeg`/`Window`/`Dynp`, exposing only `$eval()` and `$get_params()` -- no segmentation logic.

```r
set.seed(1)
tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)))

cf = costFactory$new(costFunc$new("L2"))
cf$fit(tsMat)
cf$eval(0, 100)
```
<pre>
[1] 79.86945
</pre>

```r
cf$get_params(0, 100)
```
<pre>
$mean
[1] 0.1088874
</pre>

As with the segmentation classes, `$costFunc` is an active binding: reassigning it after `$fit()` automatically re-fits the underlying module against the same data.

```r
cf$costFunc = costFunc$new("SIGMA")
```
<pre>
`costFunc` has been updated. Re-fitting the model.
</pre>

```r
cf$eval(0, 100)
```
<pre>
[1] -22.47755
</pre>

```r
cf$get_params(0, 100)
```
<pre>
$mean
[1] 0.1088874

$cov
          [,1]
[1,] 0.7986955
</pre>

`$get_params()`'s return shape depends on the cost function: `median`/`mean` for `"L1"`/`"L2"`, `mean` and `cov` for `"SIGMA"`, `coef` (intercept first) for `"VAR"`/`"LinearL2"`/`"LinearL1"`, `coef` and `cov` for `"LinearSIGMA"`, and `params` (whatever `paramFun` returns) for `"Custom"`.

## Future development

- Improve the `"L1"` cost module, potentially allowing queries in `O(log(n))` time using data structures such as a persistent segment tree with `O(nlog(n))` precomputation.
- Clean and enhance the existing object-oriented interface for improved efficiency, robustness, and accessibility (see https://github.com/edelweiss611428/R6BinSeg/tree/main for an idea).
- Implement additional cost functions (e.g., `"Poisson"`). 
- Implement other offline change-point detection classes (e.g., `BottomUp`).
- Improve `$plot()` method for models involving both dependent and independent variables.

## Contributing

We welcome all contributions, whether big or small. If you encounter a bug or have a feature request, please open an issue to let us know. 

Feel free to fork the repository and make your changes. For significant updates, it’s best to discuss them with us first. When your changes are ready, submit a pull request.

Thanks for helping us improve this project!

## License

This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) License. 


## References

- Hocking, T. D. (2024). *Finite Sample Complexity Analysis of Binary Segmentation*. arXiv preprint arXiv:2410.08654. 
- Truong, C., Oudre, L., & Vayatis, N. (2020). *Selective review of offline change point detection methods*. Signal Processing, 167, 107299. 
- Killick, R., Fearnhead, P., & Eckley, I. A. (2012). *Optimal detection of change points with a linear computational cost*. Journal of the American Statistical Association, 107(500), 1590–1598. 



