# Segmentation methods

This chapter describes the four segmentation classes, the methods they
share, and how refitting works. Each class takes a `costFunc` object and
a time series and returns change-points; they differ in the search they
use, from greedy (`binSeg`, `Window`) to exact (`PELT` for a penalty,
`Dynp` for a number of change-points).

``` r

library(rupturesRcpp)
```

## Classes

After initialising a `costFunc` object (see [Cost
functions](https://edelweiss611428.github.io/rupturesRcpp/articles/cost-functions.md)),
create a segmentation object such as `binSeg`, `Window`, `PELT`, or
`Dynp`.

| **R6 Class** | **Method** | **Description** | **Parameters/active bindings** |
|----|----|----|----|
| `binSeg` | Binary Segmentation | Recursively splits the signal at points that minimise the cost. | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates` |
| `Window` | Slicing Window | Detects change-points using local gains over sliding windows. | `minSize`, `jump`, `radius`, `costFunc`,`tsMat`, `covariates` |
| `PELT` | Pruned Exact Linear Time | Optimal segmentation with pruning for linear-time performance (`costFunc` must be `PELT`-compatible). | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates` |
| `Dynp` | Exact Dynamic Programming | Globally optimal segmentation for a specified number of change-points, via a full dynamic-programming table rather than `binSeg`’s greedy search. | `minSize`, `jump`, `nBkpsMax`, `costFunc`, `tsMat`, `covariates` |

The `covariates` argument is optional and only required for models
involving both dependent and independent variables (e.g., `"LinearL2"`,
`"LinearSIGMA"`, `"LinearL1"`). If not provided, the model is
force-fitted using only an intercept term (i.e., a column of ones).

A `PELT` object, for example, can be initialised as follows:

``` r

detectionObj = PELT$new(minSize = 1L, jump = 1L, costFunc = costFunc$new("L2"))
```

The [case
studies](https://edelweiss611428.github.io/rupturesRcpp/articles/case-studies.md)
run these classes from data to change-points: `binSeg` in [Change in
mean and
variance](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-sigma.md)
and [Change in autoregressive
dynamics](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-var.md),
and `PELT` in [Tuning the
penalty](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-pen.md)
and [Custom cost
functions](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-custom.md).

## Methods

All segmentation objects (`binSeg`, `Window`, `PELT`, `Dynp`) implement
the following methods:

- `$describe(printConfig)`: Views the (current) configurations of the
  object.
- `$fit(tsMat, covariates)`: Constructs a `C++` detection module
  corresponding to the current configurations.
- `$predict(pen, nBkps)`: Performs change-point detection given a linear
  penalty value, or a target number of change-points via `nBkps` (which
  takes precedence over `pen` when both are supplied).
- `$eval(a,b)`: Evaluates the cost of a segment (a,b\].
- `$segments()`: Returns the cost and parameter estimates of each
  segment from the latest `$predict()` (see [Segment costs and
  parameters](https://edelweiss611428.github.io/rupturesRcpp/articles/segment-costs-and-parameters.md)).
- `$plot(d, endPts,...)`: Plots change-point segmentation in `ggplot`
  style.

`binSeg`, `Window`, and `Dynp` additionally implement:

- `$getHistory()`: Returns a `data.frame` of the cost after
  `0, 1, 2, ...` change-points; for `binSeg`/`Window`, also which
  breakpoint was added at each step (see [Model
  selection](https://edelweiss611428.github.io/rupturesRcpp/articles/model-selection.md)
  for why `Dynp`’s version omits that column).
- `$plotElbow(maxK)`: Plots `$getHistory()`’s cost trajectory against
  the number of change-points, for choosing `nBkps` via the “elbow
  method” instead of tuning `pen` directly.

`Dynp` additionally implements `$costPath()`, the raw numeric vector of
exact minimal costs that `$getHistory()` wraps into a `data.frame`.

## Active bindings and refitting

Active bindings (such as `minSize` or `tsMat`) can be modified at any
time, before or after the object is created, via the `$` operator. For
consistency, if the object has already been fitted, modifying any active
bindings will automatically trigger the re-fitting process.

``` r

detectionObj$minSize = 2L #Before fitting
detectionObj$fit(a_time_series_matrix) #Fitted
detectionObj$minSize = 1L #After fitting - automatically trigger `$fit()`
```
