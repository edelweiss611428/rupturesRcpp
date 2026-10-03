# Segment costs and parameters

This chapter shows how to get the cost and parameter estimates of each
segment, either from a fitted segmentation (`$segments()`) or for any
segment you choose, without running a detection algorithm
(`costFactory`).

``` r

library(rupturesRcpp)
```

## Costs and parameters of each segment

After `$predict()`, `$segments()` returns a list with one element per
segment (available on `binSeg`, `Window`, `PELT` and `Dynp`; the example
below uses `PELT`). Each element is a list with `Start`, `End`, `Cost`
and `Params`, where the segment is `(Start, End]` with the same 0-based
convention as `$eval(a, b)`. `Cost` equals `$eval(Start, End)`, and
`Params` is the same named list `costFactory`’s `$get_params()` returns
(see below).

``` r

set.seed(1)
tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)))

PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("SIGMA"))
PELTObj$fit(tsMat)
PELTObj$predict(pen = 50)
#> [1] 100 200
```

``` r

segs = PELTObj$segments()
segs[[2]]
#> $Start
#> [1] 100
#> 
#> $End
#> [1] 200
#> 
#> $Cost
#> [1] 312.2758
#> 
#> $Params
#> $Params$mean
#> [1] 4.81096
#> 
#> $Params$cov
#>          [,1]
#> [1,] 22.70893
```

For `PELT`, the segment costs add up to the optimal penalised cost minus
`pen` times the number of change-points.

``` r

sapply(segs, `[[`, "Cost")
#> [1] -22.47755 312.27581
```

Re-fitting (including through an active binding such as `$minSize` or
`$costFunc`) clears the segmentation saved by the last `$predict()`, so
`$predict()` must be run again before calling `$segments()`.

## Costs without a detection algorithm

Sometimes you don’t need a detection algorithm at all, only fast cost
evaluation and parameter estimation for segments whose boundaries you
already have (e.g., cross-validating a `pen` value against known
change-points, or just querying a segment’s fitted parameters).
`costFactory` wraps the same `C++` cost modules used internally by
`PELT`/`binSeg`/`Window`/`Dynp`, exposing only `$eval()`,
`$get_params()` and `$segments()`, with no segmentation logic.

``` r

set.seed(1)
tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)))

cf = costFactory$new(costFunc$new("L2"))
cf$fit(tsMat)
cf$eval(0, 100)
#> [1] 79.86945
```

``` r

cf$get_params(0, 100)
#> $mean
#> [1] 0.1088874
```

As with the segmentation classes, `$costFunc` is an active binding:
reassigning it after `$fit()` automatically re-fits the underlying
module against the same data.

``` r

cf$costFunc = costFunc$new("SIGMA")
#> `costFunc` has been updated. Re-fitting the model.
```

``` r

cf$eval(0, 100)
#> [1] -22.47755
```

``` r

cf$get_params(0, 100)
#> $mean
#> [1] 0.1088874
#> 
#> $cov
#>           [,1]
#> [1,] 0.7986955
```

`$get_params()`’s return shape depends on the cost function:
`median`/`mean` for `"L1"`/`"L2"`, `mean` and `cov` for `"SIGMA"`,
`coef` (intercept first) for `"VAR"`/`"LinearL2"`/`"LinearL1"`, `coef`
and `cov` for `"LinearSIGMA"`, and `params` (whatever `paramFun`
returns) for `"Custom"`. Both `cov`s are biased estimates, divided by
the segment length `n` with no bias or degrees-of-freedom correction,
and include `epsilon` on the diagonal if `addSmallDiag = TRUE`.

`$segments(endPts)` runs `$eval()` and `$get_params()` over every
segment of a given segmentation, returning the same list as the
segmentation classes’ `$segments()` (see above). `endPts` must end at
`n`, so the output of any `$predict()` can be passed as is:

``` r

segs = cf$segments(c(100, 200))
sapply(segs, `[[`, "Cost")
#> [1] -22.47755 312.27581
```
