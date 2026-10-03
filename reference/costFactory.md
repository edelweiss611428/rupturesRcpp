# `costFactory` class

An `R6` class for fast segment-cost evaluation and parameter estimation,
without a detection algorithm.

## Details

`costFactory` builds the `C++` cost module selected by `$costFunc` at
`$fit()`, and keeps it, so every query reuses its precomputations, as
`PELT`, `binSeg` and `Window` do. `$eval()` and `$get_params()` call the
module directly: segments are `(a, b]` with 0-based `a`, and all checks
are done in `C++`.

`$new()` only stores the `costFunc` object; `$fit()` validates and
attaches the data. `$costFunc` is an active binding, so it can be
inspected or replaced after construction – if data has already been
supplied, replacing it automatically triggers `$fit()` again.

`$get_params()` returns:

- `"L1"`: `median`.

- `"L2"`: `mean`.

- `"SIGMA"`: `mean` and `cov` (plus `epsilon` on the diagonal if
  `addSmallDiag = TRUE`).

- `"VAR"`, `"LinearL2"` and `"LinearL1"`: `coef`, intercept first.

- `"LinearSIGMA"`: `coef` and `cov`.

- `"Custom"`: `params`, the output of `paramFun`; an empty list if
  `paramFun` is `NULL`.

## Methods

- `$new()`:

  Initialises a `costFactory` object.

- `$fit()`:

  Constructs the `C++` cost module.

- `$eval()`:

  Evaluates the cost of a segment.

- `$get_params()`:

  Returns the parameter estimates of a segment.

- `$segments()`:

  Returns the cost and parameter estimates of each segment of a given
  segmentation.

- `$clone()`:

  Clones the `R6` object.

## Author

Minh Long Nguyen <edelweiss611428@gmail.com>  
Huy Nhat Minh Nguyen <sleepysnorlax0115@gmail.com>

## Active bindings

- `costFunc`:

  `R6` object of class `costFunc`. Can be accessed or modified via
  `$costFunc`. Modifying `costFunc` will automatically trigger `$fit()`
  if a `tsMat` has already been fitted.

## Methods

### Public methods

- [`costFactory$new()`](#method-costFactory-new)

- [`costFactory$fit()`](#method-costFactory-fit)

- [`costFactory$eval()`](#method-costFactory-eval)

- [`costFactory$get_params()`](#method-costFactory-get_params)

- [`costFactory$segments()`](#method-costFactory-segments)

- [`costFactory$clone()`](#method-costFactory-clone)

------------------------------------------------------------------------

### Method `new()`

Initialises a `costFactory` object. Does not build the `C++` cost
module; call `$fit()` for that.

#### Usage

    costFactory$new(costFunc)

#### Arguments

- `costFunc`:

  A `R6` object of class `costFunc`. Should be created via
  `costFunc$new()` to avoid error. Default: `costFunc$new("L2")`.

#### Returns

Invisibly returns `NULL`.

------------------------------------------------------------------------

### Method `fit()`

Validates the supplied data and constructs the `C++` cost module
selected by `$costFunc`.

#### Usage

    costFactory$fit(tsMat = NULL, covariates = NULL)

#### Arguments

- `tsMat`:

  Numeric matrix. Time series of size \\n \times p\\. If `NULL`, the
  method will use the previously assigned `tsMat` (i.e., from a prior
  `$fit(tsMat)`). Default: `NULL`.

- `covariates`:

  Numeric matrix with `n` rows, used by `"LinearL2"`, `"LinearSIGMA"`
  and `"LinearL1"`. If `NULL` and no prior `covariates` were set, the
  model is force-fitted with only an intercept. Default: `NULL`.

#### Details

This method constructs the `C++` cost module and sets `private$.fitted`
to `TRUE`, enabling the use of `$eval()` and `$get_params()`.

#### Returns

Invisibly returns `NULL`.

------------------------------------------------------------------------

### Method [`eval()`](https://rdrr.io/r/base/eval.html)

Evaluates the cost of the segment (a, b\] with the module's
[`eval()`](https://rdrr.io/r/base/eval.html).

#### Usage

    costFactory$eval(a, b)

#### Arguments

- `a`:

  Integer. Start index (exclusive, 0-based).

- `b`:

  Integer. End index (inclusive).

#### Returns

The segment cost.

------------------------------------------------------------------------

### Method `get_params()`

Returns the parameter estimates of the segment (a, b\] with the module's
`get_params()`.

#### Usage

    costFactory$get_params(a, b)

#### Arguments

- `a`:

  Integer. Start index (exclusive, 0-based).

- `b`:

  Integer. End index (inclusive).

#### Returns

A named list; see Details.

------------------------------------------------------------------------

### Method [`segments()`](https://rdrr.io/r/graphics/segments.html)

Returns the cost and parameter estimates of each segment of a given
segmentation.

#### Usage

    costFactory$segments(endPts)

#### Arguments

- `endPts`:

  Integer vector. Segment end-points, e.g. from `$predict()` of `PELT`,
  `binSeg`, `Window` or `Dynp`. Sorted internally; must be unique, at
  least `1`, and end at `n`.

#### Returns

A list with one element per segment \\(Start, End\]\\, in the same
format as the `$segments()` of the segmentation classes. Each element is
a list with `Start` (exclusive, 0-based), `End` (inclusive), `Cost` (as
returned by `$eval(Start, End)`) and `Params` (as returned by
`$get_params(Start, End)`).

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    costFactory$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(1)
tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)))
cf = costFactory$new(costFunc$new("L2"))
cf$fit(tsMat)
cf$eval(0, 100)
#> [1] 79.86945
cf$get_params(0, 100)
#> $mean
#> [1] 0.1088874
#> 
cf$segments(c(100, 200))
#> [[1]]
#> [[1]]$Start
#> [1] 0
#> 
#> [[1]]$End
#> [1] 100
#> 
#> [[1]]$Cost
#> [1] 79.86945
#> 
#> [[1]]$Params
#> [[1]]$Params$mean
#> [1] 0.1088874
#> 
#> 
#> 
#> [[2]]
#> [[2]]$Start
#> [1] 100
#> 
#> [[2]]$End
#> [1] 200
#> 
#> [[2]]$Cost
#> [1] 2270.892
#> 
#> [[2]]$Params
#> [[2]]$Params$mean
#> [1] 4.81096
#> 
#> 
#> 

# `costFunc` is an active binding: swapping it re-fits automatically.
cf$costFunc = costFunc$new("SIGMA")
#> `costFunc` has been updated. Re-fitting the model.
cf$eval(0, 100)
#> [1] -22.47755
```
