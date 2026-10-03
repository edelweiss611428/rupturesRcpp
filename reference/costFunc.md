# `costFunc` class

An `R6` class specifying a cost function

## Details

Creates an instance of `costFunc` `R6` class, used in initialisation of
change-point detection modules. Currently supports the following cost
functions:

- **L1 cost function**: \$\$c\_{L_1}(y\_{(a+1)...b}) := \sum\_{t =
  a+1}^{b} \\ y_t - \tilde{y}\_{(a+1)...b} \\\_1\$\$ where
  \\\tilde{y}\_{(a+1)...b}\\ is the coordinate-wise median of the
  segment. If \\a \ge b - 1\\, return 0.

- **L2 cost function**: \$\$c\_{L_2}(y\_{(a+1)...b}) := \sum\_{t =
  a+1}^{b} \\ y_t - \bar{y}\_{(a+1)...b} \\\_2^2\$\$ where
  \\\bar{y}\_{(a+1)...b}\\ is the empirical mean of the segment. If \\a
  \ge b - 1\\, return 0.

- **SIGMA cost function**: \$\$c\_{\sum}(y\_{(a+1)...b}) := (b - a)\log
  \det \hat{\Sigma}\_{(a+1)...b}\$\$ where \\\hat{\Sigma}\_{(a+1)...b}\\
  is the empirical covariance matrix of the segment without Bessel's
  correction. Here, if `addSmallDiag = TRUE`, a small bias `epsilon` is
  added to the diagonal of estimated covariance matrices to improve
  numerical stability.  
    
  By default, `addSmallDiag = TRUE` and `epsilon = 1e-6`. In case
  `addSmallDiag = TRUE`, if the computed determinant of covariance
  matrix is either 0 (singular) or smaller than `p*log(epsilon)` - the
  lower bound, return `(b - a)*p*log(epsilon)`, otherwise, output an
  error message.

- **VAR(r) cost function**: \$\$c\_{\mathrm{VAR}}(y\_{(a+1)...b}) :=
  \sum\_{t = a+r+1}^{b} \left\\ y_t - \sum\_{j=1}^r \hat A_j y\_{t-j}
  \right\\\_2^2\$\$ where \\\hat A_j\\ are the estimated VAR
  coefficients, commonly estimated via the OLS criterion. If system is
  singular, \\a-b \< p\*r+1\\ (i.e., not enough observations), or \\a
  \ge n-p\\ (where `n` is the time series length), return 0.

- **"LinearL2"** for piecewise linear regression process with **constant
  noise variance** \$\$c\_{\text{LinearL2}}(y\_{(a+1):b}) :=
  \sum\_{t=a+1}^b \\ y_t - X_t \hat{\beta} \\\_2^2\$\$ where
  \\\hat{\beta}\\ are OLS estimates on segment \\(a+1):b\\. If segment
  is shorter than the minimum number of points needed for OLS, return 0.

- **"LinearSIGMA"** for piecewise linear regression process with
  **varying noise covariance** \$\$c\_{\text{LinearSIGMA}}(y\_{(a+1):b})
  := (b-a)\log \det \hat\Sigma\_{(a+1):b}\$\$ where
  \\\hat\Sigma\_{(a+1):b}\\ is the empirical covariance matrix of OLS
  residuals \\y - X\hat{\beta}\\ on segment \\(a+1):b\\, estimated the
  same way as in the SIGMA cost function (including the
  `addSmallDiag`/`epsilon` stabilisation and lower-bound fallback).

- **"LinearL1"** for piecewise linear regression process under **L1
  (least absolute deviations) loss**
  \$\$c\_{\text{LinearL1}}(y\_{(a+1):b}) := \sum\_{t=a+1}^b \\ y_t - X_t
  \hat{\beta} \\\_1\$\$ where \\\hat{\beta}\\ is fit column-by-column
  via Iteratively Reweighted Least Squares (IRLS), iterated until the
  change in the fit's cost is within `tol` or `maxIter` iterations are
  reached. Unlike the other regression cost functions, this has no
  \\O(1)\\-per-segment closed form, since IRLS must be re-run on each
  queried segment.

- **"Custom"** for a user-defined cost function supplied from R
  \$\$c\_{\text{Custom}}(y\_{(a+1):b}) := \text{evalFun}(y\_{(a+1):b},
  a, b)\$\$ where `evalFun` is a user-supplied function called on the
  raw segment matrix, plus the segment's own `(a,b]` bounds – the latter
  let `evalFun` align the segment against any externally-captured,
  position-indexed data (e.g. a weight vector or exogenous series
  `evalFun` closes over) without the package needing to know that data
  exists. Because each call crosses back into R, this is substantially
  slower per call than the built-in costs above; see `$evalFun` and
  `$paramFun`.

If active binding `$costFunc` is modified (via assignment operator), the
default parameters will be used.

## Methods

- `$new()`:

  Initialises a `costFunc` object.

- `$pass()`:

  Describes the `costFunc` object.

- `$clone()`:

  Clones the `costFunc` object.

## Author

Minh Long Nguyen <edelweiss611428@gmail.com>

## Active bindings

- `costFunc`:

  Character. Cost function. Can be accessed or modified via `$costFunc`.
  If `costFunc` is modified and required parameters are missing, the
  default parameters are used.

- `pVAR`:

  Integer. Vector autoregressive order. Can be accessed or modified via
  `$pVAR`.

- `addSmallDiag`:

  Logical. Whether to add a bias value to the diagonal of estimated
  covariance matrices to stabilise matrix operations. Can be accessed or
  modified via `$addSmallDiag`.

- `epsilon`:

  Double. A bias value added to the diagonal of estimated covariance
  matrices to stabilise matrix operations. Can be accessed or modified
  via `$epsilon`.

- `intercept`:

  Logical. Whether to include the intercept in regression problems. Can
  be accessed or modified via `$intercept`.

- `tol`:

  Double. IRLS convergence tolerance: iteration stops once the change in
  the fit's cost falls below `tol`. Can be accessed or modified via
  `$tol`.

- `maxIter`:

  Integer. Maximum number of IRLS iterations. Can be accessed or
  modified via `$maxIter`.

- `evalFun`:

  Function. Required for `costFunc = "Custom"`. A user-defined cost
  function, called as `evalFun(segment, a, b)`, where `segment` is the
  numeric matrix of rows `(a+1):b` for the queried segment `(a,b]`
  (0-indexed, same convention as `$eval(a, b)`). `a` and `b` let
  `evalFun` align `segment` against externally-captured,
  position-indexed data it closes over (e.g. `externalSeries[(a+1):b]`),
  which the package itself never needs to see. Must return a single
  numeric value. Can be accessed or modified via `$evalFun`.

- `paramFun`:

  Function or `NULL`. Optional for `costFunc = "Custom"`. A user-defined
  function called as `paramFun(segment, a, b)` (same convention as
  `evalFun`), used by `$get_params()` to report segment-level estimates.
  If `NULL` (default), `$get_params()` returns an empty list for
  `"Custom"`. Can be accessed or modified via `$paramFun`.

## Methods

### Public methods

- [`costFunc$new()`](#method-costFunc-new)

- [`costFunc$pass()`](#method-costFunc-pass)

- [`costFunc$clone()`](#method-costFunc-clone)

------------------------------------------------------------------------

### Method `new()`

Initialises a `costFunc` object.

#### Usage

    costFunc$new(costFunc, ...)

#### Arguments

- `costFunc`:

  Character. Cost function. Supported values include `"L2"`, `"VAR"`,
  and `"SIGMA"`. Default: `L2`.

- `...`:

  Optional named parameters required by specific cost functions.  
  If any required parameters are missing or null, default values will be
  used.

  For `"L1"` and `"L2"`, there is no extra parameter.

  For `"SIGMA"`, supported parameters are:

  `addSmallDiag`

  :   Logical. If `TRUE`, add a small value to the diagonal of estimated
      covariance matrices to stabilise matrix operations. Default:
      `TRUE`.

  `epsilon`

  :   Double. If `addSmallDiag = TRUE`, a small positive value added to
      the diagonal of estimated covariance matrices to stabilise matrix
      operations. Default: `1e-6`.

  For `"VAR"`, `pVAR` is required:

  `pVAR`

  :   Integer. Vector autoregressive order. Must be a positive integer.
      Default: `1L`.

  For `"LinearL2"`, `intercept` is required:

  `intercept`

  :   Logical. Whether to include the intercept in regression problems.
      Default: `TRUE`.

  For `"LinearSIGMA"`, supported parameters are:

  `intercept`

  :   Logical. Whether to include the intercept in regression problems.
      Default: `TRUE`.

  `addSmallDiag`

  :   Logical. If `TRUE`, add a small value to the diagonal of estimated
      residual covariance matrices to stabilise matrix operations.
      Default: `TRUE`.

  `epsilon`

  :   Double. If `addSmallDiag = TRUE`, a small positive value added to
      the diagonal of estimated residual covariance matrices to
      stabilise matrix operations. Default: `1e-6`.

  For `"LinearL1"`, supported parameters are:

  `intercept`

  :   Logical. Whether to include the intercept in regression problems.
      Default: `TRUE`.

  `tol`

  :   Double. IRLS convergence tolerance: iteration stops once the
      change in the fit's cost falls below `tol`. Default: `1e-6`.

  `maxIter`

  :   Integer. Maximum number of IRLS iterations. Default: `1000L`.

  For `"Custom"`, supported parameters are:

  `evalFun`

  :   Function. Required. See `$evalFun` for details.

  `paramFun`

  :   Function or `NULL`. Optional. See `$paramFun` for details.

------------------------------------------------------------------------

### Method `pass()`

Returns a list of configuration parameters to initialise `detection`
modules.

#### Usage

    costFunc$pass()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    costFunc$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## L2 costFunc (default)
costFuncObj = costFunc$new()
costFuncObj$pass()
#> $costFunc
#> [1] "L2"
#> 
## SIGMA costFunc
costFuncObj = costFunc$new(costFunc = "SIGMA")
costFuncObj$pass()
#> $costFunc
#> [1] "SIGMA"
#> 
#> $addSmallDiag
#> [1] TRUE
#> 
#> $epsilon
#> [1] 1e-06
#> 
# Modify active bindings
costFuncObj$epsilon = 10^-5
costFuncObj$pass()
#> $costFunc
#> [1] "SIGMA"
#> 
#> $addSmallDiag
#> [1] TRUE
#> 
#> $epsilon
#> [1] 1e-05
#> 
costFuncObj$costFunc = "VAR"
costFuncObj$pass()
#> $costFunc
#> [1] "VAR"
#> 
#> $pVAR
#> [1] 1
#> 
```
