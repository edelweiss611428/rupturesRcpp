# Getting started

This page installs `rupturesRcpp` and runs a first detection, from data
to plot.

## What the package does

`rupturesRcpp` detects change-points in multivariate time series. A
change-point is a time at which the behaviour of a series changes: its
mean, its variance, its autocorrelation, or its relationship with other
variables. Detection is offline: the whole series is available at once,
and the goal is the set of change-points that best explains it.

The cost functions and search methods follow the Python library
[ruptures](https://centre-borelli.github.io/ruptures-docs/), whose
methods are reviewed by [Truong, Oudre and Vayatis
(2020)](https://doi.org/10.1016/j.sigpro.2019.107299). `rupturesRcpp`
implements them in C++ with Rcpp and RcppArmadillo, wraps them in R6
classes, and also accepts cost functions written in R. The package was
created during Google Summer of Code 2025 for the R Project for
Statistical Computing.

## Installation

Choose how to install. This site documents the development version.

- Development version (recommended)
- From source
- CRAN release

r-universe builds the development version from the `main` branch on
GitHub, with binaries for Windows, macOS and Linux, so no compiler is
needed.

``` r

install.packages("rupturesRcpp",
                 repos = c("https://edelweiss611428.r-universe.dev",
                           "https://cloud.r-project.org"))
```

Building from GitHub needs a C++ toolchain: Rtools on Windows, the Xcode
command line tools on macOS (`xcode-select --install`), or a C++
compiler on Linux (for example `r-base-dev` on Debian and Ubuntu).

``` r

# install.packages("remotes")
remotes::install_github("edelweiss611428/rupturesRcpp")
```

``` r

install.packages("rupturesRcpp")
```

The CRAN release (1.0.3) has `binSeg`, `Window` and `PELT` with the
`"L1"`, `"L2"`, `"SIGMA"`, `"VAR"` and `"LinearL2"` costs, but not the
rest of what this site uses:

- `Dynp` and `costFactory`;
- the `"LinearSIGMA"`, `"LinearL1"` and `"Custom"` costs;
- `$predict(nBkps = ...)`, `$getHistory()`, `$plotElbow()` and
  `$segments()`.

## A first detection

Load the package and check that the version is 2.0.0 or later:

``` r

library(rupturesRcpp)
packageVersion("rupturesRcpp")
#> [1] '2.0.0'
```

Both series below change in mean and variance at t = 100. `tsMat` has
one row per time point and one column per feature.

``` r

set.seed(1)
tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)),
              c(rnorm(100, 0), rnorm(100, 5, 5)))
```

The variance changes as well as the mean, so the `"SIGMA"` cost, which
models both, suits this data. Binary segmentation (`binSeg`) is a fast
search method.

``` r

binSegObj = binSeg$new(costFunc = costFunc$new("SIGMA"))
binSegObj$fit(tsMat)
binSegObj$predict(pen = 100)
#> [1] 100 200
```

`$predict()` returns the end of each segment, so the last value is
always the number of observations. Here it finds the one change-point at
t = 100. The penalty `pen` sets how many change-points are kept.

``` r

binSegObj$plot(d = 1:2)
```

![Two simulated series split into two shaded segments at the detected
change-point, t =
100.](getting-started_files/figure-html/unnamed-chunk-7-1.png)

## Next steps

The
[Documentation](https://edelweiss611428.github.io/rupturesRcpp/articles/documentation.md)
explains each part of a detection: the cost function, the segmentation
method and how many change-points to keep. [Case
studies](https://edelweiss611428.github.io/rupturesRcpp/articles/case-studies.md)
has longer worked examples.
