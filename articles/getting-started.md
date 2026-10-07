# Getting started

This page installs `rupturesRcpp` and runs a first detection, from data
to plot.

## What the package does

`rupturesRcpp` detects change-points in multivariate time series. A
change-point is a time at which the behaviour of a series changes: its
mean, its variance, its autocorrelation, or its relationship with other
variables. Detection is offline: the whole series is available at once,
and the goal is the set of change-points that best explains it.

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

Building from GitHub requires a C++ toolchain: Rtools on Windows; the
Xcode Command Line Tools on macOS (`xcode-select --install`) and a GNU
Fortran compiler (see the [CRAN macOS
tools](https://cran.r-project.org/bin/macosx/tools/) page); or a C++
compiler on Linux (for example, `r-base-dev` on Debian and Ubuntu).

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
