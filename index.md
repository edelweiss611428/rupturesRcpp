# Welcome to rupturesRcpp

[![R-CMD-check](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/edelweiss611428/rupturesRcpp/graphs/commit-activity)
[![rupturesRcpp status
badge](https://edelweiss611428.r-universe.dev/rupturesRcpp/badges/version)](https://edelweiss611428.r-universe.dev/rupturesRcpp)
[![CRAN
Version](https://www.r-pkg.org/badges/version/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp)
[![codecov](https://codecov.io/gh/edelweiss611428/rupturesRcpp/branch/main/graph/badge.svg)](https://app.codecov.io/gh/edelweiss611428/rupturesRcpp)

`rupturesRcpp` is an R package for offline change-point detection in
multivariate time series. Given a whole series, it finds the times at
which its behaviour changes: the mean, the variance, the
autocorrelation, or the relationship with covariates. The algorithms run
in C++, behind an object-oriented R6 interface.

## Features

- Cost functions and search methods implemented in C++ with Rcpp and
  RcppArmadillo.
- One interface for every method: create an object, `$fit()` the data,
  `$predict()` the change-points and `$plot()` the result.
- Change-points chosen by a penalty (`pen`) or by their number
  (`nBkps`), with elbow plots to help choose.
- The cost and parameter estimates of each segment (`$segments()`), or
  of any segment without running a detection (`costFactory`).
- Cost functions written in R, for changes the built-in costs do not
  cover.

## Supported methods

| Cost function | Detects changes in |
|----|----|
| `"L1"`, `"L2"` | the mean (`"L1"` is robust to outliers) |
| `"SIGMA"` | the mean and covariance |
| `"VAR"` | vector autoregressive dynamics |
| `"LinearL2"`, `"LinearL1"` | a linear regression on covariates (`"LinearL1"` is robust to outliers) |
| `"LinearSIGMA"` | a linear regression and its noise covariance |
| `"Custom"` | anything written as an R function |

| Segmentation method | Search |
|----|----|
| `binSeg` | binary segmentation: greedy and fast |
| `Window` | sliding window: local gains, fast |
| `PELT` | optimal for a penalty, with pruning |
| `Dynp` | optimal for a number of change-points, by dynamic programming |

## About the project

`rupturesRcpp` ports the Python library
[ruptures](https://centre-borelli.github.io/ruptures-docs/) to R. It was
created during Google Summer of Code 2025 for The R Project for
Statistical Computing by
[@edelweiss611428](https://github.com/edelweiss611428), with mentors
[@tdhock](https://github.com/tdhock) and
[@deepcharles](https://github.com/deepcharles). The project archive is
on the [gsoc-2025
branch](https://github.com/edelweiss611428/rupturesRcpp/blob/gsoc-2025/README.md).

## Installation

``` r

# Development version, documented on this site
install.packages("rupturesRcpp",
                 repos = c("https://edelweiss611428.r-universe.dev",
                           "https://cloud.r-project.org"))

# CRAN release (1.0.3)
install.packages("rupturesRcpp")
```

[Getting
started](https://edelweiss611428.github.io/rupturesRcpp/articles/getting-started.md)
explains the difference between the two and runs a first detection.
