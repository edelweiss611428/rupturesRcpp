# Welcome to rupturesRcpp

[![R-CMD-check](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/edelweiss611428/rupturesRcpp/graphs/commit-activity) [![rupturesRcpp status badge](https://edelweiss611428.r-universe.dev/rupturesRcpp/badges/version)](https://edelweiss611428.r-universe.dev/rupturesRcpp) [![CRAN Version](https://www.r-pkg.org/badges/version/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/rupturesRcpp)](https://CRAN.R-project.org/package=rupturesRcpp) [![codecov](https://codecov.io/gh/edelweiss611428/rupturesRcpp/branch/main/graph/badge.svg)](https://app.codecov.io/gh/edelweiss611428/rupturesRcpp)

`rupturesRcpp` is an R package for **high-performance offline change-point detection in multivariate time series**. It provides a unified, object-oriented R6 interface to change-point detection methods implemented efficiently in C++. It is inspired by the highly popular Python library [ruptures](https://centre-borelli.github.io/ruptures-docs/).

The package detects **changes in the underlying structure of a time series**, including changes in mean, covariance, autoregressive dynamics, and regression relationships with covariates. It supports a range of cost functions, and segmentation algorithms, from fast heuristic methods to exact dynamic-programming and penalised optimisation methods.


Documentation and case studies are available on the [package website](https://edelweiss611428.github.io/rupturesRcpp/). See [Getting started](https://edelweiss611428.github.io/rupturesRcpp/articles/getting-started.html) for the differences between the two versions and a first detection.

## Supported Methods

### Cost Functions

| Cost function | Detects changes in |
|---------------|--------------------|
| `"L1"`, `"L2"` | Mean structure; `"L1"` provides greater robustness to outliers |
| `"SIGMA"` | Mean and covariance structure |
| `"VAR"` | Vector autoregressive dynamics |
| `"LinearL2"`, `"LinearL1"` | Linear regression relationships with covariates; `"LinearL1"` provides greater robustness to outliers |
| `"LinearSIGMA"` | Regression relationships and residual covariance |
| `"Custom"` | User-defined structural changes through an R cost function |

### Segmentation Methods

| Method | Search strategy |
|---------|-----------------|
| `binSeg` | Binary segmentation; a fast greedy search |
| `Window` | Sliding-window search based on local changes |
| `PELT` | Penalised optimisation with pruning |
| `Dynp` | Dynamic programming for a specified number of change-points |


## Installation

```r
# Development version, documented on the package website
install.packages("rupturesRcpp",
                 repos = c("https://edelweiss611428.r-universe.dev",
                           "https://cloud.r-project.org"))

# CRAN release (2.0.0)
install.packages("rupturesRcpp")
```


## Acknowledgement

The package was developed during **Google Summer of Code 2025** for [The R Project for Statistical Computing](https://www.r-project.org/) by [@edelweiss611428](https://github.com/edelweiss611428), under the mentorship of [@tdhock](https://github.com/tdhock) and [@deepcharles](https://github.com/deepcharles). The original GSoC project archive is available on the [gsoc-2025 branch](https://github.com/edelweiss611428/rupturesRcpp/blob/gsoc-2025/README.md).
