# Welcome to rupturesRcpp
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/edelweiss611428/rupturesRcpp/graphs/commit-activity)
## Description

<p align="justify"> The R package provides an efficient, object-oriented R6 interface for offline change point detection, implemented in C++ for high performance. This was created as part of the Google Summer of Code 2025 program. </p>


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

- **Cost function**(`costFunc`)
- **Segmentation method**(`binSeg`, `Window`, `PELT`)
- **Linear penalty threshold**

Each `component` is implemented using an R6-based object-oriented design for clarity and flexibility.

### Cost functions

Create a `costFunc` object by specifying the desired (supported) cost function and optional parameters:

```r
library("rupturesRcpp")
costFuncObj = costFunc$new("L2")
```
The following table shows the list of supported cost functions.

| **Cost function** | **Description**                                                                                  | **Parameters/active bindings**           | **Dimension** | **Time complexity** |
|-------------------|--------------------------------------------------------------------------------------------------|------------------------------------------|----------------|----------------------|
| `"L1"`            | Sum of `L1` distances to the segment-wise median; robust to outliers.                            | `costFunc`                               | `multi`        | O(n)                 |
| `"L2"`            | Sum of squared `L2` distances to the segment-wise mean; faster but less robust than `L1`.        | `costFunc`                               | `multi`        | O(1)                 |
| `"SIGMA"`         | Log-determinant of empirical covariance; models varying mean&variance.                           | `costFunc`, `addSmallDiag`, `epsilon`    | `multi`        | O(1)                 |
| `"VAR"`           | Sum of squared residuals from a vector autoregressive model with constant noise variance.        | `costFunc`, `pVAR`                       | `multi`        | O(1)                 |
| `"LinearL2"`      | Sum of squared residuals from a linear regression model with constant noise variance.            | `costFunc`, `intercept`                  | `multi`        | O(1)                 |

### Segmentation methods

After initialising a `costFunc` object, create a segmentation object such as `binSeg`, `Window`, or `PELT`.

| **R6 Class**     | **Method**                | **Description**                                                                | **Parameters/active bindings**                                |
|------------------|---------------------------|--------------------------------------------------------------------------------|---------------------------------------------------------------|
| `binSeg`         | Binary Segmentation       | Recursively splits the signal at points that minimise the cost.                | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates`          |
| `Window`         | Window Slicing            | Detects change-points using local gains over sliding windows.                  | `minSize`, `jump`, `radius`, `costFunc`,`tsMat`, `covariates` |
| `PELT`           | Pruned Exact Linear Time  | Optimal segmentation with pruning for linear-time performance.                 | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates`          |

The `covariates` argument is optional and only required for models involving both dependent and independent variables (e.g., `"LinearL2"`). If not provided, the model is force-fitted using only 
an intercept term (i.e., a column of ones).

A `PELT` object, for example, can be initialised as follows:
```r
detectionObj = PELT$new(minSize = 1L, jump = 1L, costFunc = costFuncObj)
```

All segmentation objects implement the following methods:

- `$describe(printConfig)`: Views the (current) configurations of the object.
- `$fit(tsMat, covariates)`: Constructs a `C++` detection module corresponding to the current configurations.
- `$predict(pen)`: Performs change-point detection given a linear penalty value.
- `$eval(a,b)`: Evaluates the cost of a segment (a,b].
- `$plot(d, endPts,...)`: Plots change-point segmentation in `ggplot` style.

Active bindings (such as `minSize` or `tsMat`) can be modified at any time—either before or after the object is created via the `$` operator. 
For consistency, if the object has already been fitted, modifying any active bindings will automatically re-trigger the fitting process.

```r
detectionObj$minSize = 2L #Before fitting
detectionObj$fit(a_time_series_matrix) #Fitted
detectionObj$minSize = 1L #After fitting - automatically re-trigger `$fit()`
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
<img width="2492" height="872" alt="image" src="https://github.com/user-attachments/assets/57a7fb68-b5cb-4189-92da-a76b36d0d8b4" />


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


## Future development

- Provide additional cost functions (e.g., `"Poisson"`, `"Linear-L1"`, and `"Linear-L2"`). 
- Implement other offline change-point detection classes (e.g., `Opt` and `Dynp`).
- Enhance the existing object-oriented interface for improved efficiency, robustness, and accessibility.
- Improve `$plot()` method for models involving both dependent and independent variables.
- Provide instructions for future contributors.
- (Tentative) Implement online change-point detection classes.

## References

- Hocking, T. D. (2024). *Finite Sample Complexity Analysis of Binary Segmentation*. arXiv preprint arXiv:2410.08654. [https://arxiv.org/abs/2410.08654](https://arxiv.org/abs/2410.08654)

- Truong, C., Oudre, L., & Vayatis, N. (2020). *Selective review of offline change point detection methods*. Signal Processing, 167, 107299. [https://doi.org/10.1016/j.sigpro.2019.107299
](https://www.sciencedirect.com/science/article/abs/pii/S0165168419303494?via%3Dihub#:~:text=https%3A//doi.org/10.1016/j.sigpro.2019.107299)

- Killick, R., Fearnhead, P., & Eckley, I. A. (2012). *Optimal detection of change points with a linear computational cost*. Journal of the American Statistical Association, 107(500), 1590–1598. [https://doi.org/10.1080/01621459.2012.737745
](https://www.tandfonline.com/doi/full/10.1080/01621459.2012.737745#:~:text=https%3A//doi.org/10.1080/01621459.2012.737745)



