# Welcome to rupturesRcpp (GSoC 2025-archive)
[![R-CMD-check](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edelweiss611428/rupturesRcpp/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/edelweiss611428/rupturesRcpp/branch/main/graph/badge.svg)](https://codecov.io/gh/edelweiss611428/rupturesRcpp)

## Abstract

<p align="justify"> This project develops an object-oriented R package for offline change-point detection. We implemented three core segmentation algorithms in C++, wrapped them with R6 classes in R, developed five multivariate cost functions, added CI/CD, and tested extensively. The package is installable and nearly ready for CRAN publication, with additional features planned for future work. </p>

## Description

<p align="justify"> This project aims to develop an object-oriented R package for offline change-point detection, based on the widely used Python library `ruptures`, which has received over 500,000 downloads as of April 8, 2025. These methods identify shifts in the underlying structure of time series and are widely used in fields such as signal processing and econometrics. Despite the popularity of `ruptures` and related R packages, their implementations do not achieve optimal runtime complexities. Several efficient C++ implementations exist, but they lack consistent interface and offer limited support for multivariate time series. 
  
To address these limitations, we re-implemented core change-point detection algorithms in `ruptures` in modern C++ and wrapped these implementations in an object-oriented R package with consistent interface. This will make state-of-the-art change-point detection tools more efficient and accessible to practitioners. </p>

<p align="justify"> This branch was created as part of the Google Summer of Code 2025 (GSoC 2025) program to document and archive the work we have done. </p>

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

## Completed work 

In this GSoC 2025 project, we have done the following:

- Implemented three popular offline change-point detection methods—Pruned Exact Linear Time (`PELT`), Binary Segmentation (`binSeg`), and Slicing Window (`Window`)—as template classes in C++ to support various cost functions with minimal overhead. (see the `Getting started/Segmentation methods` section for more details).
- Implemented five C++ classes for multivariate cost functions `"L1"`, `"L2"`, `"VAR"`, `"SIGMA"`, and `"LinearL2"`. Four of these costs can be queried in `O(1)` time given pre-computation (see the `Getting started/Cost functions` section for more details).
- Developed a modern, object-oriented R package interface via R6 classes to wrap these implementations. It is robust and well-engineered for error handling.
- Well-documented the code.
- Tested most of these R/C++ modules for robustness and correctness, achieving an overall coverage of 96% and added CI/CD via GitHub Actions to automatically test and report coverage.
  
Currently, this project can be installed and is close to being ready for practical usage, including potential CRAN publication.

## Outstanding and future work

This project is still ongoing, and there are several areas for improvement and expansion, which are open to future contributors:

1. Testing and validation
   - Increase testing for robustness and correctness of existing modules (e.g., mathematical correctness, time complexities).
2. Documentation and examples
   - Add examples to the package (currently instructions are online only).
   - Provide vignettes demonstrating practical usage, such as tuning/selecting linear penalty thresholds.
   - Provide clearer instructions for future contributors, as the current setup is somewhat complex.
3. Optimisation
   - Improve the `"L1"` cost module, potentially allowing queries in `O(log(n))` time using structures such as a persistent segment tree with `O(nlog(n))` precomputation.
4. Enhancing OOP interface
   - Clean and improve the existing segmentation classes for robustness, including better `$clone()` method (the current implementation may not fully copy all C++ class attributes to a new object), new `$save()` and `$load()` methods to save and reconstruct R6 objects, and an enhanced `$plot()` method for models involving both independent and dependent variables.
5. Code refactoring
   - Further refactor existing C++ and R code for readability and modularity. Each module currently exists in a single file and could be reorganised for easier integration with other programming languages (e.g., separating Rcpp and RcppArmadillo code and reducing Rcpp usage).
6. Feature expansion
   - Add additional cost functions (e.g., `"Poisson"`, `"Linear-L1"`).
   - Implement other offline change-point detection classes (e.g., `Opt`, `BottomUp`).
   - Develop a `costFactory` class for users focused solely on fast cost computation.
   - Develop an R6 class for penalty methods, including potential non-linear penalties.
     
## Contact

For inquiries regarding this GSoC project, the rupturesRcpp package, GSoC in general, or to report issues, please contact: `nguyem39@qut.edu.au`.

---

## Installation

To install the newest version of the package, use the following R code: 

```r
library(devtools)
install_github("edelweiss611428/rupturesRcpp") 
```

## Getting started

To detect change-points using `rupturesRcpp` you need three main components:

- **Cost function** (`costFunc`)
- **Segmentation method** (`binSeg`, `Window`, `PELT`)
- **Linear penalty threshold**

Each `component` is implemented using an R6-based object-oriented design for clarity and flexibility.

### Cost functions

Create a `costFunc` object by specifying the desired (supported) cost function and optional parameters:

```r
library("rupturesRcpp")
costFuncObj = costFunc$new("L2")
costFunc$pass() #output attributes corresponding to the specified cost function.
```
<pre>
$costFunc
[1] "L2"
</pre>

The following table shows the list of supported cost functions. Here, `n` is segment length.

| **Cost function** | **Description**                                                                                  | **Parameters/active bindings**           | **Dimension** | **Time complexity** |
|-------------------|--------------------------------------------------------------------------------------------------|------------------------------------------|----------------|----------------------|
| `"L1"`            | Sum of `L1` distances to the segment-wise median; robust to outliers.                            | `costFunc`                               | `multi`        | O(nlog(n))                 |
| `"L2"`            | Sum of squared `L2` distances to the segment-wise mean; faster but less robust than `L1`.        | `costFunc`                               | `multi`        | O(1)                 |
| `"SIGMA"`         | Log-determinant of empirical covariance; models varying mean&variance.                           | `costFunc`, `addSmallDiag`, `epsilon`    | `multi`        | O(1)                 |
| `"VAR"`           | Sum of squared residuals from a vector autoregressive model with constant noise variance.        | `costFunc`, `pVAR`                       | `multi`        | O(1)                 |
| `"LinearL2"`      | Sum of squared residuals from a linear regression model with constant noise variance.            | `costFunc`, `intercept`                  | `multi`        | O(1)                 |

If active binding `costFunc` is modified by assigning to `costFuncObj$costFunc` and the required parameters are missing, the default parameters will be used.
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

After initialising a `costFunc` object, create a segmentation object such as `binSeg`, `Window`, or `PELT`.

| **R6 Class**     | **Method**                | **Description**                                                                | **Parameters/active bindings**                                |
|------------------|---------------------------|--------------------------------------------------------------------------------|---------------------------------------------------------------|
| `binSeg`         | Binary Segmentation       | Recursively splits the signal at points that minimise the cost.                | `minSize`, `jump`, `costFunc`, `tsMat`, `covariates`          |
| `Window`         | Slicing Window            | Detects change-points using local gains over sliding windows.                  | `minSize`, `jump`, `radius`, `costFunc`,`tsMat`, `covariates` |
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

## Contributing

We welcome all contributions, whether big or small. If you encounter a bug or have a feature request, please open an issue to let us know. 

Feel free to fork the repository and make your changes. For significant updates, it’s best to discuss them with us first. When your changes are ready, submit a pull request.

Thanks for helping us improve this project!

## License

This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) License. See the [LICENSE](LICENSE.md) file for details.


## References

- Hocking, T. D. (2024). *Finite Sample Complexity Analysis of Binary Segmentation*. arXiv preprint arXiv:2410.08654. [https://arxiv.org/abs/2410.08654](https://arxiv.org/abs/2410.08654)

- Truong, C., Oudre, L., & Vayatis, N. (2020). *Selective review of offline change point detection methods*. Signal Processing, 167, 107299. [https://doi.org/10.1016/j.sigpro.2019.107299
](https://www.sciencedirect.com/science/article/abs/pii/S0165168419303494?via%3Dihub#:~:text=https%3A//doi.org/10.1016/j.sigpro.2019.107299)

- Killick, R., Fearnhead, P., & Eckley, I. A. (2012). *Optimal detection of change points with a linear computational cost*. Journal of the American Statistical Association, 107(500), 1590–1598. [https://doi.org/10.1080/01621459.2012.737745
](https://www.tandfonline.com/doi/full/10.1080/01621459.2012.737745#:~:text=https%3A//doi.org/10.1080/01621459.2012.737745)



