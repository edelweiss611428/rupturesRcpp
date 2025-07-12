# Welcome to rupturesRcpp

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

## Basic usage

To perform change point detection, `rupturesRcpp` requires three "ingredients": **cost function**, **segmentation method**, and **linear penalty threshold**. 

### `costFunc` objects

A cost object can be created by initialising a `R6` object of class `costFunc`, which specifies a supported cost function  (e.g., `"L2"`)  and its parameters. 

```r
library("rupturesRcpp")
costFuncObj = costFunc$new(costFunc = "L2")
```
The following table shows the list of supported cost functions.

| **Cost functions**| **Description**                                                                                  | **Parameters/active bindings**           | **Dimension**       |
|-------------------|--------------------------------------------------------------------------------------------------|------------------------------------------|---------------------|
| `"L1"`            | Sum of L1 distances to the segment-wise median; robust to outliers.                              | `costFunc`                               | `multi`             |
| `"L2"`            | Sum of squared L2 distances to the segment-wise mean; faster but less robust than L1.            | `costFunc`                               | `multi`             |
| `"SIGMA"`         | Log-determinant of empirical covariance; models varying variance.                                | `costFunc`, `addSmallDiag`, `epsilon`    | `multi`             |
| `"VAR"`           | Residual error from vector autoregression with constant noise.                                   | `costFunc`, `pVAR`                       | `multi`             |

By default, `costFunc$new()` creates a `"L2"` cost object. If cost-specific parameters are not specified, the default options will be used. See `?rupturesRcpp::costFunc` for more details.

### Segmentation methods

After initialising a cost object, we can specify a change-point detection object. Supported methods include `binSeg` for binary segmentation, `PELT` for pruned exact linear time, and `Window`
for window slicing. These have been implemented in `R6` classes (see the table below for more details). 

| **R6 Class**     | **Method**                | **Description**                                                                | **Parameters/active bindings**                 |
|------------------|---------------------------|--------------------------------------------------------------------------------|------------------------------------------------|
| `binSeg`         | Binary Segmentation       | Recursively splits the signal at points that minimise the cost.                | `minSize`, `jump`, `costFunc`, `tsMat`         |
| `Window`         | Window Slicing            | Detects change-points using local gains over sliding windows.                  | `minSize`, `jump`, `radius`, `costFunc`,`tsMat`|
| `PELT`           | Pruned Exact Linear Time  | Optimal segmentation with pruning for linear-time performance.                 | `minSize`, `jump`, `costFunc`, `tsMat`         |

A `PELT` object, for example, can be initialised as follows:
```r
detectionObj = PELT$new(minSize = 1L, jump = 1L, costFunc = costFuncObj)
```
If method-specific parameters are not provided, default values will be used—except for `tsMat`, which is not required at initialisation and can be set later.

The `R6` detection classes share a consistent object-oriented interface with similar methods. For example, the `PELT` class provides the following methods:

- `$new()`: Initialises a `PELT` object.
- `$describe()`: Describes the `PELT` object.
- `$fit()`: Constructs a `C++` `PELT` module corresponding to the specified parameters/active binding.
- `$predict()`: Performs change-point detection given a linear penalty value.
- `$eval()`: Evaluate the cost of a segment.
- `$plot()`: Plots change-point segmentation in `ggplot` style.
- `$clone()`: Clones the `PELT` object.

The `$fit()` method must be called before using `$predict()` or `$eval()`, as it initialises the underlying `C++` module for efficient change-point detection.

Class-specific parameters/active bindings can be modified after initialisation via assignment operator. For example, to modify the `minSize` field in the `PELT` object
to `25L`, we can simply use the `$` operator:

```r
detectionObj$minSize = 25L.
```

This will modifies the value of`private$.minSize` to `25L`. We can also use the `$` operator - `detectionObj$minSize` - to extract `minSize`.

Whenever an active binding is set or modified, internal diagnostics or re-fitting may be triggered automatically to ensure consistency. For example,
if a `C++` object has been created for `minSize = 1L`, modifying `minSize` will automatically trigger `self$fit()`.


### Simulated data example
  
To demonstrate the package usage, we consider a simple 2d time series with two piecewise Gaussian regimes and varying variance.

```r
set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)),
              c(rnorm(100,0), rnorm(100,5,5)))
```
![image](https://github.com/user-attachments/assets/73687865-b52e-4a6a-b8fd-5cf700a9be7a)

### Segmentation

#### Creating a `costFunc` object

As our example involves regimes with varying variance, a suitable `costFunc` option is `"SIGMA"`. 

```r
SIGMAObj = costFunc$new("SIGMA", addSmallDiag = TRUE, epsilon = 1e-6)
```
For `"SIGMA"`, we need to specify `addSmallDiag` and `epsilon`. If `addSmallDiag = TRUE`, a small `epsilon` is added to the diagonal of estimated covariance matrices, which stabilises matrix operations. 

#### Initialising a `binSeg` object

As the detection classes' interfaces are similar, it is sufficient to demonstrate the usage of `binSeg` only.

A `binSeg` object can be initialised as follows:

```r
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = SIGMAObj) 
```
Then, we construct a `C++` module for `binSeg` via the `$fit()` method, which requires a `tsMat`. Once `binSeg` is fitted, we will be able to use `$predict()` and `$eval()`.  

```r
binSegObj$fit(tsMat) 
```
To view the configurations of the `binSeg` object, we can use the `$describe()` method. This method also invisibly returns a list containing several fields of the `binSeg` object for further extraction.
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

#### Linear penalty

To obtain an estimated segmentation, we can use the `$predict()` method and specify a non-negative penalty value `pen`. This returns a sorted integer vector of end-points, which includes the number of observations by design. The parameter `pen` should be properly tuned. Here, we set `pen = 100`.

```r
binSegObj$predict(pen = 100)
```
<pre>
[1] 100 200
</pre>

After running `$predict()`, the segmentation output is temporarily saved to the `binSeg` object, allowing users to plot the segmentation results via the `$plot()` method, based on `ggplot2::facet_wrap`. Users can also use the layout operators `|` and `/` from `patchwork` to stack plots.

```r
binSegObj$plot(d = 1:2, 
               main = "method: binSeg; costFunc: SIGMA; pen: 100")
```
![image](https://github.com/user-attachments/assets/d5d31c3d-ced1-4667-8de5-e9ad0cdc84ec)

#### Active bindings

An implicit way to modify or set the fields of a `binSeg` object is through its active bindings. We can modify an existing `binSeg` object by assigning new values to its active bindings instead of creating a new object. 
To demonstrate this, we consider a piecewise vector autoregressive example with constant noise variance.

```r
set.seed(1)
tsMat = matrix(c(filter(rnorm(100), filter = 0.9, method = "recursive"), 
                 filter(rnorm(100), filter = -0.9, method = "recursive")))
```
![image](https://github.com/user-attachments/assets/de7a44ef-bc21-4348-b77d-3c05f1ff9c0a)

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
![image](https://github.com/user-attachments/assets/37e45eaa-43f0-492b-ba2f-2b3a7ab6c72a)


## Future development

- Provide additional cost functions (e.g., `"Poisson"`, `"Linear-L1"`, and `"Linear-L2"`). 
- Implement other change-point detection classes (e.g., `Opt` and `Dynp`).  

## References

- Hocking, T. D. (2024). *Finite Sample Complexity Analysis of Binary Segmentation*. arXiv preprint arXiv:2410.08654. [https://arxiv.org/abs/2410.08654](https://arxiv.org/abs/2410.08654)

- Truong, C., Oudre, L., & Vayatis, N. (2020). *Selective review of offline change point detection methods*. Signal Processing, 167, 107299. [https://doi.org/10.1016/j.sigpro.2019.107299
](https://www.sciencedirect.com/science/article/abs/pii/S0165168419303494?via%3Dihub#:~:text=https%3A//doi.org/10.1016/j.sigpro.2019.107299)

- Killick, R., Fearnhead, P., & Eckley, I. A. (2012). *Optimal detection of change points with a linear computational cost*. Journal of the American Statistical Association, 107(500), 1590–1598. [https://doi.org/10.1080/01621459.2012.737745
](https://www.tandfonline.com/doi/full/10.1080/01621459.2012.737745#:~:text=https%3A//doi.org/10.1080/01621459.2012.737745)



