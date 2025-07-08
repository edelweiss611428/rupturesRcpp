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

### Simulated data example
  
To demonstrate the package usage, we consider a simple 2d time series with two piecewise Gaussian regimes and varying variance.

```r
set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)),
              c(rnorm(100,0), rnorm(100,5,5)))
```
![image](https://github.com/user-attachments/assets/73687865-b52e-4a6a-b8fd-5cf700a9be7a)
### Segmentation

To perform change point detection, our package requires three "ingredients": **cost function**, **segmentation method**, and **linear penalty threshold**. 

#### Cost function

An `R6` object of class `costFunc` can be obtained via `costFunc$new()`. As our example involves regimes with varying variance, a suitable `costFunc` option is `"SIGMA"`. Supported methods include `"L2"`, `"VAR"`, and `"SIGMA"`.

```r
library("rupturesRcpp")
SIGMAObj = costFunc$new("SIGMA")
```
Each cost function may have some optional parameters (see `?costFunc` for more details). For `"SIGMA"`, we need to specify `addSmallDiag` and `epsilon`. Here, if `addSmallDiag = TRUE`, a small `epsilon` is added to the diagonal of estimated covariance matrices, which 
stabilises matrix operations. If any are not specified, the default options will be used.

#### Segmentation method

Our package currently implements two R6 classes for offline change point detection, namely `binSeg` for binary segmentation and `PELT` for pruned exact linear time. Their interfaces are similar. Thus, it is sufficient to demonstrate the usage of only `binSeg`.

A `binSeg` object can be initialised as follows:

```r
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = SIGMAObj) 
```
Here, `minSize` is the minimum segment length, and `jump` defines a search grid for potential change points. To construct a `C++` module for `binSeg`, we can use the `$fit()` method. Once `binSeg` is fitted, we can use the methods `$predict()` and `$eval()`.  

```r
binSegObj$fit(tsMat) 
```
To print the configurations of the `binSeg` object, we can use the `$describe()` method. This method also invisibly returns a list containing several fields of the `binSeg` object for extraction purpose.
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

An implicit way to modify the fields of a binSeg object is through its active bindings, which can be used not only to extract key fields but also to update them.

The `R6` class `binSeg` has 4 active bindings, namely `minSize`, `jump`, `costFunc` and `tsMat`. We can modify an existing `binSeg` object by assigning new values to its active bindings. To demonstrate this, we
consider a piecewise vector autoregressive example with constant noise variance.

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
Modifying `tsMat` will automatically trigger `$fit()`. 

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
- Implement other change point detection classes (e.g., `Opt` and `Win`).  

## References

- Hocking, T. D. (2024). *Finite Sample Complexity Analysis of Binary Segmentation*. arXiv preprint arXiv:2410.08654. [https://arxiv.org/abs/2410.08654](https://arxiv.org/abs/2410.08654)

- Truong, C., Oudre, L., & Vayatis, N. (2020). *Selective review of offline change point detection methods*. Signal Processing, 167, 107299. [https://doi.org/10.1016/j.sigpro.2019.107299
](https://www.sciencedirect.com/science/article/abs/pii/S0165168419303494?via%3Dihub#:~:text=https%3A//doi.org/10.1016/j.sigpro.2019.107299)

- Killick, R., Fearnhead, P., & Eckley, I. A. (2012). *Optimal detection of change points with a linear computational cost*. Journal of the American Statistical Association, 107(500), 1590–1598. [https://doi.org/10.1080/01621459.2012.737745
](https://www.tandfonline.com/doi/full/10.1080/01621459.2012.737745#:~:text=https%3A//doi.org/10.1080/01621459.2012.737745)



