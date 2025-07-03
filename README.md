# Welcome to rupturesRcpp

## Description

<p align="justify"> The R package provides an efficient, object-oriented R6 interface for offline change point detection, implemented in C++ for high performance, created as part of the Google Summer of Code 2025 program. </p>


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
![image](https://github.com/user-attachments/assets/66e19844-511c-4e4f-a937-26990455402d)

### Segmentation

To perform change point detection, our package requires three "ingredients": **cost function**, **segmentation method**, and **linear penalty threshold**. 

#### Cost function

A `costFun` object can be obtained via `createCostFunc()`. For example, since our example involves regimes with varying variance, a suitable `costFunc` option is `"SIGMA"`. Supported methods include `"L2"`, `"VAR"`, and `"SIGMA"`.

```r
library("rupturesRcpp")
SIGMAObj = createCostFunc(costFunc = "SIGMA")
```
Each cost function may have additional parameters (see `?createCostFunc` for more details). For `"SIGMA"`, `addSmallDiag` and `epsilon` are required. If not specified, the default options will be used.

#### Segmentation method

Our package currently implements two R6 classes for offline change point detection, namely `binSeg` for binary segmentation and `PELT` pruned exact linear time. Their interfaces are similar. Thus, it is sufficient to demonstrate only the usage of `binSeg`.

A `binSeg` object can be initialised as follows:

```r
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFuncObj = SIGMAObj) 
```
Here, `minSize` is the minimum segment length, and `jump` defines a search grid for potential change points. To input a time series matrix and perform binary segmentation for the maximum number of regimes, we can use the `$fit()` method. This is also required for `$predict()`.  

```r
binSegObj$fit(tsMat) 
```
To print the configurations of the `binSeg` object, we can use the `$describe()` method. This method also invisibly returns a list containing several fields of the `binSeg` object. 
```r
binSegObj$describe() 
```
<pre>
Binary Segmentation (binSeg)
minSize      : 1L
jump         : 1L
costFunc.    : "SIGMA"
addSmallDiag : TRUE
epsilon      : 1e-06
fitted       : TRUE
n            : 200L
p            : 2L
</pre>

#### Linear penalty

To obtain an estimated segmentation, we can use the `$predict()` method and specify a non-negative penalty value `pen`. This returns a sorted integer vector of end points, which includes the number of observations by design. The parameter `pen` should be properly tuned. Here, we set `pen = 100`.

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

![image](https://github.com/user-attachments/assets/fe16bf47-fdc5-47cc-9040-456002543b4a)

## Future development

- Provide additional cost functions (e.g., `"Poisson"`, `"Linear-L1"`, and `"Linear-L2"`). 
- Implement other change point detection classes (e.g., `Opt` and `Win`).  
- Develop features for tuning the linear penalty.


