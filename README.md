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


### How to cite?

If you use `ruptures` in a scientific publication, we would appreciate citations to the following paper:

- Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299. [[journal]](https://doi.org/10.1016/j.sigpro.2019.107299) [[pdf]](http://www.laurentoudre.fr/publis/TOG-SP-19.pdf)


## Installation

To download the newest version of the package, use the following R code: 

```
library(devtools)
install_github("edelweiss611428/rupturesRcpp") 
```

## Basic usage

### Data-generating process  
  
We consider a simple 2d time series with two piecewise Gaussian regimes and constant spherical covariance.


```r
set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5)),
              c(rnorm(100,0), rnorm(100,5)))
```


### Binary segmentation (binSeg)

To demonstrate the usage of our package, we will use the binSeg class as an example. The PELT class and future classes will follow the same interface. 

First, we initialise a binSeg object with minSize = 1L, jump = 1L, and costFunc = "L2". Currently, "L2" is the only supported cost function.

```r
library("rupturesRcpp")
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2") 
```
To input a time series matrix and perform binSeg with the maximum number of regimes, we can use $fit(). This is also needed for $predict(). 

```r
binSegObj$fit(tsMat) 
```
To print the configurations of the binSeg object, we can use $describe(). This will invisibly return a list containing several fields of the binSeg object.

```r
binSegObj$describe() 
```
<pre>
Binary Segmentation (binSeg)
minSize  : 1L
jump     : 1L
costFunc : "L2"
fitted   : TRUE
n        : 200L
p        : 2L
</pre>

To obtain an estimated segmentation, we can run $predict() and specify a non-negative penalty value for each additional change point. This returns a sorted integer vector of end points, which includes the number of observations by design. Future development will implement methods for tuning the linear penalty.

After running $predict(), a temporary segmentation result is saved to the object, which allows us to plot the segmentation results by dimension without explicitly specifying the segmentation results, although that option is viable. The $plot() method is based on facet_wrap from ggplot2, allowing users to specify the number of columns in the layout. Users can also use the layout operators | and / from patchwork to stack plots.

#### pen = 1 

```r
binSegObj$predict(pen = 1)
```

```r
pen1 = binSegObj$plot(d = 1:2, main = "binSeg: pen = 1", ncol = 2L)
pen1
```

![image](https://github.com/user-attachments/assets/603095ed-71a0-4c22-8dc5-af38c5269e1a)


#### pen = 25
```r
binSegObj$predict(pen = 25) 
pen25 = binSegObj$plot(d = 1:2, main = "binSeg: pen = 25", ncol = 2L)
pen25
```
![image](https://github.com/user-attachments/assets/dea22543-e057-44ae-a174-197d9c704aa5)


#### Vertically stacked
```r

pen1/pen25
```
![image](https://github.com/user-attachments/assets/66c21f4a-bc9d-4096-8be8-464b16340952)

#### User-provided endPts

Here, endPts must include the number of observations. This is the typical output of $predict().

```r
pred25 = binSegObj$predict(pen = 25) 
pen25 = binSegObj$plot(d = 1:2, endPts = pred25, main = "binSeg: pen = 25", ncol = 2L)
pen25
```

![image](https://github.com/user-attachments/assets/1ef897b6-9756-40c7-97fc-8819c6632101)


