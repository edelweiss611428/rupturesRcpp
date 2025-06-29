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

First, we initialise a binSeg object with minSize = 1L, jump = 1L, and costFunc = "L2". Currently, "L2", "SIGMA", and "VAR" are the only supported cost functions.

```r
library("rupturesRcpp")
binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2") 
```
To input a time series matrix and perform binSeg with the maximum number of regimes, we can use $fit(). This is also needed for $predict(). 

```r
binSegObj$fit(tsMat) 
```
To print the configurations of the binSeg object, we can use $describe(). This method invisibly returns a list containing several fields of the binSeg object. Since all attributes are stored in private fields to prevent direct modification, this feature may be useful for extracting relevant information.

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
![image](https://github.com/user-attachments/assets/28d6bc36-3b6c-45a0-9dee-bdd36b44e5bd)

#### User-provided endPts

Here, endPts must include the number of observations. This is the typical output of $predict().

```r
pred25 = binSegObj$predict(pen = 25) 
pen25 = binSegObj$plot(d = 1:2, endPts = pred25, main = "binSeg: pen = 25", ncol = 2L)
pen25
```

![image](https://github.com/user-attachments/assets/1ef897b6-9756-40c7-97fc-8819c6632101)

## Future development

- Provide additional cost functions (e.g., L1, Poisson).  
- Implement other change point detection methods (e.g., Opt and Win).  
- Develop methods for tuning the linear penalty.


