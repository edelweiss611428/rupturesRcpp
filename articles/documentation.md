# Documentation

The documentation covers the three parts of every detection: the cost
function, the segmentation method and model selection. This page lists
the basic functions and shows how the package is organised; the sidebar
links to each part.

## Basic functions

Every segmentation class (`binSeg`, `Window`, `PELT` and `Dynp`) has the
same methods, so code written for one runs with the others.

| Function | What it does |
|----|----|
| `costFunc$new()` | Chooses the cost function, which sets the kind of change to look for. |
| `binSeg$new()`, `Window$new()`, `PELT$new()`, `Dynp$new()` | Chooses the segmentation method, which searches for the change-points. |
| `$fit(tsMat)` | Attaches the time series: a matrix with one row per time point. |
| `$predict(pen = ...)` | Returns the change-points for a penalty per change-point. |
| `$predict(nBkps = ...)` | Returns the change-points for a given number of change-points. |
| `$plot()` | Plots the series with the detected segments. |
| `$segments()` | Gives the cost and parameter estimates of each segment. |
| `$eval(a, b)` | Gives the cost of the segment `(a, b]`. |

## How the package is organised

A detection with K change-points t_1 \< \dots \< t_K, where t_0 = 0 and
t\_{K+1} = n, minimises the total cost of the segments plus a penalty
\lambda per change-point:

\sum\_{k=0}^{K} c\left(y\_{(t_k+1):t\_{k+1}}\right) + \lambda K.

The cost function sets c, the segmentation method searches over the
change-points, and model selection sets \lambda, or fixes K instead. The
package follows the same split.

`rupturesRcpp`

- [Cost
  functions](https://edelweiss611428.github.io/rupturesRcpp/articles/cost-functions.md)

  What kind of change to look for

  - [`costFunc$new()`](https://edelweiss611428.github.io/rupturesRcpp/reference/costFunc.md):
    `"L1"`, `"L2"`, `"SIGMA"`, `"VAR"`, `"LinearL2"`, `"LinearSIGMA"`,
    `"LinearL1"`, `"Custom"`
  - [`costFactory`](https://edelweiss611428.github.io/rupturesRcpp/articles/segment-costs-and-parameters.md):
    costs without a detection

- [Segmentation
  methods](https://edelweiss611428.github.io/rupturesRcpp/articles/segmentation-methods.md)

  How to search for change-points

  - [`binSeg`](https://edelweiss611428.github.io/rupturesRcpp/reference/binSeg.md):
    binary segmentation
  - [`Window`](https://edelweiss611428.github.io/rupturesRcpp/reference/Window.md):
    sliding window
  - [`PELT`](https://edelweiss611428.github.io/rupturesRcpp/reference/PELT.md):
    exact for a penalty
  - [`Dynp`](https://edelweiss611428.github.io/rupturesRcpp/reference/Dynp.md):
    exact for a number of change-points

- [Model
  selection](https://edelweiss611428.github.io/rupturesRcpp/articles/model-selection.md)

  How many change-points to keep

  - `$predict(pen = ...)`: a penalty per change-point
  - `$predict(nBkps = ...)`: a fixed number
  - `$getHistory()` and `$plotElbow()`: the elbow method
