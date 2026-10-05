# Case studies

Each case study is a complete example on its own page, from data to
change-points. Pick the one closest to your problem.

| Case study | Cost function | Method | Also shows |
|----|----|----|----|
| [Change in mean and variance](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-sigma.md) | `"SIGMA"` | `binSeg` | `$describe()`, `$plot()` and the elbow method |
| [Change in autoregressive dynamics](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-var.md) | `"VAR"` | `binSeg` | Changing an object through its active bindings, and the elbow method |
| [Tuning the penalty](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-pen.md) | `"L2"` | `PELT` | Choosing `pen` from a grid of values |
| [Custom cost functions](https://edelweiss611428.github.io/rupturesRcpp/articles/case-study-custom.md) | `"Custom"` | `PELT`, `binSeg` | Data from outside `$fit()` |

Case studies on real data are to be added.
