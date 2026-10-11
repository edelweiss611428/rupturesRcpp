# `$plot(type = ...)`: fitted values, residuals and residual diagnostics

newDetectors = list(
  PELT   = function(cf) PELT$new(minSize = 5L, costFunc = cf),
  binSeg = function(cf) binSeg$new(minSize = 5L, costFunc = cf),
  Window = function(cf) Window$new(minSize = 5L, radius = 20L, costFunc = cf),
  Dynp   = function(cf) Dynp$new(minSize = 5L, nBkpsMax = 5L, costFunc = cf)
)

# Draws a plot on a null device, so that errors raised at render time are caught too
renderPlot = function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(p)
  invisible(TRUE)
}

# Fitted values of segment (a, b] from `lm()`, for comparison
lmFitted = function(Y, Z, a, b) {
  rows = (a + 1):b
  sapply(seq_len(ncol(Y)), function(j) fitted(lm(Y[rows, j] ~ Z[rows, , drop = FALSE] - 1)))
}

set.seed(1)
xReg = cbind(rnorm(200))
yReg = cbind(c(2 * xReg[1:100], -2 * xReg[101:200]) + rnorm(200, sd = 0.5))

test_that("`type` is validated", {

  PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("LinearL2"))
  PELTObj$fit(yReg, xReg)
  PELTObj$predict(pen = 30)

  expect_error(PELTObj$plot(type = "a"), "must be one of")
  expect_error(PELTObj$plot(type = c("fitted", "residuals")), "must be one of")
  expect_error(PELTObj$plot(type = NA), "must be one of")
  expect_error(PELTObj$plot(type = NA_character_), "must be one of")
  expect_error(PELTObj$plot(type = NULL), "must be one of")
  expect_error(PELTObj$plot(type = 1), "must be one of")
  expect_error(PELTObj$plot(type = "Fitted"), "must be one of")

})

test_that("Every type works for every segmentation class", {

  for (nm in names(newDetectors)) {

    obj = newDetectors[[nm]](costFunc$new("LinearL2"))
    obj$fit(yReg, xReg)
    expect_equal(obj$predict(pen = 30), c(100, 200), info = nm)

    pSeries = suppressMessages(obj$plot())
    expect_s3_class(pSeries, "ggplot")
    expect_null(pSeries$labels$subtitle)

    for (type in c("series", "fitted", "residuals")) {
      p = suppressMessages(obj$plot(type = type))
      expect_s3_class(p, "ggplot")
      expect_true(renderPlot(p))
    }

    pDiag = suppressMessages(obj$plot(type = "diagnostics"))
    expect_s3_class(pDiag, "patchwork")
    expect_s3_class(pDiag, "ggplot")
    expect_true(renderPlot(pDiag))

  }

})

test_that("The fitted values drawn are those of the segment parameters", {

  PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("LinearL2"))
  PELTObj$fit(yReg, xReg)
  PELTObj$predict(pen = 30)

  p = suppressMessages(PELTObj$plot(type = "fitted"))
  fitLayer = p$layers[[3]]$data
  expected = rbind(lmFitted(yReg, cbind(1, xReg), 0, 100), lmFitted(yReg, cbind(1, xReg), 100, 200))

  expect_equal(fitLayer$time, 1:200)
  expect_equal(fitLayer$fitted, c(expected), ignore_attr = TRUE)
  expect_equal(p$labels$subtitle, "Fitted values (LinearL2)")

  pRes = suppressMessages(PELTObj$plot(type = "residuals"))
  expect_equal(pRes$data$residual, c(yReg - expected), ignore_attr = TRUE)
  expect_equal(pRes$labels$subtitle, "Residuals (LinearL2)")

  pDiag = suppressMessages(PELTObj$plot(type = "diagnostics"))
  expect_equal(pDiag$patches$annotation$subtitle, "Residual diagnostics (LinearL2)")

})

test_that("Fitted values of the regression costs match `lm()` on each segment", {

  set.seed(2)
  n = 120
  X = cbind(rnorm(n), runif(n))
  Y = cbind(X %*% c(1, -2) + rnorm(n), X %*% c(-1, 3) + rnorm(n, sd = 2))
  endPts = c(30L, 75L, 120L)
  starts = c(0L, 30L, 75L)

  for (cf in c("LinearL2", "LinearSIGMA")) {

    for (intercept in c(TRUE, FALSE)) {

      facObj = costFactory$new(costFunc$new(cf, intercept = intercept))
      facObj$fit(Y, X)
      fitted = .segmentFitted(costFunc$new(cf, intercept = intercept)$pass(), Y, X, endPts,
                              function(a, b) facObj$get_params(a, b))

      Z = if (intercept) cbind(1, X) else X
      expected = do.call(rbind, lapply(1:3, function(k) lmFitted(Y, Z, starts[k], endPts[k])))
      expect_equal(fitted, expected, ignore_attr = TRUE, info = paste(cf, intercept))

    }
  }

  # LinearL1: X %*% coef, with `coef` from `$segments()`
  PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("LinearL1"))
  PELTObj$fit(Y, X)
  PELTObj$predict(pen = 20)
  p = suppressMessages(PELTObj$plot(d = 1:2, type = "fitted"))

  expected = do.call(rbind, lapply(PELTObj$segments(), function(s) {
    cbind(1, X[(s$Start + 1):s$End, , drop = FALSE]) %*% s$Params$coef
  }))
  expect_equal(p$layers[[3]]$data$fitted, c(expected))

})

test_that("Regression costs fitted without covariates use the intercept-only model", {

  set.seed(3)
  Y = cbind(c(rnorm(50, 0), rnorm(50, 4)))

  PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("LinearL2"))
  expect_warning(PELTObj$fit(Y), "No `covariates` found")
  PELTObj$predict(pen = 20)

  p = suppressMessages(PELTObj$plot(type = "fitted"))
  segMeans = unlist(lapply(PELTObj$segments(), function(s) rep(mean(Y[(s$Start + 1):s$End]), s$End - s$Start)))
  expect_equal(p$layers[[3]]$data$fitted, segMeans)

})

test_that("Fitted values of `\"VAR\"` are the one-step predictions from the lagged series", {

  set.seed(4)
  n = 150
  Y = matrix(rnorm(2 * n), ncol = 2)
  for (t in 3:n) {
    Y[t, ] = Y[t, ] + 0.5 * Y[t - 1, ] - 0.2 * Y[t - 2, ]
  }
  pVAR = 2L
  endPts = c(1L, 60L, 150L)

  facObj = costFactory$new(costFunc$new("VAR", pVAR = pVAR))
  facObj$fit(Y)
  fitted = .segmentFitted(costFunc$new("VAR", pVAR = pVAR)$pass(), Y, NULL, endPts,
                          function(a, b) facObj$get_params(a, b))

  # The first `pVAR` observations have no lags; segment (0, 1] has no fit at all
  expect_true(all(is.na(fitted[1:pVAR, ])))

  Z = cbind(1, rbind(NA, Y[-n, ]), rbind(NA, NA, Y[-c(n - 1, n), ]))
  for (k in 2:3) {
    a = c(0L, endPts)[k]
    b = endPts[k]
    rows = (max(a, pVAR) + 1):b
    expect_equal(fitted[rows, ], lmFitted(Y, Z, max(a, pVAR), b), ignore_attr = TRUE)
  }

  # The plot leaves out the segment without a fit and says so
  binSegObj = binSeg$new(minSize = 1L, costFunc = costFunc$new("VAR", pVAR = pVAR))
  suppressWarnings(binSegObj$fit(Y))
  expect_message(p <- binSegObj$plot(endPts = endPts, dimNames = c("a", "b"), d = 1:2, type = "fitted"),
                 "Segment\\(s\\) 1 have too few observations")
  expect_equal(p$labels$subtitle, "Fitted values (VAR(2))")
  expect_equal(nrow(p$layers[[3]]$data), 2 * (n - pVAR))
  expect_true(renderPlot(p))

})

test_that("Fitted values of `\"L1\"`, `\"L2\"` and `\"SIGMA\"` are segment medians and means", {

  set.seed(5)
  Y = cbind(c(rnorm(40), rnorm(60, 3)), rexp(100))
  endPts = c(40L, 100L)

  for (cf in c("L1", "L2", "SIGMA")) {

    facObj = costFactory$new(costFunc$new(cf))
    facObj$fit(Y)
    fitted = .segmentFitted(costFunc$new(cf)$pass(), Y, NULL, endPts,
                            function(a, b) facObj$get_params(a, b))

    centre = if (cf == "L1") function(v) apply(v, 2, median) else colMeans
    expected = rbind(matrix(centre(Y[1:40, ]), 40, 2, byrow = TRUE),
                     matrix(centre(Y[41:100, ]), 60, 2, byrow = TRUE))
    expect_equal(fitted, expected, info = cf)

  }

})

test_that("Error handling for the fitted, residual and diagnostic plots", {

  # Not available for custom cost functions
  customCost = costFunc$new("Custom", evalFun = function(segment, start, end) sum((segment - mean(segment))^2))
  PELTObj = PELT$new(minSize = 5L, costFunc = customCost)
  PELTObj$fit(yReg)
  PELTObj$predict(pen = 30)
  expect_no_error(suppressMessages(PELTObj$plot()))
  for (type in c("fitted", "residuals", "diagnostics")) {
    expect_error(suppressMessages(PELTObj$plot(type = type)), "not \"Custom\"")
  }

  # `$fit()` must have been run
  PELTObj = PELT$new(minSize = 5L, costFunc = costFunc$new("L2"))
  PELTObj$tsMat = yReg
  expect_no_error(suppressMessages(PELTObj$plot(endPts = c(100, 200))))
  expect_error(suppressMessages(PELTObj$plot(endPts = c(100, 200), type = "fitted")), "must be run before")

  # No segment long enough to estimate the regression
  X3 = matrix(rnorm(600), ncol = 3)
  PELTObj = PELT$new(minSize = 1L, costFunc = costFunc$new("LinearL2"))
  PELTObj$fit(yReg, X3)
  expect_error(suppressMessages(PELTObj$plot(endPts = seq(2, 200, by = 2), type = "residuals")),
               "No fitted values to plot")

  # `covariates` that do not match the fitted coefficients
  expect_error(.segmentFitted(list(costFunc = "LinearL2", intercept = TRUE), yReg, cbind(xReg, xReg), c(100L, 200L),
                              function(a, b) list(coef = matrix(0, 2, 1))),
               "do not match the fitted model")
  expect_error(.segmentFitted(list(costFunc = "LinearL2", intercept = TRUE), yReg, xReg[1:50, , drop = FALSE],
                              c(100L, 200L), function(a, b) list(coef = matrix(0, 2, 1))),
               "do not match")

})
