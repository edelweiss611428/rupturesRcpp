#Dynp module
set.seed(12345)

#Total cost of a segmentation `bkps` under a fitted segmentation object (mirrors PELT's `segObj`)
segCost = function(obj, bkps) {
  starts = c(0, head(bkps, -1))
  sum(mapply(function(s, e) obj$eval(s, e), starts, bkps))
}


test_that("Dynp_L1 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L1")
  DynpObj = Dynp$new(costFunc = costFuncObj)
  DynpObj$fit(tsMat)

  expect_equal( DynpObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Dynp_L2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L2")
  DynpObj = Dynp$new(costFunc = costFuncObj)
  DynpObj$fit(tsMat)

  expect_equal( DynpObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Dynp_SIGMA works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("SIGMA")
  DynpObj = Dynp$new(costFunc = costFuncObj)
  DynpObj$fit(tsMat)

  expect_equal( DynpObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Dynp_VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("VAR")
  DynpObj = Dynp$new(costFunc = costFuncObj)

  expect_warning(DynpObj$fit(tsMat), "Some systems seem singular!") #Warning once feature
  expect_equal(DynpObj$predict(pen = 0.1), seq(50,150,50))
  expect_warning(expect_false(DynpObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(DynpObj$eval(0,50), 0)), "seems singular")

})


test_that("Dynp_LinearL2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("LinearL2")
  DynpObj = Dynp$new(costFunc = costFuncObj)

  #when no covariate matrix is provided
  expect_warning(DynpObj$fit(tsMat), "an intercept")
  expect_equal(DynpObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(DynpObj$eval(0,51) == 0)
  expect_true(all.equal(DynpObj$eval(0,50), 0))
})


test_that("Active binding `minSize` works as intended", {

  set.seed(12345)
  #Shift not at the midpoint, so minSize = 1's unconstrained optimum differs from the exact
  #midpoint that minSize = 50 structurally forces
  tsMat = matrix(c(rnorm(23,0), rnorm(77,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(Dynp$new(costFunc = costFuncObj, minSize = "a"))
  expect_error(Dynp$new(costFunc = costFuncObj, minSize = NULL))
  expect_error(Dynp$new(costFunc = costFuncObj, minSize = 1:2))
  expect_error(Dynp$new(costFunc = costFuncObj, minSize = 0))

  #Getter
  DynpObj = Dynp$new(costFunc = costFuncObj)
  expect_equal(DynpObj$minSize, 1L)

  #Setter
  DynpObj$minSize = 5L
  expect_equal(DynpObj$minSize, 5L)
  expect_error(DynpObj$minSize <- "a")
  expect_error(DynpObj$minSize <- NULL)
  expect_error(DynpObj$minSize <- 1:2)
  expect_error(DynpObj$minSize <- 0)

  #Modifying `minSize` triggers refitting if fitted. With `nBkps = 1`, the single change-point
  #must respect `minSize` on both sides, so minSize = 50 on n = 100 forces it to fall exactly at 50
  DynpObj$minSize = 1L
  DynpObj$fit(tsMat)
  ms1Seg = DynpObj$predict(nBkps = 1)

  DynpObj$minSize = 50L
  ms50Seg = DynpObj$predict(nBkps = 1)
  expect_equal(ms50Seg, c(50, 100))
  expect_false(identical(ms1Seg, ms50Seg))

})

test_that("Active binding `jump` works as intended", {

  set.seed(12345)
  #Shift not on a multiple of 5, so jump = 5's grid-snapped optimum differs from jump = 1's
  tsMat = matrix(c(rnorm(23,0), rnorm(77,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(Dynp$new(costFunc = costFuncObj, jump = "a"))
  expect_error(Dynp$new(costFunc = costFuncObj, jump  = NULL))
  expect_error(Dynp$new(costFunc = costFuncObj, jump  = 1:2))
  expect_error(Dynp$new(costFunc = costFuncObj, jump  = 0))
  expect_no_error(Dynp$new())

  #Getter
  DynpObj = Dynp$new(costFunc = costFuncObj)
  expect_equal(DynpObj$jump, 1L)

  #Setter
  DynpObj$jump  = 5L
  expect_equal(DynpObj$jump, 5L)
  expect_error(DynpObj$jump <- "a")
  expect_error(DynpObj$jump <- NULL)
  expect_error(DynpObj$jump <- 1:2)
  expect_error(DynpObj$jump <- 0)

  #Modifying `jump` triggers refitting if fitted; the single change-point must land on the
  #admissible {jump, 2*jump, ...} grid
  DynpObj$jump  = 1L
  DynpObj$fit(tsMat)
  j1Seg = DynpObj$predict(nBkps = 1)

  DynpObj$jump  = 5L
  j5Seg = DynpObj$predict(nBkps = 1)
  expect_equal(j5Seg[1] %% 5, 0)
  expect_false(identical(j1Seg, j5Seg))

})

test_that("Active binding `nBkpsMax` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(15,0), rnorm(15,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: must be NULL or a single non-negative integer
  expect_error(Dynp$new(costFunc = costFuncObj, nBkpsMax = "a"))
  expect_error(Dynp$new(costFunc = costFuncObj, nBkpsMax = 1:2))
  expect_error(Dynp$new(costFunc = costFuncObj, nBkpsMax = -1))
  expect_no_error(Dynp$new(costFunc = costFuncObj, nBkpsMax = NULL))
  expect_no_error(Dynp$new(costFunc = costFuncObj, nBkpsMax = 0))

  #Getter: default is NULL (auto-resolved at `$fit()`)
  DynpObj = Dynp$new(costFunc = costFuncObj)
  expect_null(DynpObj$nBkpsMax)

  #Setter
  DynpObj$nBkpsMax = 5L
  expect_equal(DynpObj$nBkpsMax, 5L)
  expect_error(DynpObj$nBkpsMax <- "a")
  expect_error(DynpObj$nBkpsMax <- 1:2)
  expect_error(DynpObj$nBkpsMax <- -1)
  expect_no_error(DynpObj$nBkpsMax <- NULL)

  #Modifying `nBkpsMax` triggers refitting if fitted
  DynpObj$nBkpsMax = 3L
  DynpObj$fit(tsMat)
  expect_equal(DynpObj$describe()$resolvedNBkpsMax, 3L)
  expect_length(DynpObj$costPath(), 4L) #k = 0,...,3

  DynpObj$nBkpsMax = 5L
  expect_equal(DynpObj$describe()$resolvedNBkpsMax, 5L)
  expect_length(DynpObj$costPath(), 6L) #k = 0,...,5

  #`nBkpsMax = 0` -> only k = 0 is available; requesting more is capped (message), not an error
  DynpObj$nBkpsMax = 0L
  expect_length(DynpObj$costPath(), 1L)
  expect_equal(DynpObj$predict(nBkps = 0), 30L)
  expect_equal(DynpObj$predict(pen = 0), 30L)
  expect_message(bkpsCapped <- DynpObj$predict(nBkps = 1), "exceeds `nBkpsMax`")
  expect_equal(bkpsCapped, 30L)

})

test_that("`nBkpsMax` resolution: auto (message) when NULL, capped (warning) when infeasible", {

  set.seed(1)
  tsMat = matrix(rnorm(60))

  #Auto-resolve to min(20, floor(n/minSize) - 1)
  DynpObj = Dynp$new(costFunc = costFunc$new("L2"), minSize = 1L)
  expect_message(DynpObj$fit(tsMat), "`nBkpsMax` not set")
  expect_equal(DynpObj$describe()$resolvedNBkpsMax, 20L) #min(20, 59)
  expect_null(DynpObj$nBkpsMax) #user setting is left untouched

  #Auto-resolve when the natural maximum is below the 20 cap
  DynpObj2 = Dynp$new(costFunc = costFunc$new("L2"), minSize = 25L)
  expect_message(DynpObj2$fit(tsMat), "`nBkpsMax` not set")
  expect_equal(DynpObj2$describe()$resolvedNBkpsMax, 1L) #floor(60/25) - 1 = 1

  #Explicit `nBkpsMax` above the natural maximum is capped, with a warning
  DynpObj3 = Dynp$new(costFunc = costFunc$new("L2"), minSize = 20L, nBkpsMax = 10L)
  expect_warning(DynpObj3$fit(tsMat), "exceeds the maximum feasible value")
  expect_equal(DynpObj3$describe()$resolvedNBkpsMax, 2L) #floor(60/20) - 1 = 2

  #Explicit `nBkpsMax` within range is used as-is, with neither message nor warning
  DynpObj4 = Dynp$new(costFunc = costFunc$new("L2"), minSize = 1L, nBkpsMax = 5L)
  expect_no_message(DynpObj4$fit(tsMat))
  expect_no_warning(DynpObj4$fit(tsMat))
  expect_equal(DynpObj4$describe()$resolvedNBkpsMax, 5L)

})

test_that("Active binding `tsMat` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  tsMat2 = matrix(c(rnorm(50,0), rnorm(50,5)))
  tsNA = tsMat
  tsNA[1] = NA
  DynpObj = Dynp$new()

  #Wrong input: `tsMat` must be a numeric matrix without any NA
  expect_error(DynpObj$fit(tsNA))
  expect_error(DynpObj$fit(c("a", "b")))
  expect_error(DynpObj$fit(NULL))
  expect_error(DynpObj$fit(as.vector(tsMat)))
  expect_no_error(DynpObj$fit(tsMat))

  #Getter
  DynpObj$fit(tsMat)
  expect_equal(DynpObj$tsMat, tsMat)

  #Setter
  expect_error(DynpObj$tsMat <- tsNA)
  expect_error(DynpObj$tsMat <- c("a", "b"))
  expect_error(DynpObj$tsMat <- NULL)
  expect_error(DynpObj$tsMat <- as.vector(tsMat))

  DynpObj$tsMat = tsMat2
  expect_equal(DynpObj$tsMat, tsMat2)

  #Modify `tsMat` triggers refitting if fitted

  DynpObj$fit(tsMat)
  tsMat1Eval = DynpObj$eval(0,100)
  expect_equal(tsMat1Eval, sum((tsMat - mean(tsMat))^2))

  DynpObj$fit(tsMat2)
  tsMat2Eval = DynpObj$eval(0,100)
  expect_equal(tsMat2Eval, sum((tsMat2 - mean(tsMat2))^2))

})


test_that("Active binding `covariates` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  covariateMat = as.matrix(tsMat/2 + rnorm(100))
  covariateNA = covariateMat
  covariateNA[1] = NA
  costFuncObj = costFunc$new("LinearL2")
  DynpObj = Dynp$new(costFunc = costFuncObj)

  #Wrong input: `covariate` must be a numeric matrix without any NA
  expect_error(DynpObj$fit(tsMat, as.matrix(covariateMat[-1]))) #nrows do not match
  expect_error(DynpObj$fit(tsMat, covariateNA))
  expect_error(DynpObj$fit(tsMat, c("a", "b")))
  expect_error(DynpObj$fit(tsMat, as.vector(covariateMat)))
  expect_no_error(DynpObj$fit(tsMat, covariateMat))
  expect_no_error(DynpObj$fit(tsMat, NULL))

  #Getter
  DynpObj$fit(tsMat, covariateMat)
  expect_equal(DynpObj$covariates, covariateMat)

  #Setter
  expect_error(DynpObj$covariates <- as.matrix(covariateMat[-1])) #nrows do not match
  expect_error(DynpObj$covariates <- covariateNA)
  expect_error(DynpObj$covariates <- c("a", "b"))
  expect_error(DynpObj$covariates <- NULL)
  expect_error(DynpObj$covariates <- as.vector(covariateMat))

  #Modify `tsMat` triggers refitting if fitted

  set.seed(100)
  covariateMat2 = as.matrix(tsMat/2 + rnorm(100))

  DynpObj$fit(tsMat, covariateMat)
  cM1err = DynpObj$eval(0,100)
  DynpObj$covariates = covariateMat2
  cM2err = DynpObj$eval(0,100)

  expect_equal(cM1err, sum(lm(tsMat~covariateMat)$residuals^2))
  expect_equal(cM2err, sum(lm(tsMat~covariateMat2)$residuals^2))

})

test_that("Active binding `costFunc` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(Dynp$new(costFunc = 1L))
  expect_error(Dynp$new(costFunc = list(costFunc = "L2")))
  expect_error(Dynp$new(costFunc = NULL))
  expect_no_error(Dynp$new())

  #Getter
  DynpObj = Dynp$new(costFunc = costFuncObj)
  expect_equal(DynpObj$costFunc$pass()$costFunc, "L2")

  #Setter
  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(DynpObj$costFunc <- list(costFunc = "L2"))
  expect_error(DynpObj$costFunc <- NULL)
  expect_error(DynpObj$costFunc <- 1L)
  expect_no_error(DynpObj$costFunc <- costFunc$new("VAR"))

  #Modifying `costFunc` triggers refitting if fitted

  DynpObj = Dynp$new() #L2
  DynpObj$fit(tsMat)
  expect_equal(DynpObj$eval(0, 100), sum((tsMat - mean(tsMat))^2))

  DynpObj$costFunc = costFunc$new("L1")
  expect_equal(DynpObj$eval(0, 100), sum(abs(tsMat - median(tsMat))))

  DynpObj$costFunc = costFunc$new("SIGMA")
  expect_equal(DynpObj$eval(0, 100), 100*log(det(var(tsMat)*99/100+10^-6)))

})


test_that("Test that `describe()` method works properly)", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")
  DynpObj = Dynp$new(costFunc = costFuncObj)
  expect_message(DynpObj$fit(tsMat), "`nBkpsMax` not set")

  #Wrong input
  expect_error(DynpObj$describe(NULL))
  expect_error(DynpObj$describe(1:2))
  expect_error(DynpObj$describe("a"))

  #`describe()` returns the expected outputs
  expect_no_error(DynpObj$describe(T))
  expect_no_error(DynpObj$describe(F))
  expect_no_error(DynpObj$describe())

  expect_equal(DynpObj$describe(T), DynpObj$describe(F))
  expect_equal(DynpObj$describe(), list(minSize = 1L, jump = 1L, nBkpsMax = NULL,
                                          resolvedNBkpsMax = 20L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat, covariates = NULL,
                                          n = 100, p = 1))

  #Modify non-`costFunc` active bindings

  set.seed(123)
  tsMat2 = cbind(x1 = c(rnorm(50,0), rnorm(50,5)),
                 x2 = c(rnorm(50,0), rnorm(50,5)))
  DynpObj$tsMat = tsMat2
  DynpObj$minSize = 2L
  DynpObj$jump = 2L
  DynpObj$nBkpsMax = 10L

  expect_equal(DynpObj$describe(), list(minSize = 2L, jump = 2L, nBkpsMax = 10L,
                                          resolvedNBkpsMax = 10L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  #Modify `costFunc``
  costFuncObj  = costFunc$new(costFunc = "VAR", pVAR = 1L)
  DynpObj$costFunc = costFuncObj

  expect_equal(DynpObj$describe(T), DynpObj$describe(F))
  expect_equal(DynpObj$describe(), list(minSize = 2L, jump = 2L, nBkpsMax = 10L,
                                          resolvedNBkpsMax = 10L,
                                          costFunc = costFuncObj,
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, pVAR = costFuncObj$pVAR))

  costFuncObj = costFunc$new(costFunc = "SIGMA", pVAR = 1L)
  DynpObj$costFunc = costFuncObj

  expect_equal(DynpObj$describe(T), DynpObj$describe(F))
  expect_equal(DynpObj$describe(), list(minSize = 2L, jump = 2L, nBkpsMax = 10L,
                                          resolvedNBkpsMax = 10L,
                                          costFunc = costFunc$new(costFunc = "SIGMA"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2,
                                          addSmallDiag = costFuncObj$addSmallDiag,
                                          epsilon = costFuncObj$epsilon))

  expect_true(is.null(DynpObj$describe()$pVAR))


  costFuncObj = costFunc$new(costFunc = "L1", pVAR = 1L)
  DynpObj$costFunc = costFuncObj

  expect_equal(DynpObj$describe(T), DynpObj$describe(F))
  expect_equal(DynpObj$describe(), list(minSize = 2L, jump = 2L, nBkpsMax = 10L,
                                          resolvedNBkpsMax = 10L,
                                          costFunc = costFunc$new(costFunc = "L1"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  expect_true(is.null(DynpObj$describe()$pVAR))
  expect_true(is.null(DynpObj$describe()$epsilon))
  expect_true(is.null(DynpObj$describe()$addSmallDiag))


  costFuncObj = costFunc$new(costFunc = "LinearL2", pVAR = 1L)
  expect_warning(DynpObj$costFunc <- costFuncObj)

  expect_equal(DynpObj$describe(T), DynpObj$describe(F))
  expect_equal(DynpObj$describe(), list(minSize = 2L, jump = 2L, nBkpsMax = 10L,
                                          resolvedNBkpsMax = 10L,
                                          costFunc = costFunc$new(costFunc = "LinearL2"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, intercept = TRUE))

  expect_true(is.null(DynpObj$describe()$pVAR))

})

test_that("Error handling for `eval()` works properly", {

  #Cost-specific tests are in other files

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new()
  expect_error(DynpObj$eval(0, 10)) #Not fitted

  DynpObj$fit(tsMat)
  expect_error(DynpObj$eval(NULL, 10))
  expect_error(DynpObj$eval(0, NULL))

})

test_that("Test that error handling for `predict()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new()
  expect_error(DynpObj$predict()) #Not fitted

  DynpObj$fit(tsMat)
  expect_error(DynpObj$predict(NULL))
  expect_error(DynpObj$predict(-1))
  expect_error(DynpObj$predict(c(1:2)))
  expect_error(DynpObj$predict("a"))
  expect_error(DynpObj$predict(T))

  #`nBkps`-specific validation
  expect_error(DynpObj$predict(nBkps = -1), "non-negative integer")
  expect_error(DynpObj$predict(nBkps = 1.5), "non-negative integer")
  expect_error(DynpObj$predict(nBkps = "a"), "non-negative integer")
  expect_error(DynpObj$predict(nBkps = c(1,2)), "non-negative integer")
  expect_no_error(DynpObj$predict(nBkps = 0))

})

test_that("`$predict(nBkps=)` exceeding `nBkpsMax` is capped (message), not an error", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new(nBkpsMax = 5L)
  DynpObj$fit(tsMat)

  expect_message(bkpsCapped <- DynpObj$predict(nBkps = 5 + 50), "exceeds `nBkpsMax`")
  expect_length(bkpsCapped, 6L)
  expect_identical(bkpsCapped, DynpObj$predict(nBkps = 5))

  #`nBkps` takes precedence over `pen` when both are supplied
  expect_identical(DynpObj$predict(pen = 999999, nBkps = 2), DynpObj$predict(nBkps = 2))

})

test_that("`$predict(nBkps=)`/`$costPath()` reflect infeasibility when `jump` is large relative to `minSize`", {

  set.seed(1)
  tsMat = matrix(rnorm(100))
  #grid = {0, 50, 100}: a single interior admissible point -> at most 1 change-point is ever feasible
  DynpObj = Dynp$new(minSize = 1L, jump = 50L)
  expect_message(DynpObj$fit(tsMat), "`nBkpsMax` not set")

  cp = DynpObj$costPath()
  expect_true(is.finite(cp[1])) #k = 0 always feasible
  expect_true(is.finite(cp[2])) #k = 1 feasible: the only admissible split, at 50
  expect_true(all(is.infinite(cp[-(1:2)]))) #k >= 2 infeasible, given only one interior grid point

  expect_equal(DynpObj$predict(nBkps = 1), c(50, 100))
  expect_error(DynpObj$predict(nBkps = 2), "No valid segmentation exists")
  #`pen`-based selection silently skips the infeasible (Inf-cost) entries, never erroring
  expect_equal(DynpObj$predict(pen = 0), c(50, 100))

})

test_that("Error handling for `plot()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new()

  #Without providing endPts
  DynpObj$fit(tsMat)
  expect_error(DynpObj$plot()) #No tmpEndPts available
  DynpObj$predict(nBkps = 1)
  expect_no_error(DynpObj$plot())

  ## Invalid `main`
  expect_error(DynpObj$plot(main = c(1:3)))
  expect_error(DynpObj$plot(main = T))
  expect_error(DynpObj$plot(main = NA))
  expect_error(DynpObj$plot(main = NULL))
  expect_error(DynpObj$plot(main = c("a", "b")))

  ## Invalid `xlab`
  expect_error(DynpObj$plot(xlab = c(1:3)))
  expect_error(DynpObj$plot(xlab = T))
  expect_error(DynpObj$plot(xlab = NA))
  expect_error(DynpObj$plot(xlab = NULL))
  expect_error(DynpObj$plot(xlab = c("a", "b")))

  ## Invalid `d`
  expect_error(DynpObj$plot(d = 2)) #d > p
  expect_error(DynpObj$plot(d = "a"))
  expect_error(DynpObj$plot(d = T))
  expect_error(DynpObj$plot(d = NULL))
  expect_error(DynpObj$plot(d = NA))

  ## Invalid `dimNames`
  expect_error(DynpObj$plot(dimNames = c("X1", "X2")))
  expect_error(DynpObj$plot(dimNames = T))
  expect_error(DynpObj$plot(dimNames = NA))
  expect_error(DynpObj$plot(dimNames = NULL))
  expect_error(DynpObj$plot(main = 123))

  #endPts provided
  expect_error(DynpObj$plot(endPts = "a"))
  expect_error(DynpObj$plot(endPts = NA))
  expect_error(DynpObj$plot(endPts = NULL))
  expect_error(DynpObj$plot(endPts = 1:99)) #endPts not include n
  expect_error(DynpObj$plot(endPts = 0:99)) ##min endPts < 1
  expect_error(DynpObj$plot(endPts = 0:101)) ##max endPts > n
  expect_error(DynpObj$plot(endPts = c(50,50,100))) ##duplicated endpts

})


test_that("Error handling for C++ module Dynp_L2 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_, int nBkpsMax_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::DynpCpp_L2, tsMat, 0, 1, 5)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_L2, tsMat, 1, 0, 5)) #jump_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_L2, tsMat, 10, 3, 5)) #segment too short
  expect_error(new(rupturesRcpp:::DynpCpp_L2, tsMat2, 10, 3, -1)) #nBkpsMax_ < 0
  expect_no_error(new(rupturesRcpp:::DynpCpp_L2, tsMat2, 10, 3, 5)) #len = minLen here

})


test_that("Error handling for C++ module Dynp_L1 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_, int nBkpsMax_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::DynpCpp_L1_cwMed, tsMat, 0, 1, 5)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_L1_cwMed, tsMat, 1, 0, 5)) #jump_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_L1_cwMed, tsMat, 10, 3, 5)) #segment too short
  expect_error(new(rupturesRcpp:::DynpCpp_L1_cwMed, tsMat2, 10, 3, -1)) #nBkpsMax_ < 0
  expect_no_error(new(rupturesRcpp:::DynpCpp_L1_cwMed, tsMat2, 10, 3, 5)) #len = minLen here

})


test_that("Error handling for C++ module Dynp_VAR works as intended", {
  #constructor: const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int nBkpsMax_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  #pVAR = 1
  expect_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat, 1, 0, 1, 5)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat, 1, 1, 0, 5)) #jump_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat, 1, 10, 3, 5)) #segment too short
  expect_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat2, 1, 10, 3, -1)) #nBkpsMax_ < 0
  expect_no_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat2, 1, 10, 3, 5)) #len = minLen here

  set.seed(123)
  tsMat3 = cbind(x1 = rnorm(5), x2 = rnorm(5), x3 = rnorm(5))
  expect_error(new(rupturesRcpp:::DynpCpp_VAR, tsMat3, 2, 1, 1, 5)) #segment too short for fitting VAR(2)
})


test_that("Error handling for C++ module Dynp_SIGMA works as intended", {
  #constructor: const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_, int nBkpsMax_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::DynpCpp_SIGMA, tsMat, T, 10^-6, 0, 1, 5)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_SIGMA, tsMat, T, 10^-6, 1, 0, 5)) #jump_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_SIGMA, tsMat, T, 10^-6, 10, 3, 5)) #segment too short
  expect_error(new(rupturesRcpp:::DynpCpp_SIGMA, tsMat2, T, 10^-6, 10, 3, -1)) #nBkpsMax_ < 0
  expect_no_error(new(rupturesRcpp:::DynpCpp_SIGMA, tsMat2, T, 10^-6, 10, 3, 5)) #len = minLen here

})


test_that("Error handling for C++ module Dynp_LinearL2 works as intended", {
  #constructor: const arma::mat& tsMat,  const arma::mat& covariates, bool intercept_, int minSize_, int jump_, int nBkpsMax_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))
  covariateMat = as.matrix(rep(1,23))
  covariateMat2 = as.matrix(rep(1,24))

  expect_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat, covariateMat, T,  0, 1, 5)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat, covariateMat, T,  1, 0, 5)) #jump_ = 0
  expect_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat, covariateMat, T, 10, 3, 5)) #segment too short
  expect_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat2, covariateMat2, T, 10, 3, -1)) #nBkpsMax_ < 0
  expect_no_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat2, covariateMat2, T, 10, 3, 5)) #len = minLen here

  set.seed(1234)
  tsMat3 = as.matrix(rnorm(5))
  covariateMat3 = cbind(rnorm(5),rnorm(5),rnorm(5),rnorm(5),rnorm(5))
  expect_error(new(rupturesRcpp:::DynpCpp_LinearL2, tsMat3, covariateMat3, T, 1, 1, 2)) #too short for fitting linear regression model
})

test_that("Some additional tests", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new()
  DynpObj$fit(tsMat)
  expect_equal(DynpObj$predict(999999), 100) #Too large `pen` -> only return `n`

})

test_that("`$getHistory()`/`$costPath()` return a well-formed, non-increasing cost history", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new(minSize = 2L, nBkpsMax = 10L)

  expect_error(DynpObj$getHistory(), "must be run before") #not fitted yet
  expect_error(DynpObj$costPath(), "must be run before")

  DynpObj$fit(tsMat)
  hist = DynpObj$getHistory()
  cp = DynpObj$costPath()

  expect_s3_class(hist, "data.frame")
  expect_identical(names(hist), c("k", "cost")) #no `added_bkp` column, unlike binSeg/Window
  expect_equal(hist$k, 0:10)
  expect_equal(hist$cost, cp)
  expect_length(cp, 11L) #resolvedNBkpsMax + 1

  #cost is non-increasing as k grows: more change-points can never hurt
  expect_true(all(diff(hist$cost) <= 1e-8))

  #cost at k=0 is the whole-series cost
  expect_equal(hist$cost[1], DynpObj$eval(0, 100))

  #the cost reported for each k actually matches re-evaluating $predict(nBkps = k)'s segmentation
  for(k in 0:10){
    bkps = DynpObj$predict(nBkps = k)
    expect_equal(segCost(DynpObj, bkps), cp[k + 1], tolerance = 1e-8)
  }

})

test_that("`$plotElbow()` returns a ggplot object and respects `maxK`", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new(minSize = 2L, nBkpsMax = 10L)

  expect_error(DynpObj$plotElbow(), "must be run before") #not fitted yet

  DynpObj$fit(tsMat)

  expect_no_error(DynpObj$plotElbow())
  p = DynpObj$plotElbow()
  expect_s3_class(p, "ggplot")

  pCapped = DynpObj$plotElbow(maxK = 2)
  expect_equal(nrow(pCapped$data), 3) #k = 0,1,2

  expect_error(DynpObj$plotElbow(maxK = 0), "positive integer")
  expect_error(DynpObj$plotElbow(maxK = "a"), "positive integer")

})

test_that("`$predict(nBkps=)` returns a valid, independently-exact segmentation for every k", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  DynpObj = Dynp$new(minSize = 2L, nBkpsMax = 5L)
  DynpObj$fit(tsMat)

  for(k in 0:5){
    bkps = DynpObj$predict(nBkps = k)
    expect_length(bkps, k + 1L) #k change-points + n
    expect_equal(tail(bkps, 1), 100)
    expect_equal(bkps, sort(unique(bkps))) #strictly ascending, no duplicates
    expect_true(all(diff(c(0, bkps)) >= 2)) #respects minSize = 2
  }

  #Existing `pen`-only behaviour is unaffected
  expect_equal(DynpObj$predict(999999), 100)

  #Validation
  expect_error(DynpObj$predict(nBkps = -1), "non-negative integer")
  expect_error(DynpObj$predict(nBkps = 1.5), "non-negative integer")
  expect_error(DynpObj$predict(nBkps = "a"), "non-negative integer")
  expect_no_error(DynpObj$predict(nBkps = 0))

})


# ========================================================
#  Dynp is *exact*: brute-force verification for small n
# ========================================================
# Unlike binSeg (greedy) or Window (local maxima of a sliding gain), Dynp's whole purpose is to
# be exact for a specified number of change-points. For small n, we can verify this directly by
# exhaustively enumerating every (minSize, jump)-admissible segmentation and confirming Dynp's
# reported minimal cost, and the cost of the segmentation it actually returns, match the
# brute-force minimum -- keeping n and nBkpsMax small here so the exhaustive search (and DP table)
# stay fast, per Dynp's O(nBkpsMax * M^2) complexity.

bruteForceMinCost = function(DynpObj, n, minSize, jump, k) {

  if(k == 0){
    return(DynpObj$eval(0, n))
  }

  admissible = seq(jump, n - 1, by = jump)
  admissible = admissible[admissible >= minSize]

  if(length(admissible) < k){
    return(Inf)
  }

  combos = combn(admissible, k)
  best = Inf

  for(col in seq_len(ncol(combos))){
    bkps = c(combos[, col], n)
    starts = c(0, head(bkps, -1))

    if(any(bkps - starts < minSize)){
      next
    }

    cost = segCost(DynpObj, bkps)
    if(cost < best){
      best = cost
    }
  }

  best
}

test_that("Dynp matches brute-force search over all admissible segmentations (L1, L2)", {

  set.seed(2024)
  n = 14
  nBkpsMaxTest = 3 #keeps choose(n-1, k) and the DP table cheap

  configs = list(
    list(minSize = 1L, jump = 1L),
    list(minSize = 2L, jump = 1L),
    list(minSize = 1L, jump = 2L)
  )

  for(cfName in c("L1", "L2")){

    tsMat = matrix(rnorm(n))

    for(cfg in configs){

      DynpObj = Dynp$new(minSize = cfg$minSize, jump = cfg$jump, nBkpsMax = nBkpsMaxTest,
                          costFunc = costFunc$new(cfName))
      DynpObj$fit(tsMat)
      cp = DynpObj$costPath()

      for(k in 0:nBkpsMaxTest){

        bf = bruteForceMinCost(DynpObj, n, cfg$minSize, cfg$jump, k)
        expect_equal(cp[k + 1], bf, tolerance = 1e-8)

        if(is.finite(bf)){
          #The segmentation `$predict(nBkps = k)` returns must itself attain that cost
          bkps = DynpObj$predict(nBkps = k)
          expect_equal(segCost(DynpObj, bkps), bf, tolerance = 1e-8)
        }
      }
    }
  }

})

test_that("Dynp achieves cost <= binSeg's greedy solution for the same `nBkps`", {

  set.seed(99)
  tsMat = matrix(c(rnorm(40,0), rnorm(40,3), rnorm(40,-2)))
  kMax = 6L

  DynpObj = Dynp$new(nBkpsMax = kMax)
  DynpObj$fit(tsMat)

  binSegObj = binSeg$new()
  binSegObj$fit(tsMat)

  for(k in 1:kMax){
    dynpCost = segCost(DynpObj, DynpObj$predict(nBkps = k))
    binSegCost = segCost(binSegObj, binSegObj$predict(nBkps = k))

    expect_lte(dynpCost, binSegCost + 1e-8)
  }

})

test_that("`$segments()` returns the cost and params of each segment", {

  set.seed(1)
  X = matrix(c(rnorm(100, 0), rnorm(100, 5)))
  DynpObj = Dynp$new(nBkpsMax = 5L) #Explicit `nBkpsMax`: no "not set" message
  DynpObj$fit(X)
  expect_error(DynpObj$segments(), "Must run") #No `$predict()` yet

  bkps = DynpObj$predict(pen = 10)
  segs = DynpObj$segments()

  expect_length(segs, length(bkps))
  expect_equal(sapply(segs, `[[`, "Start"), c(0, head(bkps, -1)))
  expect_equal(sapply(segs, `[[`, "End"), bkps)

  for (s in segs) {
    Xe = X[(s$Start + 1):s$End, , drop = FALSE]
    expect_equal(s$Cost, sum((Xe - mean(Xe))^2))
    expect_equal(s$Params$mean, colMeans(Xe))
  }

  #Costs sum to the exact minimal cost for the number of change-points returned, for both `$predict()` modes
  cp = DynpObj$costPath()
  expect_equal(sum(sapply(segs, `[[`, "Cost")), cp[length(bkps)])

  for (k in 0:3) {
    bkps = DynpObj$predict(nBkps = k)
    segs = DynpObj$segments()
    expect_equal(sapply(segs, `[[`, "End"), bkps)
    expect_equal(sum(sapply(segs, `[[`, "Cost")), cp[k + 1])
  }

  #The intercept-only force-fit path returns early from `$fit()`, but must still clear stale end points
  linObj = Dynp$new(nBkpsMax = 5L, costFunc = costFunc$new("LinearL2"))
  expect_warning(linObj$fit(X), "an intercept")
  linObj$predict(pen = 10)
  expect_warning(linObj$fit(), "an intercept")
  expect_error(linObj$segments(), "Must run")

})

test_that("`$segments()` agrees with `costFactory` for every cost function", {

  set.seed(2)
  X = matrix(c(rnorm(100, 0), rnorm(100, 5, 3)))
  covariates = matrix(rnorm(200))
  DynpObj = Dynp$new(minSize = 5L, nBkpsMax = 5L)
  DynpObj$fit(X)
  DynpObj$covariates = covariates #`$fit()` only keeps `covariates` for Linear costs

  for (cf in c("L1", "L2", "SIGMA", "VAR", "LinearL2", "LinearSIGMA")) {

    suppressMessages(DynpObj$costFunc <- costFunc$new(cf))
    expect_error(DynpObj$segments(), "Must run") #Refitting clears stale end points

    DynpObj$predict(pen = 10)
    facObj = costFactory$new(costFunc$new(cf))
    facObj$fit(X, covariates)

    for (s in DynpObj$segments()) {
      expect_equal(s$Cost, facObj$eval(s$Start, s$End))
      expect_equal(s$Params, facObj$get_params(s$Start, s$End))
    }
  }

})
