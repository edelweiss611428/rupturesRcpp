# costFactory

set.seed(12345)
n = 120; p = 2
tsMat = cbind(as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive")),
              as.numeric(stats::filter(rnorm(n), -0.5, method = "recursive")))
X = cbind(rnorm(n), rnorm(n))
Ylin = X %*% matrix(c(1, -2, 0.5, 3), 2) + rnorm(2 * n)          # Y = X B + noise, p = 2

L2eval = function(segment, a, b) sum(sweep(segment, 2, colMeans(segment))^2)
L2params = function(segment, a, b) colMeans(segment)

# costFunc, data, covariates, get_params() fields
cases = list(
  list(costFunc$new("L1"), tsMat, NULL, "median"),
  list(costFunc$new("L2"), tsMat, NULL, "mean"),
  list(costFunc$new("SIGMA"), tsMat, NULL, c("mean", "cov")),
  list(costFunc$new("VAR", pVAR = 2), tsMat, NULL, "coef"),
  list(costFunc$new("LinearL2"), Ylin, X, "coef"),
  list(costFunc$new("LinearSIGMA"), Ylin, X, c("coef", "cov")),
  list(costFunc$new("LinearL1"), Ylin, X, "coef"),
  list(costFunc$new("Custom", evalFun = L2eval, paramFun = L2params), tsMat, NULL, "params")
)

# random segments of length >= 20 (full rank for every cost)
segs = t(replicate(8, { a = sample(0:(n - 20), 1); c(a, a + 19L + sample.int(n - a - 19L, 1)) }))

test_that("costFactory validates inputs", {
  expect_error(costFactory$new("L2", tsMat), "must be a `R6` object")
  expect_error(costFactory$new(costFunc$new(), "a"), "numeric time series matrix")
  expect_error(costFactory$new(costFunc$new(), cbind(c(1, NA))), "contains NAs")
  expect_error(costFactory$new(costFunc$new()), "must be provided")
  expect_error(costFactory$new(costFunc$new("LinearL2"), tsMat, X[1:10, ]), "do not match")
  expect_error(costFactory$new(costFunc$new("Custom"), tsMat), "requires `evalFun`")
  expect_named(costFactory$new(tsMat = tsMat)$get_params(0, 10), "mean")                   # default costFunc: L2
  cf = costFactory$new(costFunc$new(), tsMat)
  expect_error(cf$eval(-1, 10), "0 <= start")
  expect_error(cf$eval(0, n + 1), "0 < end")
  expect_error(cf$eval(10, 10), "smaller than b")
  expect_error(cf$get_params(10, 10), "smaller than b")
})

test_that("$eval() and $get_params() match the PELT module for every cost function", {
  for (cs in cases) {
    cf = costFactory$new(cs[[1]], cs[[2]], cs[[3]])
    P = PELT$new(costFunc = cs[[1]]); suppressWarnings(P$fit(cs[[2]], cs[[3]]))
    mod = P$.__enclos_env__$private$.PELTModule
    for (i in seq_len(nrow(segs))) {
      a = segs[i, 1]; b = segs[i, 2]
      expect_equal(cf$eval(a, b), P$eval(a, b))
      expect_identical(cf$get_params(a, b), mod$get_params(a, b))
      expect_named(cf$get_params(a, b), cs[[4]])
    }
  }
})

test_that("Linear costs force-fit an intercept without covariates", {
  for (cfn in c("LinearL2", "LinearSIGMA", "LinearL1")) {
    expect_warning(cf <- costFactory$new(costFunc$new(cfn), Ylin), "Force-fitting")
    expect_equal(dim(cf$get_params(0, 50)$coef), c(1L, p))
  }
  expect_warning(cf <- costFactory$new(costFunc$new("LinearL2"), Ylin), "Force-fitting")
  expect_equal(c(cf$get_params(0, 50)$coef), colMeans(Ylin[1:50, ]))
})

test_that("VAR: segments shorter than the number of coefficients give NA coef and zero cost", {
  cf = costFactory$new(costFunc$new("VAR", pVAR = 2), tsMat)
  e = cf$get_params(0, 3)                                          # n < J = 5
  expect_true(all(is.na(e$coef)))
  expect_equal(dim(e$coef), c(5L, 2L))
  expect_equal(cf$eval(0, 3), 0)
})

test_that("Custom cost without paramFun returns an empty list", {
  cf = costFactory$new(costFunc$new("Custom", evalFun = L2eval), tsMat)
  expect_identical(cf$get_params(0, 50), list())
})
