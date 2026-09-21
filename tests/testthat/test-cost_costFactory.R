# costFactory + segSummary

set.seed(12345)
n = 120; p = 2
tsMat = cbind(as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive")),
              as.numeric(stats::filter(rnorm(n), -0.5, method = "recursive")))
X = cbind(rnorm(n), rnorm(n))
Ylin = X %*% matrix(c(1, -2, 0.5, 3), 2) + rnorm(2 * n)          # Y = X B + noise, p = 2

# random segments of length >= 20 (full rank for every cost)
segs = t(replicate(8, { a = sample(0:(n - 20), 1); c(a, sample((a + 20):n, 1)) }))

test_that("costFactory validates inputs", {
  expect_error(costFactory$new("L2", tsMat), "must be a `R6` object")
  expect_error(costFactory$new(costFunc$new(), "a"), "numeric time series matrix")
  expect_error(costFactory$new(costFunc$new(), cbind(c(1, NA))), "contains NAs")
  expect_error(costFactory$new(costFunc$new()), "must be provided")
  expect_error(costFactory$new(costFunc$new("LinearL2"), tsMat, X[1:10, ]), "do not match")
  cf = costFactory$new(costFunc$new(), tsMat)
  expect_error(cf$eval(-1, 10), "0 <= start")
  expect_error(cf$eval(0, n + 1), "0 < end")
  expect_error(cf$eval(10, 10), "smaller than b")
})

test_that("costFactory$eval matches PELT$eval", {
  for (cfn in c("L2", "VAR")) {
    cf = costFactory$new(costFunc$new(cfn), tsMat)
    P = PELT$new(costFunc = costFunc$new(cfn)); P$fit(tsMat)
    for (i in seq_len(nrow(segs))) expect_equal(cf$eval(segs[i, 1], segs[i, 2]), P$eval(segs[i, 1], segs[i, 2]))
  }
})

test_that("L1 / L2 / SIGMA estimates decompose the cost", {
  cf1 = costFactory$new(costFunc$new("L1"), tsMat)
  cf2 = costFactory$new(costFunc$new("L2"), tsMat)
  cfS = costFactory$new(costFunc$new("SIGMA"), tsMat)
  cfS0 = costFactory$new(costFunc$new("SIGMA", addSmallDiag = FALSE), tsMat)
  for (i in seq_len(nrow(segs))) {
    a = segs[i, 1]; b = segs[i, 2]; len = b - a
    e1 = cf1$estimate(a, b); expect_equal(sum(e1$scale) * len, cf1$eval(a, b))
    e2 = cf2$estimate(a, b); expect_equal(sum(e2$var) * len, cf2$eval(a, b))
    eS = cfS$estimate(a, b)
    expect_equal(len * log(det(eS$cov + diag(1e-6, p))), cfS$eval(a, b), tolerance = 1e-8)
    expect_equal(len * log(det(eS$cov)), cfS0$eval(a, b), tolerance = 1e-8)
    expect_named(e2$mean, c("X1", "X2"))
  }
})

test_that("VAR estimates decompose the cost, including the first segment", {
  cf = costFactory$new(costFunc$new("VAR", pVAR = 2), tsMat)
  for (seg in c(list(c(0, 40)), split(segs, seq_len(nrow(segs))))) {
    a = seg[1]; b = seg[2]
    e = cf$estimate(a, b)
    nEff = b - max(a, 2)
    expect_equal(sum(e$var) * nEff, cf$eval(a, b))
    expect_equal(rownames(e$coef), c("(Intercept)", "X1.lag1", "X2.lag1", "X1.lag2", "X2.lag2"))
    expect_equal(colnames(e$coef), c("X1", "X2"))
  }
  e = cf$estimate(0, 3)                                            # n < J = 5
  expect_true(all(is.na(e$coef))); expect_equal(dim(e$coef), c(5L, 2L)); expect_equal(cf$eval(0, 3), 0)
})

test_that("LinearL2 estimates decompose the cost, with and without covariates", {
  cf = costFactory$new(costFunc$new("LinearL2"), Ylin, X)
  for (i in seq_len(nrow(segs))) {
    a = segs[i, 1]; b = segs[i, 2]
    e = cf$estimate(a, b)
    expect_equal(sum(e$var) * (b - a), cf$eval(a, b))
    expect_equal(rownames(e$coef), c("(Intercept)", "Z1", "Z2"))
  }
  expect_warning(cf0 <- costFactory$new(costFunc$new("LinearL2"), Ylin), "Force-fitting")
  e0 = cf0$estimate(0, 50)
  expect_equal(rownames(e0$coef), "(Intercept)")
  expect_equal(drop(e0$coef), colMeans(Ylin[1:50, ]), ignore_attr = TRUE)
  expect_equal(sum(e0$var) * 50, cf0$eval(0, 50))
})

test_that("summary() returns a printable, extractable segSummary", {
  cf = costFactory$new(costFunc$new("SIGMA"), tsMat)
  endPts = c(40L, 90L, n)
  s = cf$summary(endPts)
  expect_s3_class(s, "segSummary")
  expect_length(s, 3L)
  expect_named(s[[1]], c("start", "end", "n", "cost", "mean", "cov"))
  expect_equal(vapply(s, function(x) x$end, integer(1)), endPts)
  expect_equal(s[[2]]$cost, cf$eval(40, 90))
  expect_identical(s[[2]][c("mean", "cov")], cf$estimate(40, 90))
  expect_error(cf$summary(c(40, 90)), "must be `n`")
  expect_error(cf$summary(c(0, n)), "less than 1")
  expect_error(cf$summary(c(40, 40, n)), "duplicated")
  expect_output(print(s), "Segment 2: \\(40, 90\\]")
  expect_output(print(s), "epsilon = 1e-06")
  df = as.data.frame(s)
  expect_equal(nrow(df), 3L)
  expect_equal(df$cost, vapply(s, function(x) x$cost, numeric(1)))
  expect_equal(df$mean.X1[2], s[[2]]$mean[["X1"]])
  expect_true("cov.X1.X2" %in% names(df))
  dfv = as.data.frame(costFactory$new(costFunc$new("VAR"), tsMat)$summary(endPts))
  expect_true("coef.(Intercept).X1" %in% names(dfv))
})
