# Level offsets: the cumulative sums are built on centred data, so costs and estimates must not lose
# precision when a series sits far from zero.

set.seed(1)
n = 2000; a = 1000; b = 1100; off = 1e7
x = cbind(as.numeric(stats::filter(rnorm(n), 0.5, method = "recursive")), rnorm(n))
X = matrix(seq_len(n) / n)
Y = cbind(1 + 2 * X + rnorm(n, sd = 0.1), -1 + X + rnorm(n, sd = 0.1))

test_that("costs are invariant to a level offset", {
  for (cfn in c("L2", "SIGMA", "VAR")) {
    expect_equal(costFactory$new(costFunc$new(cfn), x + off)$eval(a, b),
                 costFactory$new(costFunc$new(cfn), x)$eval(a, b), tolerance = 1e-6)
  }
  for (cfn in c("LinearL2", "LinearSIGMA")) {
    expect_equal(costFactory$new(costFunc$new(cfn), Y + off, X + off)$eval(a, b),
                 costFactory$new(costFunc$new(cfn), Y, X)$eval(a, b), tolerance = 1e-6)
    expect_equal(suppressWarnings(costFactory$new(costFunc$new(cfn), Y + off))$eval(a, b),   # force-fit: ones column
                 suppressWarnings(costFactory$new(costFunc$new(cfn), Y))$eval(a, b), tolerance = 1e-6)
  }
})

test_that("estimates are returned on the original scale", {
  s0 = costFactory$new(costFunc$new("SIGMA"), x)$get_params(a, b)
  s1 = costFactory$new(costFunc$new("SIGMA"), x + off)$get_params(a, b)
  expect_equal(s1$mean, s0$mean + off)
  expect_equal(s1$cov, s0$cov, tolerance = 1e-6)
  B0 = costFactory$new(costFunc$new("LinearL2"), Y, X)$get_params(a, b)$coef
  B1 = costFactory$new(costFunc$new("LinearL2"), Y + off, X + off)$get_params(a, b)$coef
  expect_equal(B1[2, ], B0[2, ], tolerance = 1e-6)                        # slopes
  expect_equal(B1[1, ], B0[1, ] + off - off * B0[2, ], tolerance = 1e-6)  # intercept of y + c on x + d: b0 + c - d * b1
  V0 = costFactory$new(costFunc$new("VAR"), x)$get_params(a, b)$coef
  V1 = costFactory$new(costFunc$new("VAR"), x + off)$get_params(a, b)$coef
  expect_equal(V1[-1, ], V0[-1, ], tolerance = 1e-6)                      # lag coefficients
})
