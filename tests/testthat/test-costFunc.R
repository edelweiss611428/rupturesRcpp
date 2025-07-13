#test-costFunc.R
set.seed(12345)
test_that("Test for error", {

  expect_error(costFunc$new("Sydney"), regexp = "not supported")
  expect_error(costFunc$new("SIGMA", addSmallDiag = 1111), regexp = "single boolean")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = T),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = -1),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = c(T,T), epsilon = -1),
               regexp = "single boolean")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = 0),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = c(1,1)),
               regexp = "single positive")
  expect_error(costFunc$new("VAR", pVAR = 0.5),
               regexp = "single positive integer")

})


test_that("Expect no error", {

  expect_no_error(costFunc$new()) #L2 by default
  expect_no_error(costFunc$new("SIGMA"))
  expect_no_error(costFunc$new("VAR"))

})


test_that("Correct params (L2)", {

  L2obj = costFunc$new()
  expect_equal(L2obj$costFunc, "L2")

})


test_that("Correct params (VAR)", {

  VARobj = costFunc$new("VAR")  #default
  expect_equal(VARobj$costFunc, "VAR")
  expect_equal(VARobj$pVAR, 1L)
  expect_null(VARobj$epsilon)
  expect_null(VARobj$addSmallDiag)

  VARobj = costFunc$new("VAR", pVAR = 5L, epsilon = 1e-6)
  expect_equal(VARobj$costFunc, "VAR")
  expect_equal(VARobj$pVAR, 5L)
  expect_null(VARobj$epsilon)
  expect_null(VARobj$addSmallDiag)

})


test_that("Correct params (SIGMA)", {

  SIGMAobj = costFunc$new("SIGMA") #default
  expect_equal(SIGMAobj$costFunc, "SIGMA")
  expect_equal(SIGMAobj$addSmallDiag, T)
  expect_equal(SIGMAobj$epsilon, 1e-6)
  expect_null(SIGMAobj$VAR)

  SIGMAobj = costFunc$new("SIGMA", addSmallDiag = F, epsilon = 1e-5, pVAR = 2L)
  expect_equal(SIGMAobj$costFunc, "SIGMA")
  expect_equal(SIGMAobj$addSmallDiag, F)
  expect_equal(SIGMAobj$epsilon, 1e-5)
  expect_null(SIGMAobj$pVAR)

})


test_that("Correct params (LinearL2)", {

  LinearL2obj = costFunc$new("LinearL2") #default
  expect_equal(LinearL2obj$costFunc, "LinearL2")
  expect_equal(LinearL2obj$intercept, T)
  expect_null(LinearL2obj$VAR)

  LinearL2obj = costFunc$new("LinearL2", intercept = FALSE)
  expect_equal(LinearL2obj$costFunc, "LinearL2")
  expect_equal(LinearL2obj$intercept, F)
  expect_null(LinearL2obj$pVAR)

})



