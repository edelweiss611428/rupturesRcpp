#test-createCostFunc.R

test_that("Test for error", {
  expect_error(createCostFunc("Sydney"), regexp = "not supported")
  expect_error(createCostFunc("SIGMA", addSmallDiag = 1111), regexp = "single boolean")
  expect_error(createCostFunc("SIGMA", addSmallDiag = T, epsilon = T),
               regexp = "single non-negative")
  expect_error(createCostFunc("SIGMA", addSmallDiag = T, epsilon = -1),
               regexp = "single non-negative")
  expect_error(createCostFunc("SIGMA", addSmallDiag = c(T,T), epsilon = -1),
               regexp = "single boolean")
  expect_error(createCostFunc("SIGMA", addSmallDiag = T, epsilon = -1),
               regexp = "single non-negative")
  expect_error(createCostFunc("SIGMA", addSmallDiag = T, epsilon = c(1,1)),
               regexp = "single non-negative")
  expect_error(createCostFunc("VAR", pVAR = 0.5),
               regexp = "single positive integer")
})


test_that("Expect no error", {
  expect_no_error(createCostFunc()) #L2 by default
  expect_no_error(createCostFunc("SIGMA"))
  expect_no_error(createCostFunc("VAR"))
})


test_that("Correct params (L2)", {
  L2obj = createCostFunc()
  expect_equal(L2obj$costFunc, "L2")
})


test_that("Correct params (VAR)", {
  VARobj = createCostFunc("VAR")  #default
  expect_equal(VARobj$costFunc, "VAR")
  expect_equal(VARobj$pVAR, 1L)

  VARobj = createCostFunc("VAR", pVAR = 5L, epsilon = 1e-6)
  expect_equal(VARobj$costFunc, "VAR")
  expect_equal(VARobj$pVAR, 5L)
  expect_null(VARobj$epsilon)
})


test_that("Correct params (SIGMA)", {
  SIGMAobj = createCostFunc("SIGMA") #default
  expect_equal(SIGMAobj$costFunc, "SIGMA")
  expect_equal(SIGMAobj$addSmallDiag, T)
  expect_equal(SIGMAobj$epsilon, 1e-6)

  SIGMAobj = createCostFunc("SIGMA", addSmallDiag = F, epsilon = 1e-5, pVAR = 2L)
  expect_equal(SIGMAobj$costFunc, "SIGMA")
  expect_equal(SIGMAobj$addSmallDiag, F)
  expect_equal(SIGMAobj$epsilon, 1e-5)
  expect_null(SIGMAobj$pVAR)
})



