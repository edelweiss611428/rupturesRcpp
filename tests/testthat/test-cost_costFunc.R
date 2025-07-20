#costFunc module
set.seed(12345)


test_that("Cost function not supported", {

  expect_error(costFunc$new("Sydney"), regexp = "not supported")
  expect_error(costFunc$new(c("L1", "L2")), regexp = "must be a single character value!")

})


test_that("Expect that SIGMA cost module gives correct error messages", {

  #addSmallDiag must be a single boolean value

  expect_error(costFunc$new("SIGMA", addSmallDiag = "abc"), regexp = "single boolean")
  expect_error(costFunc$new("SIGMA", addSmallDiag = 1111), regexp = "single boolean")
  expect_error(costFunc$new("SIGMA", addSmallDiag = c(T,T), epsilon = -1),
               regexp = "single boolean")

  #epsilon must be a single positive value
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = "a"),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = T),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = -1),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = 0),
               regexp = "single positive")
  expect_error(costFunc$new("SIGMA", addSmallDiag = T, epsilon = c(1,1)),
               regexp = "single positive")

})

test_that("Expect that VAR cost module gives correct error messages", {

  expect_error(costFunc$new("VAR", pVAR = T),
               regexp = "single positive integer")

  expect_warning(expect_error(costFunc$new("VAR", pVAR = "a"),
                              regexp = "single positive integer")) #warning due to as.integer("a")

  expect_error(costFunc$new("VAR", pVAR = 0.5),
               regexp = "single positive integer")

})

test_that("Expect that LinearL2 cost module gives correct error messages", {


  expect_error(costFunc$new("LinearL2", intercept = "a"),
               regexp = "single boolean")

  expect_error(costFunc$new("LinearL2", intercept = 1L),
               regexp = "single boolean")

})


test_that("Expect default modules give no error when initialising", {

  expect_no_error(costFunc$new())
  expect_no_error(costFunc$new("SIGMA"))
  expect_no_error(costFunc$new("VAR"))
  expect_no_error(costFunc$new("LinearL2"))

})


test_that("Correct params (L2)", {

  costFuncObj = costFunc$new()
  expect_equal(costFuncObj$costFunc, "L2")
  expect_equal(costFuncObj$pass()$costFunc, "L2")

  #pass() only gives required params, not others
  costFuncObj = costFunc$new("L2", addSmallDiag = T, epsilon = 10^-5, pVAR = 1L)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), "costFunc"))

})

test_that("Modifying active binding `costFunc` gives the expected default configurations", {

  costFuncObj = costFunc$new()
  costFuncObj$costFunc = "VAR"
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "pVAR")))
  expect_equal(costFuncObj$pVAR, 1L)

  costFuncObj$costFunc = "SIGMA"
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "addSmallDiag", "epsilon")))
  expect_equal(costFuncObj$addSmallDiag, T)
  expect_equal(costFuncObj$epsilon, 10^-6)

  costFuncObj$costFunc = "LinearL2"
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "intercept")))
  expect_equal(costFuncObj$intercept, T)

})


test_that("Correct params (VAR)", {

  costFuncObj = costFunc$new("VAR")  #default configurations
  expect_equal(costFuncObj$costFunc, "VAR")
  expect_equal(costFuncObj$pVAR, 1L)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "pVAR")))

  costFuncObj = costFunc$new("VAR", pVAR = 5L, epsilon = 1e-6)
  expect_equal(costFuncObj$costFunc, "VAR")
  expect_equal(costFuncObj$pVAR, 5L)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "pVAR")))

  #Modifying active binding `pVAR`
  costFuncObj = costFunc$new("VAR")
  costFuncObj$pVAR = 10L
  expect_equal(costFuncObj$pVAR, 10L)

})


test_that("Correct params (SIGMA)", {

  costFuncObj = costFunc$new("SIGMA") #default configurations
  expect_equal(costFuncObj$costFunc, "SIGMA")
  expect_equal(costFuncObj$addSmallDiag, T)
  expect_equal(costFuncObj$epsilon, 1e-6)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "epsilon", "addSmallDiag")))

  costFuncObj = costFunc$new("SIGMA", addSmallDiag = F, epsilon = 1e-5, pVAR = 2L)
  expect_equal(costFuncObj$costFunc, "SIGMA")
  expect_equal(costFuncObj$addSmallDiag, F)
  expect_equal(costFuncObj$epsilon, 1e-5)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "epsilon", "addSmallDiag")))

  #Modifying active bindings
  costFuncObj = costFunc$new("SIGMA")
  costFuncObj$epsilon = 0.5
  expect_equal(costFuncObj$epsilon, 0.5)
  costFuncObj$addSmallDiag = F
  expect_equal(costFuncObj$addSmallDiag, F)

})


test_that("Correct params (LinearL2)", {

  costFuncObj = costFunc$new("LinearL2") #default configurations
  expect_equal(costFuncObj$costFunc, "LinearL2")
  expect_equal(costFuncObj$intercept, T)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "intercept")))


  costFuncObj = costFunc$new("LinearL2", intercept = FALSE)
  expect_equal(costFuncObj$costFunc, "LinearL2")
  expect_equal(costFuncObj$intercept, F)
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "intercept")))

  #Modifying active binding `intercept`
  costFuncObj = costFunc$new("LinearL2")
  costFuncObj$intercept = F
  expect_equal(costFuncObj$intercept, F)

})



