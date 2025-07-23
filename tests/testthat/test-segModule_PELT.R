#PELT module

X_constantSeg = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

test_that("PELT_L2 modules in `rupturesRcpp` and `ruptures` give the same results", {

  set.seed(12345)
  X_noBkp = matrix(rnorm(200))

  skip_if_not_installed("reticulate")
  reticulate = asNamespace("reticulate")

  if (!reticulate::py_module_available("numpy") ||
      !reticulate::py_module_available("ruptures")) {
    skip("Required Python modules not available")
  }

  ruptures = reticulate::import("ruptures")
  np = reticulate::import("numpy")
  npX_noBkp = np$array(X_noBkp)

  #PyPeltL2(min_size = 1L, jump = 1L)
  PyPelt = ruptures$Pelt(model = "l2", min_size = 1L, jump = 1L)
  PyPelt$fit(npX_noBkp)

  #RPeltL2(min_size = 1L, jump = 1L)
  RPelt = PELT$new(minSize = 1L, jump = 1L, costFunc = costFunc$new("L2"))
  RPelt$fit(X_noBkp)

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 3L, jump = 3L)
  PyPelt$min_size = 3L
  PyPelt$jump = 3L

  #RPeltL2(min_size = 3L, jump = 3L)
  RPelt$minSize = 3L
  RPelt$jump = 3L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 5L, jump = 5L)
  PyPelt$min_size = 5L
  PyPelt$jump = 5L

  #RPeltL2(min_size = 5L, jump = 5L)
  RPelt$minSize = 5L
  RPelt$jump = 5L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

})



test_that("PELT_L1 works for constant segments", {

  costFuncObj = costFunc$new("L1")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("PELT_L2 works for constant segments", {

  costFuncObj = costFunc$new("L2")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("PELT_SIGMA works for constant segments", {

  costFuncObj = costFunc$new("SIGMA")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})


test_that("PELT_VAR works for constant segments", {

  costFuncObj = costFunc$new("VAR")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_warning(expect_equal(PeltObj$predict(pen = 0.1), seq(50,150,50)),
                 "Some systems seem singular!")
  expect_warning(expect_false(PeltObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(PeltObj$eval(0,50), 0)), "seems singular")

})

test_that("PELT_LinearL2 works for constant segments", {

  costFuncObj = costFunc$new("LinearL2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_warning(PELTObj$fit(X_constantSeg), "No `covariates` found!")

  expect_equal(PELTObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(PELTObj$eval(0,51) == 0)
  expect_true(all.equal(PELTObj$eval(0,50), 0))

})


test_that("Test if active bindings work properly (i.e., setter/getter/input validating abilility)", {

  X1 = matrix(rep(0,100))
  X2 = matrix(1:100)
  #L2
  costFuncObj = costFunc$new("L2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  PELTObj$tsMat = X2
  expect_no_error(PELTObj$fit())
  expect_true(all.equal(PELTObj$tsMat, X2))
  expect_true(all.equal(PELTObj$eval(0,2), 0.5)) #Setting tsMat via active bindings should work
  PELTObj$tsMat = X1
  expect_true(all.equal(PELTObj$eval(0,2), 0)) #Expect refitted after PELTObj$tsMat = X1

  #VAR
  PELTObj$tsMat = X2
  costFuncObj = costFunc$new("VAR")
  PELTObj$costFunc = costFuncObj #Modify cost funcc
  expect_warning(expect_false(PELTObj$eval(0,2)== 0.5), "singular") #No longer use the original cost function

  PELTObj$minSize = 5L
  expect_true(PELTObj$minSize == 5L)

  PELTObj$jump = 5L
  expect_true(PELTObj$jump == 5L)
})


test_that("Expect error when segment is too short", {

  X1 = matrix(rep(0,10))
  costFuncObj = costFunc$new("L2")
  PELTObj = PELT$new(costFunc = costFuncObj, minSize = 10)
  expect_error(PELTObj$fit(X1)) ## nObs must be at least 2*minSize

})



test_that("Expect error when sizes mismatch", {

  X1 = matrix(rnorm(100))
  costFuncObj = costFunc$new("LinearL2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  PELTObj$tsMat = X1
  PELTObj$covariates = matrix(rnorm(99))
  expect_error(PELTObj$fit()) #Sizes not match
  PELTObj$covariates = matrix(rnorm(100))
  expect_no_error(PELTObj$fit())
})

