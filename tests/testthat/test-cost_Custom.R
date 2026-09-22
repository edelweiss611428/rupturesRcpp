#Custom (user-defined, R-callback) cost module

# ========================================================
#           (R) User-defined cost/param functions
# ========================================================

# Replicates the built-in L2 cost exactly, but as a plain R closure operating on
# the raw segment matrix -- used to check that "Custom" matches "L2" bit-for-bit.
# `evalFun`/`paramFun` are always called as f(segment, a, b); a, b are unused here.
R_customL2eval = function(segment, a, b){
  segment = as.matrix(segment)
  cm = colMeans(segment)
  sum(sweep(segment, 2, cm, FUN = "-")^2)
}

R_customL2params = function(segment, a, b){
  colMeans(as.matrix(segment))
}

# ========================================================
#                   Simulated datasets
# ========================================================

set.seed(12345)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)))
nr = nrow(tsMat)

nCases = 10
idx1 = sample.int(nr-2, nCases)
idx2 = integer(nCases)

for(i in 1:nCases){
  idx2[i] = sample((idx1[i]+1):nr, 1)
}


test_that("Expect C++ .eval()/.get_params() methods in Cost_RFunc module give the correct results", {

  RFuncModule = new(rupturesRcpp:::Cost_RFunc, tsMat, R_customL2eval, R_customL2params, TRUE)
  L2module = new(rupturesRcpp:::Cost_L2, tsMat, TRUE)

  for(i in 1:nCases){
    expect_equal(RFuncModule$eval(idx1[i],idx2[i]), L2module$eval(idx1[i],idx2[i]))
  }

  expect_equal(RFuncModule$get_params(0,nr)$params, L2module$get_params(0,nr)$mean)

  expect_error(RFuncModule$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  expect_error(RFuncModule$eval(-1,nr),
               regexp = "out of bounds") #arma::mat indexing error

  expect_no_error(RFuncModule$resetWarning(FALSE)) #Generally does nothing here

})

test_that("Cost_RFunc without a paramFun returns an empty list from get_params()", {

  RFuncModule = new(rupturesRcpp:::Cost_RFunc, tsMat, R_customL2eval, NULL, TRUE)
  expect_equal(RFuncModule$get_params(0,nr), list())

})

test_that("Cost_RFunc passes (segment, a, b) through, using the same (a,b] convention as $eval()", {

  seenArgs = list()
  recordingEval = function(segment, a, b){
    seenArgs[[length(seenArgs)+1]] <<- list(nrow = nrow(as.matrix(segment)), a = a, b = b)
    0
  }
  RFuncModule = new(rupturesRcpp:::Cost_RFunc, tsMat, recordingEval, NULL, TRUE)

  RFuncModule$eval(10L, 25L)
  expect_equal(seenArgs[[1]], list(nrow = 15L, a = 10L, b = 25L))

})

test_that("Expect $eval()/$predict() for costFunc = 'Custom' matches costFunc = 'L2' exactly, for PELT/binSeg/Window", {

  customCF = costFunc$new("Custom", evalFun = R_customL2eval, paramFun = R_customL2params)
  L2CF = costFunc$new("L2")

  for(Algo in list(PELT, binSeg, Window)){

    customObj = Algo$new(costFunc = customCF)
    L2Obj = Algo$new(costFunc = L2CF)
    customObj$fit(tsMat)
    L2Obj$fit(tsMat)

    for(i in 1:nCases){
      expect_equal(customObj$eval(idx1[i],idx2[i]), L2Obj$eval(idx1[i],idx2[i]))
    }

    expect_identical(customObj$predict(pen = 25), L2Obj$predict(pen = 25))

  }

})

test_that("costFunc = 'Custom' validates evalFun/paramFun and requires evalFun to be set", {

  expect_error(costFunc$new("Custom")$pass(),
               regexp = "requires `evalFun` to be set")

  expect_error(costFunc$new("Custom", evalFun = "not a function"),
               regexp = "`evalFun` must be a function!")

  expect_error(costFunc$new("Custom", evalFun = R_customL2eval, paramFun = "not a function"),
               regexp = "`paramFun` must be a function or NULL!")

  expect_no_error(costFunc$new("Custom", evalFun = R_customL2eval))

  customCF = costFunc$new("Custom", evalFun = R_customL2eval)
  expect_null(customCF$paramFun)
  expect_equal(customCF$pass(), list(costFunc = "Custom", evalFun = R_customL2eval, paramFun = NULL))

  expect_error(PELT$new(costFunc = costFunc$new("Custom"))$fit(tsMat),
               regexp = "requires `evalFun` to be set")

})

test_that("$get_params() for costFunc = 'Custom' uses paramFun when supplied, else an empty list", {

  withParams = PELT$new(costFunc = costFunc$new("Custom", evalFun = R_customL2eval, paramFun = R_customL2params))
  withParams$fit(tsMat)
  expect_equal(withParams$.__enclos_env__$private$.PELTModule$get_params(0L, nr)$params,
               R_customL2params(tsMat))

  noParams = PELT$new(costFunc = costFunc$new("Custom", evalFun = R_customL2eval))
  noParams$fit(tsMat)
  expect_equal(noParams$.__enclos_env__$private$.PELTModule$get_params(0L, nr), list())

})

test_that("describe() reports evalFun/paramFun for costFunc = 'Custom'", {

  PELTObj = PELT$new(costFunc = costFunc$new("Custom", evalFun = R_customL2eval, paramFun = R_customL2params))
  PELTObj$fit(tsMat)

  desc = PELTObj$describe()
  expect_identical(desc$evalFun, R_customL2eval)
  expect_identical(desc$paramFun, R_customL2params)
  expect_no_error(PELTObj$describe(TRUE))

})

test_that("A user-defined cost can close over external, position-aligned data via (a,b]", {

  # `externalSeries` is never passed to `PELT`/`tsMat` at all -- `evalFun` reaches it purely by
  # lexical closure, aligning it to the current segment using the `a`/`b` bounds Cost_RFunc passes
  # through. This is the main capability plain segment-only costs can't offer.
  set.seed(42)
  externalSeries = as.matrix(rnorm(nr))

  # Cost: RSS of regressing `tsMat`'s segment on the *aligned* slice of `externalSeries`.
  externalRegCost = function(segment, a, b){
    x = externalSeries[(a+1):b, , drop = FALSE]
    sum(lm(segment ~ x)$residuals^2)
  }

  customCF = costFunc$new("Custom", evalFun = externalRegCost)
  customObj = PELT$new(costFunc = customCF)
  customObj$fit(tsMat)

  # Cross-check against the built-in LinearL2 cost fitted with `externalSeries` as `covariates`
  # (which *is* told about the alignment directly) -- the two should agree, on segments large
  # enough (n >> p = 2 here) that the OLS fit is well-determined either way it's computed.
  linCF = costFunc$new("LinearL2", intercept = TRUE)
  linObj = PELT$new(costFunc = linCF)
  linObj$fit(tsMat, externalSeries)

  for(bounds in list(c(0,50), c(50,150), c(150,200), c(0,200))){
    expect_equal(customObj$eval(bounds[1], bounds[2]), linObj$eval(bounds[1], bounds[2]))
  }

})
