#test-PELT.R

test_that("PELT with L2/SIGMA/VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  #L2
  costFuncObj = costFunc$new("L2")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)
  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  #SIGMA
  costFuncObj = costFunc$new("SIGMA")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)
  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  #VAR
  costFuncObj = costFunc$new("VAR")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)

  expect_warning(expect_equal(PeltObj$predict(pen = 0.1), seq(50,150,50)),
                 "Some systems seem singular!")
  expect_warning(expect_false(PeltObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(PeltObj$eval(0,50), 0)), "seems singular")

})
