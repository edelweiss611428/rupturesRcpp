#test-PELT.R

test_that("PELT with L2/SIGMA/VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  costFunc = c("L2", "SIGMA", "VAR")
  for(i in 1:3){

    costFuncObj = createCostFunc(costFunc[i])
    PeltObj = PELT$new(costFuncObj = costFuncObj)
    PeltObj$fit(tsMat)
    expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  }


})

