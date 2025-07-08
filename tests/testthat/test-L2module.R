#test-L2module.R

R_L2eval = function(X, start, end){
  R_start = start+1
  R_end = end
  nc = ncol(X)

  Xe = matrix(X[R_start:R_end,], ncol = nc)
  cMXe = colMeans(Xe)

  return(sum(sweep(Xe, 2, cMXe, FUN = "-")^2))
}

set.seed(1)
tsMat = cbind(c(rnorm(10,0), rnorm(10,5,5)))

nr = nrow(tsMat)
tsMat_L2module = new(rupturesRcpp:::Cost_L2, tsMat)

binSegObj = binSeg$new() #L2 by default
binSegObj$fit(tsMat)

PELTObj = PELT$new() #L2 by default
PELTObj$fit(tsMat)

test_that("Expect equal", {

  for(i in 0:(nr-1)){
    for(j in (i+1):nr){
      gs = R_L2eval(tsMat, i, j)
      expect_equal(tsMat_L2module$eval(i,j), gs)
      expect_equal(binSegObj$eval(i,j), gs)
      expect_equal(PELTObj$eval(i,j), gs)
    }
  }

})

test_that("Expect error", {

  expect_error(tsMat_L2module$eval(0,nr+1),
               regexp = "out of bounds")
  expect_error(tsMat_L2module$eval(-1,nr),
               regexp = "out of bounds")
  expect_error(binSegObj$eval(0,0),
               regexp = "smaller") #error if start >= end
  expect_error(PELTObj$eval(0,0),
               regexp = "smaller")

})



test_that("Cpp L2 module gives 0 if start>=end+1", {

  expect_equal(tsMat_L2module$eval(0,0), 0)
  expect_equal(tsMat_L2module$eval(0,1), 0)

})
