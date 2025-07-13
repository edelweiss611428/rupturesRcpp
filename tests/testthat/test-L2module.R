#test-L2module.R
set.seed(12345)
R_L2eval = function(X, start, end){
  R_start = start+1
  R_end = end
  nc = ncol(X)

  Xe = matrix(X[R_start:R_end,], ncol = nc)
  cMXe = colMeans(Xe)

  return(sum(sweep(Xe, 2, cMXe, FUN = "-")^2))
}

set.seed(1)
tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)))

nr = nrow(tsMat)
tsMat_L2module = new(rupturesRcpp:::Cost_L2, tsMat)

binSegObj = binSeg$new() #L2 by default
binSegObj$fit(tsMat)

PELTObj = PELT$new() #L2 by default
PELTObj$fit(tsMat)

WinObj = Window$new() #L2 by default
WinObj$fit(tsMat)

test_that("Expect equal", {

  nCases = 10
  idx1 = sample.int(nr-2, nCases) #to avoid case nr-1 sample(100) wont give you 100
  idx2 = integer(nCases)

  for(i in 1:10){
    idx2[i] = sample((idx1[i]+1):nr, 1)
  }

  for(i in 1:10){

    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(tsMat_L2module$eval(idx1[i],idx2[i]), gs)
    expect_equal(binSegObj$eval(idx1[i],idx2[i]), gs)
    expect_equal(PELTObj$eval(idx1[i],idx2[i]), gs)
    expect_equal(WinObj$eval(idx1[i],idx2[i]), gs)
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
  expect_error(WinObj$eval(0,0),
               regexp = "smaller")

})



test_that("Cpp L2 module gives 0 if start>=end+1", {

  expect_equal(tsMat_L2module$eval(0,0), 0)
  expect_equal(tsMat_L2module$eval(0,1), 0)

})
