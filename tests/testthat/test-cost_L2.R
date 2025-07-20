#L2 cost module

# ========================================================
#                    (R) L2 cost function
# ========================================================

R_L2eval = function(X, start, end){
  R_start = start+1
  R_end = end
  nc = ncol(X)

  Xe = matrix(X[R_start:R_end,], ncol = nc)
  cMXe = colMeans(Xe)

  return(sum(sweep(Xe, 2, cMXe, FUN = "-")^2))
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


test_that("Expect C++ .eval() method in L2 cost module gives the correct results", {

  L2module = new(rupturesRcpp:::Cost_L2, tsMat, TRUE) # (warnOnce) = (TRUE) by default

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(L2module$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(L2module$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  expect_error(L2module$eval(-1,nr),
               regexp = "out of bounds") #arma::mat indexing error

  #Give 0 if end - start = 0 or 1
  expect_equal(L2module$eval(0,0), 0)
  expect_equal(L2module$eval(0,1), 0)

  expect_no_error(L2module$resetWarning(FALSE)) #Generally does nothing here


})

test_that("Expect $eval() method in PELT_L2 gives the correct results/error message", {


  PELTObj = PELT$new() #This uses L2 cost module by default
  PELTObj$fit(tsMat)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(PELTObj$eval(idx1[i],idx2[i]), gs)
  }


  expect_error(PELTObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(PELTObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(PELTObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(PELTObj$eval(0,1), 0)


})



test_that("Expect .eval() method in C++ PELT_L2 class gives the correct results/error message", {

  #.constructor<arma::mat, int, int>()
  PELTCppObj = new(PELTCpp_L2, tsMat, 1L, 1L)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(PELTCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(PELTCppObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(PELTCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(PELTCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(PELTCppObj$eval(0,1), 0)


})

test_that("Expect $eval() method in binSeg_L2 gives the correct results", {

  binSegObj = binSeg$new() #This uses L2 cost module by default
  binSegObj$fit(tsMat)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(binSegObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(binSegObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(binSegObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(binSegObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(binSegObj$eval(0,1), 0)

})


test_that("Expect .eval() method in C++ binSeg_L2 class gives the correct results/error message", {

  #.constructor<arma::mat, int, int>()
  binSegCppObj = new(binSegCpp_L2, tsMat, 1L, 1L)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(binSegCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(binSegCppObj$eval(-1,nr),
               regexp = "0 <= start < nSamples` must be true!")

  expect_error(binSegCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(binSegCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(binSegCppObj$eval(0,1), 0)


})



test_that("Expect $eval() method in Window_L2 gives the correct results", {


  WinObj = Window$new() #This uses L2 cost module by default
  WinObj$fit(tsMat)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(WinObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(WinObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(WinObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(WinObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(WinObj$eval(0,1), 0)

})



test_that("Expect .eval() method in C++ Window_L2 class gives the correct results/error message", {

  #.constructor<arma::mat, int, int,int>()
  WindowCppObj = new(windowCpp_L2, tsMat, 1L, 1L, 10)

  for(i in 1:nCases){
    gs = R_L2eval(tsMat, idx1[i],idx2[i])
    expect_equal(WindowCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(WindowCppObj$eval(-1,nr),
               regexp = "0 <= start < nSamples` must be true!")

  expect_error(WindowCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(WindowCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(WindowCppObj$eval(0,1), 0)

})


test_that("C++ L2 module gives 0 if start>=end+1", {

  L2module = new(rupturesRcpp:::Cost_L2, tsMat, TRUE)
  expect_equal(L2module$eval(0,0), 0)
  expect_equal(L2module$eval(0,1), 0)

})
