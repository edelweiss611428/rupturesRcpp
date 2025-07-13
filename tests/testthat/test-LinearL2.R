#test-L2module.R
set.seed(12345)
R_LinearL2eval = function(Y, X, start, end){

  subX = X[(start+1):end,]
  subY = Y[(start+1):end,]

  return(sum(lm(subY~subX)$residuals^2))
}


X = cbind(rnorm(100), rnorm(100))
Y = X*cbind(rep(c(1,5), each = 50), rep(c(1,5), each = 50)) + cbind(rnorm(100),rnorm(100))

nr = 100

tsMat_LinearL2module = new(rupturesRcpp:::Cost_LinearL2, Y, X, TRUE, TRUE) #include intercept

WinObj = Window$new(costFunc = costFunc$new("LinearL2")) #L2 by default
WinObj$fit(Y,X)


test_that("Expect equal", {

  nCases = 10
  idx1 = sample.int(nr-2, nCases) #to avoid case nr-1 sample(100) wont give you 100
  idx2 = integer(nCases)

  for(i in 1:10){
    idx2[i] = sample((idx1[i]+1):nr, 1)
  }


  for(i in 1:10){

    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(tsMat_LinearL2module$eval(idx1[i],idx2[i]), gs)
    expect_equal(WinObj$eval(idx1[i],idx2[i]), gs)
  }

})

