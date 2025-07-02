library("testthat")
set.seed(1)
n = 200
nRegimes = 5
ARcor = rbeta(nRegimes, 5, 1)*sample(c(-1,1),nRegimes,T)
X = numeric(0)
idx = 0

for(i in 1:nRegimes){
  X = c(X, filter(rnorm(n), filter = ARcor[i], method = "recursive"))
}
tsMat = as.matrix(X)
VAR5Obj = new(rupturesRcpp::Cost_VAR, tsMat, 5) #5th order VAR
J = 5*1+1

test_that("VAR gives 0 when (start, end] is smaller than J", {
  for(i in 1:(J-1)){
    expect_equal(VAR5Obj$effEvalCpp(0,i,TRUE, 10^-6), 0)
  }
})

