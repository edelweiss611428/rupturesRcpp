if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

if (!requireNamespace("changepoint", quietly = TRUE)) {
  install.packages("changepoint")
}

library("changepoint")
library("reticulate")
library("testthat")

virtualenv_create("r-reticulate")  # or use conda_create("r-reticulate") if you prefer conda
virtualenv_install("r-reticulate", packages = "ruptures")

test_that("ruptures has been installed" ,{
  expect_true(reticulate::py_module_available("ruptures"))
})

set.seed(12345)


library(testthat)
library(reticulate)

# Import ruptures python module
ruptures = import("ruptures")
np = import("numpy")

test_that("rupturesRcpp L2 cost matches ruptures python CostL2", {

  # Generate test data

  signal = c(rnorm(100,0.25), rnorm(100,-0.25))

  # Convert signal to numpy array for ruptures

  pySignal = np$array(signal)
  Rsignal = matrix(signal)

  #cost modules
  pyL2 = ruptures$costs$CostL2()$fit(pySignal)
  rL2 = new(rupturesRcpp:::Cost_L2, Rsignal)

  # Example segments
  nCases = 10
  nr = length(signal)
  idx1 = sample.int(nr-2, nCases) #to avoid case nr-1 sample(100) wont give you 100
  idx2 = integer(nCases)


  for(i in 1:10){
    idx2[i] = sample((idx1[i]+1):nr, 1)
  }

  for (i in 1:10) {
    start = idx1[i]
    end = idx2[i]

    py_eval = pyL2$error(as.integer(start), as.integer(end))
    r_eval = rL2$eval(start,end)

    expect_equal(py_eval, r_eval)
  }
})


test_that("Python ruptures PELT-L2 matches rupturesRcpp result (minSize = jump = 1)", {

  use_virtualenv("r-reticulate", required = TRUE)

  ruptures = import("ruptures")
  np = import("numpy")

  # Simulated signal
  signal = c(rnorm(100,0.25), rnorm(100,-0.25))
  pySignal = np$array(signal)

  Pyalgo = ruptures$Pelt(model = "l2", min_size = 1L, jump = 1L)$fit(pySignal)

  Rsignal = matrix(signal)
  Ralgo = rupturesRcpp::PELT$new() #The same configurations by default
  Ralgo$fit(Rsignal)


  for(i in 10^(seq(-1,1, by = 0.25))){
    Pyresult = Pyalgo$predict(pen = i)
    Rresult = Ralgo$predict(pen = i)

    #changepoint package
    CPalgo = cpt.mean(signal, method = "PELT", penalty = "Manual", pen.value = i, minseglen = 1)
    CPresult = c(cpts(CPalgo), 200)

    # Expectation
    expect_equal(Pyresult, Rresult)
    # expect_equal(CPresult, Rresult)
  }


})

test_that("Python ruptures PELT-L2 matches rupturesRcpp result (minSize 5, = jump = 1)", {

  use_virtualenv("r-reticulate", required = TRUE)

  ruptures = import("ruptures")
  np = import("numpy")

  # Simulated signal
  signal = c(rnorm(100,0.25), rnorm(100,-0.25))
  pySignal = np$array(signal)

  Pyalgo = ruptures$Pelt(model = "l2", min_size = 5L, jump = 1L)$fit(pySignal)

  Rsignal = matrix(signal)
  Ralgo = rupturesRcpp::PELT$new(minSize = 5L) #The same configurations by default
  Ralgo$fit(Rsignal)


  for(i in 10^(seq(-1,1, by = 0.25))){
    Pyresult = Pyalgo$predict(pen = i)
    Rresult = Ralgo$predict(pen = i)

    #changepoint package
    CPalgo = cpt.mean(signal, method = "PELT", penalty = "Manual", pen.value = i, minseglen = 5)
    CPresult = c(cpts(CPalgo), 200)

    print(Rresult)
    print(CPresult)
    # Expectation
    expect_equal(Pyresult, Rresult)
    expect_equal(CPresult, Rresult)
  }


})



test_that("Python ruptures PELT-L2 matches rupturesRcpp result (minSize = jump = 5)", {

  use_virtualenv("r-reticulate", required = TRUE)

  ruptures = import("ruptures")
  np = import("numpy")

  # Simulated signal
  signal = c(rnorm(100,0.25), rnorm(100,-0.25))
  pySignal = np$array(signal)

  Pyalgo = ruptures$Pelt(model = "l2", min_size = 5L, jump = 5L)$fit(pySignal)

  Rsignal = matrix(signal)
  Ralgo = rupturesRcpp::PELT$new(minSize = 5L, jump = 5L) #L2 by default
  Ralgo$fit(Rsignal)


  for(i in 10^(seq(-1,1, by = 0.25))){
    Pyresult = Pyalgo$predict(pen = i)
    Rresult = Ralgo$predict(pen = i)

    expect_equal(Pyresult, Rresult)

  }
})

