# calling context() avoids a strange bug in testthat 2.0.0
# cannot actually run brms models in tests as it takes way too long
context("Tests for bmscstan error messages")

test_that("bmsc_stan produces expected errors", {
  dat.ctrl <- data.frame(y = rnorm(100), x = rnorm(100), id = rep(1:5, 20), cond = rep(c(1,2), each = 50))
  dat.sc   <- data.frame(y = rnorm(20), x = rnorm(20), cond = rep(1:2, 20))

  # missing arguments
  expect_error(BMSC(),
               "the argument \"formula\" is not specified",
               fixed = TRUE)
  expect_error(BMSC(y~x*cond+(cond|id), data_ctrl = dat.ctrl),
               "the dataframe \"data_sc\" is not specified")
  expect_error(BMSC(y~x*cond+(cond|id), data_sc = dat.sc),
               "the dataframe \"data_ctrl\" is not specified")
  expect_error(BMSC(y~x*cond+(cond|id), data_ctrl = dat.ctrl, data_sc = dat.sc, typeprior = "gamma"),
               "Not a valid typeprior")
})
