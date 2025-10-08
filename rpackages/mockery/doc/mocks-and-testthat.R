## ----include=FALSE------------------------------------------------------------
library(mockery)

library(knitr)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
summary <- NULL

# pretend we're testing a package
Sys.setenv(TESTTHAT = "true")
Sys.setenv(TESTTHAT_PKG = "mockery")

## ----create_mock--------------------------------------------------------------
m <- mock()

## ----return_values------------------------------------------------------------
m <- mock(1, 2, 3)
m()
m()
m()

## ----return_expression--------------------------------------------------------
x <- 1
y <- 2
m <- mock(x + y)
m()

## ----cycle_no-----------------------------------------------------------------
try({
m <- mock(1, 2)
m()
m()
m()
})

## ----cycle_true---------------------------------------------------------------
m <- mock(1, 2, cycle = TRUE)
m()
m()
m()
m()

## ----cycle_expression---------------------------------------------------------
x <- 1
y <- 2
m <- mock(1, x + y, cycle = TRUE)

m()
m()

## ----cycle_expression_2nd-----------------------------------------------------
y <- 10
m()
m()

## ----return_expression_env----------------------------------------------------
x <- 1
y <- 2
e <- new.env()
m <- mock(x + y, envir = e, cycle = TRUE)

m()
e$x <- 10
m()

## ----with_mock, message=FALSE-------------------------------------------------
library(testthat)

m <- mock(1)
f <- function (x) summary(x)
with_mocked_bindings(f = m, {
  expect_equal(f(iris), 1)
})

## ----expect_called------------------------------------------------------------
m <- mock(1, 2)

m()
expect_called(m, 1)

m()
expect_called(m, 2)

## ----expect_called_error------------------------------------------------------
try({
expect_called(m, 1)
expect_called(m, 3)
})

## ----expect_call--------------------------------------------------------------
m <- mock(1)
with_mocked_bindings(summary = m, {
  summary(iris)
})

expect_call(m, 1, summary(iris))

## ----call_doesnt_match--------------------------------------------------------
try({
expect_call(m, 1, summary(x))
#> Error: expected call summary(x) does not mach actual call summary(iris).
})

## ----expect_args--------------------------------------------------------------
expect_args(m, 1, iris)

## ----expect_args_different----------------------------------------------------
try({
expect_args(m, 1, iris[-1, ])
})

## ----expect_args_named--------------------------------------------------------
try({
m <- mock(1)
with_mocked_bindings(summary = m, {
  summary(object = iris)
})

expect_args(m, 1, object = iris)
})

## ----expect_args_unnamed------------------------------------------------------
try({
expect_args(m, 1, iris)
})

