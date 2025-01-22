library(survival)
first <- subset(rhDNase, !duplicated(id)) # first row for each subject
dnase <- tmerge(first, first, id=id, tstop=as.numeric(end.dt -entry.dt))

# Subjects whose fu ended during the 6 day window are the reason for
#  this next line
temp.end <- with(rhDNase, pmin(ivstop+6, end.dt-entry.dt))
dnase <- tmerge(dnase, rhDNase, id=id,
                infect=event(ivstart),
                end=  event(temp.end))
# toss out the non-at-risk intervals, and extra variables
#  3 subjects had an event on their last day of fu, infect=1 and end=1
dnase <- subset(dnase, (infect==1 | end==0), c(id:trt, fev:infect))
dnase <- subset(dnase, !duplicated(id))
dnase$time <- dnase$tstop - dnase$tstart

# data preparation
y <- dnase$time
n <- nrow(dnase)
X <- matrix(c(rep(1, n), dnase$trt, dnase$fev), nrow = n)
delta <- dnase$infect

# priors
mu_0 <- c(4.4, 0.25, 0.04)
v_0 <- 1
alpha_0 <- 501
omega_0 <- 500

expect_list <- list(
  ELBO = c(-4857.1776),
  alpha = 744,
  omega = 674.93613,
  Sigma = matrix(c(0.03621216, -0.01027572, -0.00047336,
                   -0.01027572,  0.01987219,  0.00002210,
                   -0.00047336,  0.00002210,  0.00000827),
                 nrow = 3, ncol = 3, byrow = TRUE),
  mu = matrix(c(4.11256275, 0.41556503, 0.0214483),
              nrow = 1, ncol = 3, byrow = TRUE),
  iterations = 10
)

test_that("survregVB", {
  list <- survregVB(Surv(time, infect) ~ trt + fev, dnase,
                    alpha_0, omega_0, mu_0, v_0,
                    max_iteration = 100, threshold = 0.0005)

  expect_equal(unname(list$ELBO), expect_list$ELBO)
  expect_equal(list$alpha, expect_list$alpha)
  expect_equal(unname(list$omega), expect_list$omega)
  expect_equal(unname(list$Sigma), expect_list$Sigma)
  expect_equal(unname(list$mu), expect_list$mu)
  expect_equal(list$iterations, expect_list$iterations)
})

