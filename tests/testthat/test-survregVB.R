test_that("survregVB", {
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

  result <- survregVB(Surv(time, infect) ~ trt + fev, dnase,
                    501, 500, c(4.4, 0.25, 0.04), 1,
                    max_iteration = 100, threshold = 0.0005)
  expected <- list(
    ELBO = -4857.1776,
    alpha = 744,
    omega = 674.93613,
    Sigma = matrix(c(0.03621216, -0.01027572, -0.00047336,
                     -0.01027572, 0.01987219, 0.00002210,
                     -0.00047336, 0.00002210, 0.00000827),
                   nrow = 3, byrow = TRUE,
                   dimnames = list(c("(Intercept)", "trt", "fev"),
                                   c("(Intercept)", "trt", "fev"))),
    mu =  c("(Intercept)" = 4.11256275, "trt" = 0.41556503,
            "fev" = 0.0214483),
    iterations = 10,
    n = 645,
    call = quote(survregVB(formula = Surv(time, infect) ~ trt + fev,
                           data = dnase, alpha_0 = 501,
                           omega_0 = 500, mu_0 = c(4.4, 0.25, 0.04),
                           v_0 = 1, max_iteration = 100,
                           threshold = 5e-04))
  )
  class(expected) <- 'survregVB'
  expect_equal(result, expected)
})

