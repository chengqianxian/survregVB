test_that("survregVB", {
  library(survival)
  result <- survregVB(Surv(time, infect) ~ trt + fev, dnase,
    501, 500, c(4.4, 0.25, 0.04), 1,
    max_iteration = 100, threshold = 0.0005
  )
  expected <- list(
    ELBO = -4857.09446,
    alpha = 744,
    omega = 674.6480035,
    mu = c(
      "(Intercept)" = 4.112387227, "trt" = 0.415469071,
      "fev" = 0.02129090543
    ),
    Sigma = matrix(
      c(
        0.0362042587960, -1.027351580e-02, -4.732594311e-04,
        -0.0102735158000,  1.986778912e-02,  2.209508562e-05,
        -0.0004732594311,  2.209508562e-05,  8.265796487e-06
      ),
      nrow = 3, byrow = TRUE,
      dimnames = list(
        c("(Intercept)", "trt", "fev"),
        c("(Intercept)", "trt", "fev")
      )
    ),
    iterations = 9,
    n = 645,
    call = quote(survregVB(
      formula = Surv(time, infect) ~ trt + fev,
      data = dnase, alpha_0 = 501,
      omega_0 = 500, mu_0 = c(4.4, 0.25, 0.04),
      v_0 = 1, max_iteration = 100,
      threshold = 5e-04
    ))
  )
  class(expected) <- "survregVB"
  expect_equal(result, expected)
})
