library(survival)
set.seed(1)
x1 <- rnorm(300, 1, 0.2)
x2 <- rbinom(300, 1, 0.5)
z <- rlogis(300)
beta0 <- 0.5
beta1 <- 0.2
beta2 <- 0.8
b <- 0.8
y <- beta0 + beta1 * x1 + beta2 * x2 + b * z
T <- exp(y)

# generate censoring times
set.seed(1)
cen.time.10 <- runif(300, 0, 48)
cen.time.30 <- runif(300, 0, 17)

# obtain observed time
T.10 <- pmin(T, cen.time.10)
T.30 <- pmin(T, cen.time.30)

# obtain censoring indicator
delta <- rep(1, 300)
delta.10 <- ifelse(T == T.10, 1, 0)
delta.30 <- ifelse(T == T.30, 1, 0)

# create X matrix
X <- matrix(c(rep(1, 300), x1, x2), nrow = 300)

# priors, use non-informative priors
mu_0 <- c(0, 0, 0)
v_0 <-  0.1
alpha_0 <- 11
omega_0 <- 10

test_that("survregVB.fit", {
  result <- survregVB.fit(Surv(T, delta), X, alpha_0, omega_0, mu_0, v_0,
                         max_iteration = 100, threshold = 0.01)
  expected <- list(
    ELBO = -1723.49313,
    alpha = 311,
    omega = 250.808607,
    Sigma = matrix(c(0.18803590, -0.17354666, -0.01347363,
                     -0.17354666, 0.17261785, -0.00006534,
                     -0.01347363, -0.00006534, 0.02541609),
                   nrow = 3, ncol = 3, byrow = FALSE,
                   dimnames = list(NULL, NULL)),
    mu = c(0.05739117, 0.52733610, 0.84186673),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)

  result <- survregVB.fit(Surv(T.10, delta.10), X, alpha_0, omega_0,
                            mu_0, v_0, max_iteration = 100,
                            threshold = 0.01)
  expected <- list(
    ELBO = -1509.54706,
    alpha = 279,
    omega = 220.24649,
    Sigma = matrix(c(
      0.18618723, -0.17209913, -0.01362264,
      -0.17209913, 0.17137950, 0.00027204,
      -0.01362264, 0.00027204, 0.02529148),
      nrow = 3, ncol = 3, byrow = FALSE,
      dimnames = list(NULL, NULL)),
    mu = c(-0.02694604, 0.58702326, 0.88834972),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)

  result <- survregVB.fit(Surv(T.30, delta.30), X, alpha_0, omega_0,
                          mu_0, v_0, max_iteration = 100,
                          threshold = 0.01)
  expected <- list(
    ELBO = -1149.4328,
    alpha = 220,
    omega = 181.545475,
    Sigma = matrix(c(
      0.22110740, -0.20505044, -0.01636788,
      -0.20505044, 0.20381168, 0.00156902,
      -0.01636788, 0.00156902, 0.02919876
    ), nrow = 3, ncol = 3, byrow = FALSE,
    dimnames = list(NULL, NULL)),
    mu = c(0.14789765, 0.44725114, 0.87358818),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)
})
