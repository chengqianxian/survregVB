library(survival)
set.seed(1)
x1 <- rnorm(300, 1, 0.2)
x2 <- rbinom(300, 1, 0.5)
z <- rlogis(300)

# generate survival times
beta0 <- 0.5
beta1 <- 0.2
beta2 <- 0.8
b <- 0.8
Y <- beta0 + beta1 * x1 + beta2 * x2 + b * z
T <- exp(Y)

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
v_0 <- 0.1
alpha_0 <- 11
omega_0 <- 10

test_that("survregVB.fit", {
  result <- survregVB.fit(Surv(T, delta), X, alpha_0, omega_0, mu_0, v_0,
    max_iteration = 100, threshold = 0.01
  )
  expected <- list(
    ELBO = -1723.493128,
    alpha = 311,
    omega = 250.8086135,
    mu = c(0.05739213186, 0.5273365895, 0.8418679603),
    Sigma = matrix(
      c(
        0.18803590234, -0.1735466575414, -0.0134736286120,
        -0.17354665754,  0.1726178469305, -0.0000653350561,
        -0.01347362861, -0.0000653350561,  0.0254160913341
      ),
      nrow = 3, ncol = 3, byrow = FALSE,
      dimnames = list(NULL, NULL)
    ),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)

  result <- survregVB.fit(Surv(T.10, delta.10), X, alpha_0, omega_0,
    mu_0, v_0,
    max_iteration = 100,
    threshold = 0.01
  )
  expected <- list(
    ELBO = -1509.54707,
    alpha = 279,
    omega = 220.2464816,
    mu = c(-0.02694633027, 0.5870221907, 0.8883491399),
    Sigma = matrix(
      c(
        0.18618723857, -0.1720991407945, -0.0136226390896,
        -0.17209914079,  0.1713795066094,  0.0002720367882,
        -0.01362263909,  0.0002720367882,  0.0252914786949
      ),
      nrow = 3, ncol = 3, byrow = FALSE,
      dimnames = list(NULL, NULL)
    ),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)

  result <- survregVB.fit(Surv(T.30, delta.30), X, alpha_0, omega_0,
    mu_0, v_0,
    max_iteration = 100,
    threshold = 0.01
  )
  expected <- list(
    ELBO = -1149.432804,
    alpha = 220,
    omega = 181.5454727,
    mu = c(0.1478980612, 0.4472516245, 0.8735885567),
    Sigma = matrix(
      c(
        0.22110740073, -0.205050440078, -0.016367876916,
        -0.20505044008,  0.203811682106,  0.001569021946,
        -0.01636787692,  0.001569021946,  0.029198756794
      ),
      nrow = 3, ncol = 3, byrow = FALSE,
      dimnames = list(NULL, NULL)
    ),
    iterations = 8,
    n = 300
  )
  expect_equal(result, expected)
})
