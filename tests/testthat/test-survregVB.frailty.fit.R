library(survival)
cluster <- rep(1:15, each = 5)

set.seed(1)
x1 <- rnorm(75, 1, 0.2)
x2 <- rbinom(75, 1, 0.5)
epsilon <- rlogis(75)
cen.time.15 <- runif(75, 0, 48)

beta0 <- 0.5
beta1 <- 0.2
beta2 <- 0.8
b <- 0.8

Y <- numeric(50)
random.int <- rnorm(15, 0, 1)
for (l in 1:75) {
  Y[l] <- beta0 + beta1 * x1[l] + beta2 * x2[l] + random.int[(l - 1) %/% 5 + 1] + b * epsilon[l]
}
T <- exp(Y)

# obtain observed time
T.15 <- pmin(T, cen.time.15)

# obtain censoring indicator
delta <- rep(1, 75)
delta.15 <- ifelse(T == T.15, 1, 0)

# prepare for the data structure
X <- matrix(0, nrow = 75, ncol = 3)
X[, 1] <- 1
X[, 2] <- x1
X[, 3] <- x2

# priors, informative priors
mu_0 <- c(0, 0, 0)
v_0 <- 0.1
alpha_0 <- 3
omega_0 <- 2
lambda_0 <- 3
eta_0 <- 2

test_that("survregVB.frailty.fit", {
  result <- survregVB.frailty.fit(
    Y = Surv(T, delta), X = X, alpha_0 = alpha_0,
    omega_0 = omega_0, mu_0 = mu_0, v_0 = v_0,
    lambda_0 = lambda_0, eta_0 = eta_0,
    cluster = cluster, max_iteration = 100,
    threshold = 0.01
  )
  expected <- list(
    clustered = TRUE,
    ELBO = -312.4124906,
    alpha = 78,
    omega = 52.37919667,
    mu = c(-0.29243491, 0.80833985, 0.55749383),
    Sigma = matrix(
      c(
        0.50246837, -0.456473457, -0.038855380,
        -0.45647346,  0.445258447,  0.004346971,
        -0.03885538,  0.004346971,  0.062897163
      ),
      nrow = 3, byrow = TRUE,
      dimnames = list(NULL, NULL)
    ),
    tau = setNames(
      c(
        -0.2077060470, 0.7037876292, 0.4342901316, 0.1274945678, 1.6122896798,
        0.2853256227, 0.7071727654, 1.8326163036, -1.1214338058, -0.7639845143,
        -0.5222271141, -0.6751576177, -1.5857748947, -0.7926863618, 0.4216802745
      ),
      as.character(1:15)
    ),
    sigma = setNames(
      c(
        0.2246220256, 0.1624513933, 0.1624513933, 0.1885279119, 0.1885279119,
        0.1885279119, 0.2246220256, 0.1885279119, 0.1885598395, 0.1624513933,
        0.1624513933, 0.2778092651, 0.1624513933, 0.1885279119, 0.2335261631
        ),
      as.character(1:15)
    ),
    lambda = 10.5,
    eta = 10.0387575,
    iterations = 13,
    n = 75
  )

  expect_equal(result, expected)
})
