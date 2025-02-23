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
for (l in 1:75){
  Y[l] <- beta0 + beta1 * x1[l] + beta2 * x2[l] + random.int[(l-1) %/% 5 + 1] + b * epsilon[l]
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
  result <- survregVB.frailty.fit(Y = Surv(T, delta), X = X, alpha_0 = alpha_0,
                                  omega_0 = omega_0, mu_0 = mu_0, v_0 = v_0,
                                  lambda_0 = lambda_0, eta_0 = eta_0,
                                  cluster = cluster, max_iteration = 100,
                                  threshold = 0.01)
  expected <- list(
    clustered = TRUE,
    ELBO = -312.412483,
    alpha = 78,
    omega = 52.37919,
    mu = c(-0.292427124, 0.808340621, 0.557493622),
    Sigma = matrix(c(0.50246823, -0.45647333, -0.03885537,
                     -0.45647333, 0.44525833, 0.00434697,
                     -0.03885537, 0.00434697, 0.06289714),
                   nrow = 3, byrow = TRUE,
                   dimnames = list(NULL, NULL)),
    tau = setNames(c(-0.207712506, 0.703780589, 0.434283093, 0.127487655, 1.612282757,
                      0.285318913, 0.707166236, 1.832609362, -1.121440475, -0.763991484,
                      -0.522234151, -0.675163493, -1.585781794, -0.792693127, 0.421673759),
                    as.character(1:15)),
    sigma = setNames(c(0.224621934, 0.162451332, 0.162451332, 0.188527839, 0.188527839,
                        0.188527839, 0.224621934, 0.188527839, 0.188559766, 0.162451332,
                        0.162451332, 0.277809145, 0.162451332, 0.188527839, 0.233526067),
                      as.character(1:15)),
    lambda = 10.5,
    eta = 10.0387529,
    iterations = 13,
    n = 75
  )

  expect_equal(result, expected)
})
