## code to prepare `DATASET` dataset goes here
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
  Y[l] <- beta0 + beta1 * x1[l] + beta2 * x2[l] +
    random.int[(l - 1) %/% 5 + 1] + b * epsilon[l]
}
T <- exp(Y)

# obtain observed time
T.15 <- pmin(T, cen.time.15)

# obtain censoring indicator
delta <- rep(1, 75)
delta.15 <- ifelse(T == T.15, 1, 0)

simulation_frailty <- data.frame(x1, x2, T, T.15, delta, delta.15, cluster)

usethis::use_data(simulation_frailty, overwrite = TRUE)
