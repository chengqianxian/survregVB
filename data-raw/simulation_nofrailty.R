## code to prepare `data1` dataset goes here

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
Time <- exp(Y)

# generate censoring times
set.seed(1)
cen.time.10 <- runif(300, 0, 48)
cen.time.30 <- runif(300, 0, 17)

# obtain observed time
Time.10 <- pmin(Time, cen.time.10)
Time.30 <- pmin(Time, cen.time.30)

# obtain censoring indicator
delta <- rep(1, 300)
delta.10 <- ifelse(Time == Time.10, 1, 0)
delta.30 <- ifelse(Time == Time.30, 1, 0)

simulation_nofrailty <- data.frame(x1, x2, Time, Time.10, Time.30, delta,
                                   delta.10, delta.30)

usethis::use_data(simulation_nofrailty, overwrite = TRUE)
