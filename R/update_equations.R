library(matlib) # to get inverse matrix
get_alpha_star <- function(alpha_0, delta) return(alpha_0 + sum(delta))

get_omega_star <- function(curr_b, omega_0, delta, Y, X, curr_mu){
  res <- 0
  for (i in 1:nrow(X)) {
    C_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -5, 0,
                  ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -1.701 &
                           (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -5, 0.0426,
                         ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 0 &
                                  (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -1.701, 0.3052,
                                ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 1.702 &
                                         (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > 0, 0.6950,
                                       ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 5 &
                                                (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > 1.702, 0.9574, 1)))))
    res_i <- (delta[i] - C_i * (1 + delta[i])) * (Y[i] - sum(X[i, ] * curr_mu))
    res <- res + res_i
  }
  omega <- omega_0 - res
  return(omega)
}

get_Sigma_star <- function(curr_b, curr_mu, alpha, omega, delta, X, v_0, Y){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  p <- ncol(X)
  X_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in 1:nrow(X)) {
    B_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -5, 0,
                  ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -1.7 &
                           (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -5, 0.0189,
                         ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 1.7 &
                                  (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -1.7, 0.1138,
                                ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 5 &
                                         (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > 1.7, 0.0190, 0))))
    x_matrix <- (1 + delta[i]) * B_i * (X[i, ] %*% t(X[i, ]))
    X_matrix <- X_matrix + x_matrix
  }
  Sigma_inv <- diag(v_0, p) + 2 * expect_inverse_b_2 * X_matrix
  if(nrow(Sigma_inv) == 1) Sigma <- matrix(1 / Sigma_inv, nrow = 1)
  else(Sigma <- inv(Sigma_inv))
  return(Sigma)
}

get_mu_star <- function(curr_b, curr_mu, alpha, omega, delta, X, v_0, Sigma, mu_0, Y){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  expect_inverse_b <- alpha / omega
  p <- ncol(X)
  YX_matrix <- matrix(0, nrow = 1, ncol = p)
  for (i in 1:nrow(X)) {
    A_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -5, 0,
                  ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -1.7 &
                           (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -5, 0.1696,
                         ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 1.7 &
                                  (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -1.7, 0.5,
                                ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 5 &
                                         (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > 1.7, 0.8303, 1))))
    B_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -5, 0,
                  ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= -1.7 &
                           (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -5, 0.0189,
                         ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 1.7 &
                                  (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > -1.7, 0.1138,
                                ifelse((Y[i] - sum(X[i, ] * curr_mu)) / curr_b <= 5 &
                                         (Y[i] - sum(X[i, ] * curr_mu)) / curr_b > 1.7, 0.0190, 0))))
    yx_matrix <- (-delta[i] + (1 + delta[i]) * A_i) * expect_inverse_b * X[i, ] +
      2 * (1 + delta[i]) * B_i * expect_inverse_b_2 * Y[i] * X[i, ]
    YX_matrix <- YX_matrix + yx_matrix
  }
  YX_matrix <- v_0 * mu_0 + YX_matrix
  mu <- YX_matrix %*% Sigma
  return(mu)
}
