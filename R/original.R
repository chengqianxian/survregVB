### update Sigma which is the covariance matrix of beta (dimension: p*p)
# cluster_index is vector with length of sample size (nrow(X))
# curr_b is the posterior mean of b, same for curr_mu
Mod_1_get_Sigma_star <- function(X, Y, curr_b, curr_mu, curr_tau, alpha, omega, delta, cluster_index, v_0){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  p <- ncol(X)
  X_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in 1:nrow(X)) {
    # zeta_i is the coefficient of quadratic approximation
    # zeta_i has included the information from cluster index
    zeta_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                     ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                              (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0189,
                            ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                     (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.1138,
                                   ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                            (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.0190, 0))))
    x_matrix <- (1 + delta[i]) * zeta_i * (X[i, ] %*% t(X[i, ]))
    X_matrix <- X_matrix + x_matrix
  }
  Sigma_inv <- diag(v_0, p) + 2 * expect_inverse_b_2 * X_matrix
  if(nrow(Sigma_inv) == 1) Sigma <- matrix(1 / Sigma_inv, nrow = 1)
  else(Sigma <- inv(Sigma_inv))
  return(Sigma)
}

### update mu which is the mean vector of beta coefficients
Mod_1_get_mu_star <- function(X, Y, curr_b, curr_mu, curr_tau, alpha, omega, delta, cluster_index, v_0, mu_0, Sigma){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  expect_inverse_b <- alpha / omega
  p <- ncol(X)
  YX_matrix <- matrix(0, nrow = 1, ncol = p)
  for (i in 1:nrow(X)) {
    rho_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                    ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                             (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.1696,
                           ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                    (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.5,
                                  ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                           (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.8303, 1))))
    zeta_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                     ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                              (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0189,
                            ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                     (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.1138,
                                   ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                            (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.0190, 0))))
    yx_matrix <- (-delta[i] + (1 + delta[i]) * rho_i) * expect_inverse_b * X[i, ] +
      2 * (1 + delta[i]) * zeta_i * expect_inverse_b_2 * (Y[i] - curr_tau[cluster_index[i]]) * X[i, ]
    YX_matrix <- YX_matrix + yx_matrix
  }
  YX_matrix <- v_0 * mu_0 + YX_matrix
  mu <- YX_matrix %*% Sigma
  return(mu)
}

### update sigma_i which is the variance for the q(gamma_i)
## sigma here is the variance, not the square root of the variance
## in our manuscript, we use sigma_i^2
## we update sigma_i for i=1,...,K at the same time
## the updated results are saved as a vector, sigma, for all clusters
Mod_1_get_sigma_star <- function(X, Y, curr_b, curr_mu, curr_tau, alpha, omega, lambda, eta, delta, cluster_index){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  expect_inverse_sigma <- lambda / eta
  zeta <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    # zeta_i is the coefficient of quadratic approximation
    zeta[i] <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                      ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                               (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0189,
                             ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                      (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.1138,
                                    ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                             (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.0190, 0))))
  }
  sigma <- numeric(length(unique(cluster_index)))
  for (k in 1:length(unique(cluster_index))) {
    delta_k <- delta[cluster_index == k]
    zeta_k <- zeta[cluster_index == k]
    sigma[k] <-  expect_inverse_sigma + 2 * expect_inverse_b_2 * sum((1 + delta_k) * zeta_k) #### + to *
    sigma[k] <- 1 / sigma[k]
  }
  return(sigma)
}

### update tau_i which is the variance for the q(gamma_i), that is the update the current tau
## we update tau_i for i=1,...,K at the same time
## the updated results are saved as a vector, tau, for all clusters

Mod_1_get_tau_star <- function(X, Y, curr_b, curr_mu, curr_tau, alpha, omega, sigma, delta, cluster_index){
  expect_inverse_b_2 <- (alpha + alpha^2) / omega^2
  expect_inverse_b <- alpha / omega
  zeta <- numeric(nrow(X))
  rho <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    # zeta_i is the coefficient of quadratic approximation
    zeta[i] <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                      ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                               (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0189,
                             ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                      (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.1138,
                                    ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                             (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.0190, 0))))
    rho[i] <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                     ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -1.7 &
                              (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.1696,
                            ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.7 &
                                     (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.7, 0.5,
                                   ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                            (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.7, 0.8303, 1))))
  }
  tau <- numeric(length(unique(cluster_index)))
  for (k in 1:length(unique(cluster_index))) {
    delta_k <- delta[cluster_index == k]
    zeta_k <- zeta[cluster_index == k]
    rho_k <- rho[cluster_index == k]
    y_k <- Y[cluster_index == k]
    X_k <- X[cluster_index == k, ]
    tau_k <- 0
    for (i in 1:nrow(X_k)) {
      tau_sub <- expect_inverse_b * ((1 + delta_k[i]) * rho_k[i] - delta_k[i]) +
        2 * expect_inverse_b_2 * (1 + delta_k[i]) * zeta_k[i] * (y_k[i] - sum(X_k[i, ] * curr_mu))
      tau_k <- tau_k + tau_sub
    }
    tau[k] <- tau_k * sigma[k]
  }
  return(tau)
}

### Update omega which is the scale parameter in q(b)
Mod_1_get_omega_star <- function(X, Y, curr_b, curr_mu, curr_tau, delta, cluster_index, omega_0){
  res <- 0
  for (i in 1:nrow(X)) {
    phi_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                    ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]])/ curr_b <= -1.701 &
                             (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0426,
                           ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 0 &
                                    (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.701, 0.3052,
                                  ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.702 &
                                           (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 0, 0.6950,
                                         ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                                  (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.702, 0.9574, 1)))))
    omega_i <- (delta[i] - (1 + delta[i]) * phi_i) * (Y[i] - sum(X[i, ] * curr_mu))
    res <- res + omega_i
  }
  omega <- omega_0 - res
  return(omega)
}

### Update eta which is scale parameter in q(sigma_gamma)
Mod_1_get_eta_star <- function(curr_sigma, curr_tau, eta_0){
  eta <- eta_0 + 0.5 * sum(curr_sigma + curr_tau^2)
  return(eta)
}

# This R script file calculates the ELBO in the model with the random intercept
Mod_1_vb.aft.get.elbo <- function(Y, X, delta, mu_0, v_0, alpha_0, omega_0, lambda_0, eta_0, curr_mu, Sigma, curr_tau, curr_sigma, curr_b, alpha, omega, lambda, eta, cluster_index){
  res <- 0
  for (i in 1:nrow(X)) {
    phi_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= -5, 0,
                    ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]])/ curr_b <= -1.701 &
                             (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -5, 0.0426,
                           ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 0 &
                                    (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > -1.701, 0.3052,
                                  ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 1.702 &
                                           (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 0, 0.6950,
                                         ifelse((Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b <= 5 &
                                                  (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]]) / curr_b > 1.702, 0.9574, 1)))))
    res_i <- (delta[i] - phi_i * (1 + delta[i])) * (Y[i] - sum(X[i, ] * curr_mu) - curr_tau[cluster_index[i]])
    res <- res + res_i
  }
  r <- sum(delta)
  expect_inverse_b <- alpha / omega
  expect_log_b <- log(omega) - digamma(alpha)
  expect_log_likelihood <- -r * expect_log_b + expect_inverse_b * res # the first term in the ELBO

  ### ----------
  diff_beta <- (- v_0 * (sum(diag(Sigma)) + sum((curr_mu - mu_0) * (curr_mu - mu_0))) + log(det(Sigma))) / 2 # the second term in the ELBO

  ### ----------
  K <- length(unique(cluster_index)) # number of groups
  expect_log_sigma_gamma <- log(eta) - digamma(lambda)
  expect_sigmma_gamma_gamma <- lambda / eta * sum(curr_sigma + curr_tau^2)
  log_sigma_star_sum <- sum(log(curr_sigma))

  diff_gamma <- -0.5 * K * expect_log_sigma_gamma - 0.5 * expect_sigmma_gamma_gamma - 0.5 * log_sigma_star_sum # the third term

  ### ----------
  diff_b <- (alpha - alpha_0) * expect_log_b + (omega - omega_0) * expect_inverse_b - alpha * log(omega) # the fourth term

  ### ----------
  # the fifth term
  diff_sigmma_gamma <- (lambda - lambda_0) * expect_log_sigma_gamma + (eta - eta_0) * lambda / eta - lambda * log(eta)

  elbo <- expect_log_likelihood + diff_beta + diff_gamma + diff_b + diff_sigmma_gamma
  return(elbo)
}

Mod_1_vb_aft_final <- function(Y, X, delta, cluster_index, mu_0, v_0, alpha_0, omega_0, lambda_0, eta_0, max_iteration = 100, threshold = 0.01){
  Y <- log(Y)
  n <- nrow(X)
  p <- ncol(X)
  K <- length(unique(cluster_index))

  alpha <- alpha_0 + sum(delta) # fixed always
  lambda <- lambda_0 + K / 2 # fixed always

  # initialization
  Curr_omega <- omega_0
  mu <- mu_0
  tau <- rep(0, K)
  Curr_eta <- eta_0


  converged <- FALSE
  iteration <- 0

  Curr_b <- Curr_omega / (alpha_0 - 1)

  curr_elbo <- 0

  while (converged == FALSE & iteration < max_iteration) {
    iteration <- iteration + 1
    Curr_Sigma <- Mod_1_get_Sigma_star(X, Y, Curr_b, mu, tau, alpha, Curr_omega, delta, cluster_index, v_0)
    Curr_mu <- Mod_1_get_mu_star(X, Y, Curr_b, mu, tau, alpha, Curr_omega, delta, cluster_index, v_0, mu_0, Curr_Sigma)
    Curr_sigma <- Mod_1_get_sigma_star(X, Y, Curr_b, Curr_mu, tau, alpha, Curr_omega, lambda, Curr_eta, delta, cluster_index)
    Curr_tau <- Mod_1_get_tau_star(X, Y, Curr_b, Curr_mu, tau, alpha, Curr_omega, Curr_sigma, delta, cluster_index)
    Curr_omega <- Mod_1_get_omega_star(X, Y, Curr_b, Curr_mu, Curr_tau, delta, cluster_index, omega_0)
    Curr_eta <- Mod_1_get_eta_star(Curr_sigma, Curr_tau, eta_0)

    elbo <- Mod_1_vb.aft.get.elbo(Y, X, delta, mu_0, v_0, alpha_0, omega_0, lambda_0, eta_0, Curr_mu, Curr_Sigma, Curr_tau, Curr_sigma, Curr_b, alpha, Curr_omega, lambda, Curr_eta, cluster_index)
    converged_1 <- ifelse(abs(elbo - curr_elbo) <= threshold, TRUE, FALSE)
    converged_2 <- ifelse(sum(abs(mu - Curr_mu)) <= threshold, TRUE, FALSE)
    converged <- ifelse(converged_1 == FALSE & converged_2 == FALSE, FALSE, TRUE)
    curr_elbo <- elbo
    mu <- Curr_mu
    tau <- Curr_tau
    Curr_b <- Curr_omega / (alpha - 1)
  }
  if(iteration > max_iteration & converged == FALSE) warning("The algorithm can not converge and the max iteration has achieved")
  return(list(
    converged_ELBO = elbo,
    iteration_used = iteration,
    alpha = alpha,
    lambda = lambda,
    Sigma = Curr_Sigma,
    mu = Curr_mu,
    sigma = Curr_sigma,
    tau = Curr_tau,
    omega = Curr_omega,
    eta = Curr_eta))
}


