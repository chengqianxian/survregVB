## ELBO Calculations ===================================================

#' Calculates the approximated expectation over the density of observed data
#' \emph{D} given parameters \emph{β} and \emph{b}, \eqn{\log(p(D|\beta,b))}.
#'
#' @inheritParams elbo
#' @returns The approximated log-likelihood \eqn{\log(p(D|\beta,b))}.
#'
#' @noRd
expectation_log_likelihood <- function(y, X, delta, alpha, omega, mu,
                                       expectation_b) {
  res <- 0
  for (i in seq_len(nrow(X))) {
    bz_i <- y[i] - sum(X[i, ] * mu)
    z_i <- bz_i / expectation_b
    phi_i <- get_phi(z_i)
    res_i <- (delta[i] - phi_i * (1 + delta[i])) * bz_i
    res <- res + res_i
  }

  r <- sum(delta)
  expectation_log_b <- expectation_log_b(alpha, omega)
  expectation_inverse_b <- expectation_inverse_b(alpha, omega)

  expectation_log_likelihood <- -r * expectation_log_b +
    expectation_inverse_b * res
  expectation_log_likelihood
}

#' Calculates the difference between the expectations of \eqn{\log(p(\beta))}
#' and \eqn{\log(q^*(\beta))}.
#'
#' @inheritParams elbo
#' @returns The difference between the expectations of \eqn{\log(p(\beta))}
#'  and \eqn{\log(q^*(\beta))}.
#'
#' @noRd
diff_beta <- function(mu_0, v_0, mu, Sigma) {
  res <- sum(diag(Sigma)) + sum((mu - mu_0) * (mu - mu_0))
  diff_beta <- (-v_0 * res + log(det(Sigma))) / 2
  diff_beta
}

#' Calculate the difference between the expectations of \eqn{\log(p(b))}
#' and \eqn{\log(q^*(b))}.
#'
#' @inheritParams elbo
#' @returns The difference between the expectations of \log(p(b)) and
#'  \log(q^*(b)).
#'
#' @noRd
diff_b <- function(alpha_0, omega_0, alpha, omega) {
  expectation_log_b <- expectation_log_b(alpha, omega)
  expectation_inverse_b <- expectation_inverse_b(alpha, omega)

  alpha_res <- (alpha - alpha_0) * expectation_log_b
  omega_res <- (omega - omega_0) * expectation_inverse_b

  diff_b <- alpha_res + omega_res - alpha * log(omega)
  diff_b
}

#' Calculates the difference between the expectations of \eqn{\log(p(\gamma))}
#' and \eqn{\log(q^*(\gamma))}.
#'
#' @inheritParams elbo_cluster
#' @returns The difference between the expectations of \eqn{\log(p(\gamma))}
#'  and \eqn{\log(q^*(\gamma))}.
#'
#' @noRd
diff_gamma <- function(tau, sigma, lambda, eta, cluster) {
  K <- length(unique(cluster))
  expectation_log_sigma_gamma <- log(eta) - digamma(lambda)
  expectation_sigmma_gamma_gamma <- lambda / eta * sum(sigma + tau^2)
  log_sigma_star_sum <- sum(log(sigma))
  diff_gamma <- -0.5 * K * expectation_log_sigma_gamma -
    0.5 * expectation_sigmma_gamma_gamma - 0.5 * log_sigma_star_sum
  diff_gamma
}

#' Calculates the difference between the expectations of
#' \eqn{\log(p(\sigma^2_gamma))} and \eqn{\log(q^*(\sigma^2_gamma))}.
#'
#' @inheritParams elbo_cluster
#' @returns The difference between the expectations of
#'  \eqn{\log(p(\sigma^2_gamma))} and \eqn{\log(q^*(\sigma^2_gamma))}.
#'
#' @noRd
diff_sigma_gamma <- function(lambda_0, eta_0, lambda, eta) {
  expectation_log_sigma_gamma <- log(eta) - digamma(lambda)
  lambda_res <- (lambda - lambda_0) * expectation_log_sigma_gamma
  eta_res <- (eta - eta_0) * lambda
  diff_sigma_gamma <- (lambda_res + eta_res / eta) - lambda * log(eta)
  diff_sigma_gamma
}

#' Calculates the variational Bayes convergence criteria, evidence lower
#' bound (ELBO), optimized in `survregVB.fit`.
#'
#' @name elbo
#'
#' @inheritParams elbo_cluster
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
elbo <- function(y, X, delta, alpha_0, omega_0, mu_0, v_0, alpha, omega,
                 mu, Sigma, expectation_b) {
  expectation_log_likelihood <-
    expectation_log_likelihood(y, X, delta, alpha, omega, mu, expectation_b)
  diff_beta <- diff_beta(mu_0, v_0, mu, Sigma)
  diff_b <- diff_b(alpha_0, omega_0, alpha, omega)

  elbo <- expectation_log_likelihood + diff_beta + diff_b
  elbo
}

#' Calculates the variational Bayes convergence criteria, evidence lower
#' bound (ELBO), optimized in `survregVB.frailty.fit`.
#'
#' @name elbo_cluster
#'
#' @inheritParams survregVB.frailty.fit
#' @param y A vector of observed log-transformed survival times.
#' @param delta A binary vector indicating right censoring.
#' @param alpha The shape parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega The scale parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}, a vector of means.
#' @param Sigma  Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}, a covariance
#'  matrix.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(\gamma_i)}, a vector of
#'  means.
#' @param sigma Parameter \eqn{\sigma^{2*}_i} of \eqn{q^*(\gamma_i)}, a
#'  vector of variance.
#' @param lambda The shape parameter \eqn{\lambda^*} of \eqn{q^*(\sigma^2_\gamma)}.
#' @param eta The scale parameter \eqn{\eta^*} of \eqn{q^*(\sigma^2_\gamma)}.
#' @param expectation_b The expected value of \emph{b}.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns The evidence lower bound (ELBO).
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
elbo_cluster <- function(y, X, delta, alpha_0, omega_0, mu_0, v_0, lambda_0,
                         eta_0, alpha, omega, mu, Sigma, tau, sigma, lambda,
                         eta, expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  expectation_log_likelihood <-
    expectation_log_likelihood(
      y_cluster, X, delta, alpha, omega, mu,
      expectation_b
    )
  diff_beta <- diff_beta(mu_0, v_0, mu, Sigma)
  diff_gamma <- diff_gamma(tau, sigma, lambda, eta, cluster)
  diff_b <- diff_b(alpha_0, omega_0, alpha, omega)
  diff_sigma_gamma <- diff_sigma_gamma(lambda_0, eta_0, lambda, eta)

  elbo <- expectation_log_likelihood + diff_beta + diff_gamma + diff_b +
    diff_sigma_gamma
  elbo
}
