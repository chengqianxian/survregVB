## Parameter Calculations ==============================================

#' Calculates parameter \eqn{\alpha^*} of \eqn{q^*(b)} to optimize the evidence
#' based lower bound (ELBO) in \code{survregVB.fit} and \code{survregVB.frailty.fit}.
#'
#' @param alpha_0 Hyperparameter \eqn{\alpha_0} of the prior distribution
#' of \emph{b}.
#' @param delta A binary vector indicating right censoring.
#' @returns Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{survregVB.frailty.fit}}
alpha_star <- function(alpha_0, delta) {
  r <- sum(delta)
  alpha <- alpha_0 + r
  alpha
}

#' Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the evidence
#' based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param omega_0 Hyperparameter \eqn{\omega_0} of the prior distribution
#'  of \emph{b}.
#' @param mu Parameter \eqn{\mu^*} of the approximate posterior distribution
#'  of \emph{β}. \eqn{q}
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
omega_star <- function(y, X, delta, omega_0, mu, expectation_b) {
  res <- 0
  for (i in 1:nrow(X)) {
    bz_i <- y[i] - sum(X[i, ] * mu)
    z_i <- bz_i / expectation_b
    phi_i <- get_phi(z_i)
    res_i <- (delta[i] - (1 + delta[i]) * phi_i) * bz_i
    res <- res + res_i
  }
  omega <- omega_0 - res
  omega
}

#' Calculates parameter \eqn{\mu^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param mu_0 Hyperparameter \eqn{\mu_0} of the prior distribution of \emph{β}.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param Sigma Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
mu_star <- function(y, X, delta, mu_0, v_0, alpha, omega, mu, Sigma, expectation_b) {
  p <- ncol(X)
  expectation_inverse_b <- expectation_inverse_b(alpha, omega)
  expectation_inverse_b_2 <- expectation_inverse_b_2(alpha, omega)

  yX_matrix <- matrix(0, nrow = 1, ncol = p)
  for (i in 1:nrow(X)) {
    z_i <- (y[i] - sum(X[i, ] * mu)) / expectation_b
    rho_i <- get_rho(z_i)
    zeta_i <- get_zeta(z_i)

    res_i_a <- expectation_inverse_b *
      (-delta[i] + (1 + delta[i]) * rho_i) * X[i, ]
    res_i_b <- expectation_inverse_b_2 *
      (1 + delta[i]) * y[i] * zeta_i * X[i, ]
    yX_matrix_i <- res_i_a + 2 * res_i_b
    yX_matrix <- yX_matrix + yX_matrix_i
  }

  yX_matrix <- v_0 * mu_0 + yX_matrix
  mu <- yX_matrix %*% Sigma
  mu
}

#' Calculates parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#'
#' @importFrom matlib inv
#' @export
#' @seealso \code{\link{survregVB.fit}}
Sigma_star <- function(y, X, delta, v_0, alpha, omega, mu,
                       expectation_b) {
  p <- ncol(X)
  X_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in 1:nrow(X)) {
    z_i <- (y[i] - sum(X[i, ] * mu)) / expectation_b
    zeta_i <- get_zeta(z_i)
    X_matrix_i <- (1 + delta[i]) * zeta_i * (X[i, ] %*% t(X[i, ]))
    X_matrix <- X_matrix + X_matrix_i
  }

  expectation_inverse_b_2 <- expectation_inverse_b_2(alpha, omega)
  Sigma_inv <- diag(v_0, p) + 2 * expectation_inverse_b_2 * X_matrix

  if (nrow(Sigma_inv) == 1) {
    Sigma <- matrix(1 / Sigma_inv, nrow = 1)
  } else {
    Sigma <- inv(Sigma_inv)
  }
  Sigma
}

## With cluster index ==================================================

#' Calculates the clustered survival times by adjusting the survival times
#' to account for shared frailty.

#' @param y A vector of observed log-transformed survival times.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(sigma^2_{\gamma})} for
#'  all clusters.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#'
#' @returns A vector of log-transformed survival time adjusted for the frailty
#' of the corresponding cluster.
#'
#' @noRd
get_cluster_y <- function(y, tau, cluster) {
  for (i in 1:length(y)) {
    y[i] = y[i] - tau[cluster[i]]
  }
  y
}

#' Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the evidence
#' based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param omega_0 Hyperparameter \eqn{\omega_0} of the prior distribution
#'  of \emph{b}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(sigma^2_{\gamma})} for
#'  all clusters.
#' @param expectation_b The expected value of b.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns Parameter \eqn{\omega^*} of the approximate posterior distribution
#'  of \emph{b}.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
omega_star_cluster <- function(y, X, delta, omega_0, mu, tau, expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  res <- 0
  for (i in 1:nrow(X)) {
    z_i <- (y_cluster[i] - sum(X[i, ] * mu)) / expectation_b
    phi_i <- get_phi(z_i)
    bz_i <- y[i] - sum(X[i, ] * mu)
    res_i <- (delta[i] - (1 + delta[i]) * phi_i) * bz_i
    res <- res + res_i
  }
  omega <- omega_0 - res
  omega
}

#' Calculates parameter \eqn{\mu^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param mu_0 Hyperparameter \eqn{\mu_0} of the prior distribution of \emph{β}.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param Sigma Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(sigma^2_{\gamma})} for
#'  all clusters.
#' @param expectation_b The expected value of b.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
mu_star_cluster <- function(y, X, delta, mu_0, v_0, alpha, omega, mu, Sigma,
                            tau, expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  mu_star(y_cluster, X, delta, mu_0, v_0, alpha, omega, mu, Sigma, expectation_b)
}

#' Calculates parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(sigma^2_{\gamma})} for
#'  all clusters.
#' @param expectation_b The expected value of b.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
Sigma_star_cluster <- function(y, X, delta, v_0, alpha, omega, mu, tau,
                               expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  Sigma_star(y_cluster, X, delta, v_0, alpha, omega, mu, expectation_b)
}

#' Calculates parameter \eqn{sigma^{2*}} of \eqn{q^*(\gamma)} for all clusters
#' to optimize the evidence based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu The parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(\gamma)} for all clusters.
#' @param lambda Parameter \eqn{\lambda^*} of \eqn{q^*(sigma^2_{\gamma})}
#'  for all clusters.
#' @param eta Parameter \eqn{\eta^*} of \eqn{q^*(sigma^2_{\gamma})} for
#'  all clusters.
#' @param expectation_b The expected value of b.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns Parameter \eqn{\sigma^{2*}} of \eqn{q^*(\gamma)} for all clusters.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
sigma_star <- function(y, X, delta, alpha, omega, mu, tau,  lambda, eta,
                       expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  expectation_inverse_b_2 <- expectation_inverse_b_2(alpha, omega)
  expectation_inverse_sigma <- expectation_inverse_sigma(lambda, eta)
  zeta <- numeric(nrow(X))

  for (i in 1:nrow(X)) {
    z_i <- (y_cluster[i] - sum(X[i, ] * mu)) / expectation_b
    zeta[i] <- get_zeta(z_i)
  }

  K = length(unique(cluster))
  sigma <- numeric(K)
  for (k in 1:K) {
    delta_k <- delta[cluster == k]
    zeta_k <- zeta[cluster == k]
    sigma[k] <-  expectation_inverse_sigma + 2 *
      expectation_inverse_b_2 * sum((1 + delta_k) * zeta_k)
    sigma[k] <- 1 / sigma[k]
  }
  sigma
}

#' Calculates parameter \eqn{tau^*_i} of \eqn{q^*(\gamma)} for all clusters
#' to optimize the evidence based lower bound (ELBO) in
#' \code{survregVB.frailty.fit}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param alpha Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#' @param omega Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#' @param mu The parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#' @param tau Parameter \eqn{\tau^*} of
#'  \eqn{q^*(\gamma)} for all clusters.
#' @param Parameter \eqn{\sigma^{2*}} of \eqn{q^*(\gamma)}} for all clusters.
#' @param expectation_b The expected value of b.
#' @param cluster A numeric vector indicating the cluster assignment for
#'  each observation.
#' @returns Parameter \eqn{\tau^*} of \eqn{q^*(\gamma)} for all clusters.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
tau_star <- function(y, X, delta, alpha, omega, mu, tau, sigma,
                     expectation_b, cluster) {
  y_cluster <- get_cluster_y(y, tau, cluster)
  expectation_inverse_b <- expectation_inverse_b(alpha, omega)
  expectation_inverse_b_2 <- expectation_inverse_b_2(alpha, omega)
  n = nrow(X)
  zeta <- numeric(n)
  rho <- numeric(n)

  for (i in 1:n) {
    # zeta_i is the coefficient of quadratic approximation
    z_i <- (y_cluster[i] - sum(X[i, ] * mu)) / expectation_b
    zeta[i] <- get_zeta(z_i)
    rho[i] <- get_rho(z_i)
  }

  K = length(unique(cluster))
  tau <- numeric(K)
  for (k in 1:K) {
    delta_k <- delta[cluster == k]
    zeta_k <- zeta[cluster == k]
    rho_k <- rho[cluster == k]
    y_k <- y[cluster == k]
    X_k <- X[cluster == k, ]
    tau_k <- 0

    for (i in 1:nrow(X_k)) {
      res_i_a <- expectation_inverse_b *
        (-delta_k[i] + (1 + delta_k[i]) * rho_k[i])
      res_i_b <- expectation_inverse_b_2 * (1 + delta_k[i]) *
        zeta_k[i] * (y_k[i] - sum(X_k[i, ] * mu))
      tau_k <- tau_k + res_i_a + 2 * res_i_b
    }
    tau[k] <- tau_k * sigma[k]
  }
  tau
}

#' Calculates parameter \eqn{lambda^*} of \eqn{q^*(sigma^2_{\gamma})} for
#' all clusters to optimize the evidence based lower bound (ELBO) in
#' \code{survregVB.frailty.fit}.
#'
#' @param lambda_0 Hyperparameter \eqn{\lambda_0} of the prior distribution
#'  of \eqn{sigma^2_{\gamma}}.
#' @param K The number of clusters.
#' @return Parameter \eqn{lambda^*} of \eqn{q^*(sigma^2_{\gamma})} for all
#'  clusters.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
lambda_star <- function(lambda_0, K) {
  lambda <- lambda_0 + K / 2
  lambda
}

#' Calculates parameter \eqn{eta^*} of \eqn{q^*(sigma^2_{\gamma})} for all
#' clusters to optimize the evidence based lower bound (ELBO) in
#' \code{survregVB.frailty.fit}.
#'
#' @param eta_0 Hyperparameter \eqn{\eta_0} of the prior distribution of
#'  \eqn{sigma^2_{\gamma}}.
#' @param tau Parameter \eqn{\tau^*} of \eqn{q^*(\gamma)} for all clusters.
#' @param sigma Parameter \eqn{\sigma^{2*}} of \eqn{q^*(\gamma)} for all
#'  clusters.
#' @return Parameter \eqn{eta^*} of \eqn{q^*(sigma^2_{\gamma})} for all
#'  clusters.
#'
#' @export
#' @seealso \code{\link{survregVB.frailty.fit}}
eta_star <- function(eta_0, tau, sigma) {
  eta <- eta_0 + 0.5 * sum(sigma + tau^2)
  eta
}
