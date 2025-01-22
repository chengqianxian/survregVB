## Parameter Update Calculations =======================================

#' Calculates parameter \eqn{\alpha^*} of the approximate posterior
#' distribution of \emph{b} in order to optimize the evidence based
#' lower bound (ELBO) in `survregVB.fit`.
#'
#' @param alpha_0 Hyperparameter \eqn{\alpha_0} of the prior
#'  distribution of \emph{b}.
#' @param delta A binary vector indicating right censoring.
#' @returns Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{elbo}}
alpha_star <- function(alpha_0,
                       delta) {
  r <- sum(delta)
  alpha <- alpha_0 + r
  return(alpha)
}


#' Calculates parameter \eqn{\omega^*} of the approximate posterior
#' distribution of \emph{b} in order to optimize the evidence based
#' lower bound (ELBO) in `survregVB.fit`.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param omega_0 Hyperparameter \eqn{\omega_0} of the prior
#'  distribution of \emph{b}.
#' @param mu Parameter \eqn{\mu^*} of the approximate posterior
#'  distribution of \emph{β}.
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution of \emph{b}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{elbo}}
omega_star <- function(y,
                       X,
                       delta,
                       omega_0,
                       mu,
                       expectation_b) {
  res <- 0
  for (i in 1:nrow(X)) {
    bz_i <- y[i] - sum(X[i, ] * mu)
    z_i <- bz_i / expectation_b
    phi_i <- get_phi(z_i)
    res_i <- (delta[i] - (1 + delta[i]) * phi_i) * bz_i
    res <- res + res_i
  }
  omega <- omega_0 - res
  return(omega)
}


#' Calculates parameter \eqn{\mu^*} of the approximate posterior
#' distribution of \emph{β} in order to optimize the evidence based
#' lower bound (ELBO) in `survregVB.fit`.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param mu_0 Hyperparameter \eqn{\mu_0} of the prior distribution of
#'  \emph{β}.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of
#'  \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param curr_mu The current value of the parameter \eqn{\mu^*} of the
#'  approximate posterior distribution of \emph{β}.
#' @param Sigma Parameter \eqn{\Sigma^*} of the approximate posterior
#'  distribution of\emph{β}.
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\mu^*} of the approximate posterior
#'  distribution of \emph{β}.
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{elbo}}
mu_star <- function(y,
                    X,
                    delta,
                    mu_0,
                    v_0,
                    alpha,
                    omega,
                    curr_mu,
                    Sigma,
                    expectation_b) {
  p <- ncol(X)
  expectation_inverse_b <- expectation_inverse_b(alpha, omega)
  expectation_inverse_b_2 <- expectation_inverse_b_2(alpha, omega)

  yX_matrix <- matrix(0, nrow = 1, ncol = p)
  for (i in 1:nrow(X)) {
    z_i <- ((y[i] - sum(X[i, ] * curr_mu))) / expectation_b
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
  return(mu)
}

#' Calculates parameter \eqn{\Sigma^*} of the approximate posterior
#' distribution of \emph{β} in order to optimize the evidence based
#' lower bound (ELBO) in `survregVB.fit`.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param v_0 Hyperparameter \eqn{v_0} of the prior distribution of
#'  \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param curr_mu The current value of the parameter \eqn{\mu^*} of the
#'  approximate posterior distribution of \emph{β}.
#' @param expectation_b The expected value of b.
#' @returns Parameter \eqn{\Sigma^*} of the approximate posterior
#'  distribution of \emph{β}.
#'
#' @importFrom matlib inv
#' @export
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{elbo}}
sigma_star <- function(y,
                       X,
                       delta,
                       v_0,
                       alpha,
                       omega,
                       curr_mu,
                       expectation_b) {
  p <- ncol(X)
  X_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in 1:nrow(X)) {
    z_i <- (y[i] - sum(X[i, ] * curr_mu)) / expectation_b
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

  return(Sigma)
}



