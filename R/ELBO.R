## ELBO Calculations ===================================================

#' Calculates the approximated expectation over the density of observed
#' data \emph{D} given parameters \emph{β} and \emph{b},
#' \eqn{\log(p(D|\beta,b))}.
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param alpha Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution  of \emph{b}.
#' @param curr_mu The current value of the parameter \eqn{\mu^*} of the
#'  approximate posterior distribution of \emph{β}.
#' @param expectation_b The expected value of b.
#' @returns The approximated log-likelihood \eqn{\log(p(D|\beta,b))}.
#'
#' @noRd
expectation_log_likelihood <- function(y, X, delta, alpha, omega,
                                       curr_mu, expectation_b) {
  res <- 0
  for (i in 1:nrow(X)) {
    bz_i <- y[i] - sum(X[i, ] * curr_mu)
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

#' Calculates \eqn{\diff_\beta}, the difference between the expectations
#' of \eqn{\log(p(\beta))} and \eqn{\log(q(\beta))}.
#'
#' @param v_0 Precision hyperparameter \eqn{v_0} of the prior
#'  distribution of \emph{β}.
#' @param mu_0 Hyperparameter \eqn{\mu_0} of the prior distribution of
#'  \emph{β}.
#' @param curr_mu The current value of the parameter \eqn{\mu^*} of the
#'  approximate posterior distribution of \emph{β}.
#' @param Sigma Parameter \eqn{\Sigma^*} of the approximate posterior
#'  distribution of \emph{β}.
#' @returns The difference between the expectations of
#'  \eqn{\log(p(\beta))} and \eqn{\log(q(\beta))}.
#'
#' @noRd
diff_beta <- function(mu_0, v_0, curr_mu, Sigma) {
  res <- sum(diag(Sigma)) + sum((curr_mu - mu_0) * (curr_mu - mu_0))
  diff_beta <- (-v_0 * res + log(det(Sigma))) / 2
  diff_beta
}

#' Calculate \eqn{\text{diff}_b}, the difference between the
#'  expectations of \eqn{\log(p(b))} and \eqn{\log(q(b))}.
#'
#' @param alpha_0 Hyperparameter \eqn{\alpha_0} of the prior
#'  distribution of \emph{b}.
#' @param omega_0 Hyperparameter \eqn{\omega_0} of the prior
#'  distribution of \emph{b}.
#' @param alpha Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @returns The difference between the expectations of \log(p(b)) and
#'  \log(q(b)).
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

#' Calculates the variational Bayes convergence criteria, evidence lower
#' bound (ELBO), optimized in `survregVB.fit`.
#'
#' @name elbo
#'
#' @param y A vector of observed log-transformed survival times.
#' @param X A matrix of predictors (covariates), including an intercept.
#' @param delta A binary vector indicating right censoring.
#' @param alpha_0 Hyperparameter \eqn{\alpha_0} of the prior
#'  distribution of \emph{b}.
#' @param omega_0 Hyperparameter \eqn{\omega_0} of the prior
#'  distribution of \emph{b}.
#' @param v_0 Precision hyperparameter \eqn{v_0} of the prior
#'  distribution of \emph{β}.
#' @param mu_0 Hyperparameter \eqn{\mu_0} of the prior distribution of
#'  \emph{β}.
#' @param alpha Parameter \eqn{\alpha^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the approximate posterior
#'  distribution of \emph{b}.
#' @param curr_mu The current value of the parameter \eqn{\mu^*} of the
#'  approximate posterior distribution of \emph{β}.
#' @param Sigma Parameter \eqn{Sigma^*} of the approximate posterior
#'  distribution of\emph{β}.
#' @param expectation_b The expected value of b.
#' @returns The evidence lower bound (ELBO).
#'
#' @export
#' @seealso \code{\link{survregVB.fit}}
elbo <- function(y, X, delta, alpha_0, omega_0, mu_0, v_0, alpha, omega,
                 curr_mu, Sigma, expectation_b) {
  expectation_log_likelihood <-
    expectation_log_likelihood(y, X, delta, alpha, omega, curr_mu,
                               expectation_b)
  diff_beta <- diff_beta(mu_0, v_0, curr_mu, Sigma)
  diff_b <- diff_b(alpha_0, omega_0, alpha,  omega)

  elbo <- expectation_log_likelihood + diff_beta + diff_b
  elbo
}
