#' Calculates the credible interval for the regression coefficient
#' \emph{β}.
#'
#' @param mu The vector of means \eqn{\mu} of the posterior distribution
#'  of \emph{β}.
#' @param Sigma The covariance matrix \eqn{\Sigma} of the posterior
#'  distribution \emph{β}.
#' @param ci The significance level. (Default:0.95).
#' @returns Matrix containing the credible intervals for \emph{β}.
#'
#' @noRd
beta_ci <- function(mu, Sigma, ci = 0.95) {
  k <- length(mu)
  lower_bounds <- numeric(k)
  upper_bounds <- numeric(k)
  for (i in 1:k) {
    lower_bounds[i] <-
      qnorm((1 - ci) / 2, mu[i], sqrt(diag(Sigma)[i]))
    upper_bounds[i] <-
      qnorm(1 - (1 - ci) / 2, mu[i], sqrt(diag(Sigma)[i]))
  }
  ci_matrix <- cbind(CI.Lower = lower_bounds, CI.Upper = upper_bounds)
  return(ci_matrix)
}

#' Calculates the credible interval for the scale parameter \emph{b}
#' using the highest density intervals (HDI) for b.
#'
#' @param alpha The parameter \eqn{\alpha} of the posterior distribution
#'  of \emph{b}.
#' @param omega The parameter \eqn{\omega} of the posterior
#'  distribution \emph{b}.
#' @param ci The significance level. (Default:0.95).
#' @returns The credible interval for \emph{b}.
#'
#' @importFrom bayestestR hdi
#' @importFrom invgamma rinvgamma
#' @noRd
b_ci <- function(alpha, omega, ci = 0.95) {
  set.seed(100)
  posterior <- rinvgamma(100000, alpha, omega)
  lower <- hdi(posterior, ci = 0.95)$CI_low
  upper <- hdi(posterior, ci = 0.95)$CI_high
  return(c(lower, upper))
}
