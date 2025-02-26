## Expectation Calculations=============================================

#' Calculates the expected value of \eqn{q^*(b)}, an
#' \eqn{\text{Inverse-Gamma}(\alpha,\omega)} density function.
#'
#' @inheritParams elbo
#' @returns The expectation of b.
#' @noRd
expectation_b <- function(alpha, omega) {
  b <- omega / (alpha - 1)
  b
}

#' Calculates the expected value of the inverse of \eqn{q^*(b)}, an
#' \eqn{\text{Inverse-Gamma}(\alpha,\omega)} density function.
#'
#' @inheritParams elbo
#' @returns The expectation of \eqn{1/b}.
#' @noRd
expectation_inverse_b <- function(alpha, omega) {
  inv_b <- alpha / omega
  inv_b
}

#' Calculates the expected value of the squared inverse of \eqn{q^*(b)}, an
#' \eqn{\text{Inverse-Gamma}(\alpha,\omega)} density function.
#'
#' @inheritParams elbo
#' @return The expectation of \eqn{1/b^2}.
#' @noRd
expectation_inverse_b_2 <- function(alpha, omega) {
  inv_b_2 <- (alpha + alpha^2) / omega^2
  inv_b_2
}

#' Calculates the expected value of the log of \eqn{q^*(b)}, an
#' \eqn{\text{Inverse-Gamma}(\alpha,\omega)} density function.
#'
#' @inheritParams elbo
#' @return The expectations of \eqn{\log(b)}.
#' @noRd
expectation_log_b <- function(alpha, omega) {
  log_b <- log(omega) - digamma(alpha)
  log_b
}
