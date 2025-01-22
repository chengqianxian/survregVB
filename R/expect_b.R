## Expectation of b Calculations =======================================

#' Calculates the expected value of
#' \eqn{b\sim\text{Inverse-Gamma}(\alpha,\omega)}.
#'
#' @param alpha Parameter \eqn{\alpha^*} of the distribution of
#'  \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the distribution of
#'  \emph{b}.
#' @returns The expectation of b.
#' @noRd
expectation_b <- function(alpha,
                          omega) {
  b = omega / (alpha - 1)
  return(b)
}

#' Calculates the expected value of the inverse of
#' \eqn{b\sim\text{Inverse-Gamma}(\alpha,\omega)}.
#'
#' @param alpha Parameter \eqn{\alpha^*} of the distribution of
#'  \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the distribution of
#'  \emph{b}.
#' @returns The expectation of \eqn{1/b}.
#' @noRd
expectation_inverse_b <- function(alpha,
                                  omega) {
  inv_b = alpha / omega
  return(inv_b)
}

#' Calculates the expected value of the squared inverse of
#'\eqn{b\sim\text{Inverse-Gamma}(\alpha,\omega)}.
#'
#' @param alpha Parameter \eqn{\alpha^*} of the distribution of
#'  \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the distribution of
#'  \emph{b}.
#' @return The expectation of \eqn{1/b^2}.
#' @noRd
expectation_inverse_b_2 <- function(alpha,
                                    omega) {
  inv_b_2 = (alpha + alpha^2) / omega^2
  return(inv_b_2)
}

#' Calculates the expected value of the log of
#'\eqn{b\sim\text{Inverse-Gamma}(\alpha,\omega)}.
#'
#' @param alpha Parameter \eqn{\alpha^*} of the distribution of
#'  \emph{b}.
#' @param omega Parameter \eqn{\omega^*} of the distribution of
#'  \emph{b}.
#' @return The expectations of \eqn{\log(b)}.
#' @noRd
expectation_log_b <- function(alpha,
                              omega) {
  log_b <- log(omega) - digamma(alpha)
  return(log_b)
}
