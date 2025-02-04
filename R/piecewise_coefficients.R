## Piecewise Approximations ============================================

#' Get the piecewise approximation coefficient \emph{\phi} based on the
#' value of an input random variable \eqn{z_i}. Used in calculations for
#' \eqn{\omega^*} and the evidence lower bound (ELBO).
#'
#' @param z_i The value of the random variable \eqn{z_i} for the
#'  \emph{i}-th sample
#' @returns The value of \emph{phi}, the piecewise approximation
#'  coefficient based on the value of \eqn{z_i}.
#' @importFrom dplyr case_when
#' @noRd
get_phi <- function(z_i) {
  phi <- case_when(
    z_i <= -5 ~ 0,
    z_i <= -1.701 ~ 0.0426,
    z_i <= 0 ~ 0.3052,
    z_i <= 1.702 ~ 0.6950,
    z_i <= 5 ~ 0.9574,
    TRUE ~ 1
  )
  phi
}

#' Get the piecewise approximation coefficient \emph{\zeta} based on the
#' value of an input random variable \eqn{z_i}. Used in calculations for
#' \eqn{\mu^*}.
#'
#' @param z_i The value of the random variable \eqn{z_i} for the
#'  \emph{i}-th sample
#' @returns The value of \emph{zeta}, the piecewise approximation
#'  coefficient based on the value of \eqn{z_i}.
#' @importFrom dplyr case_when
#' @noRd
get_zeta <- function(z_i) {
  zeta <- case_when(
    z_i <= -5 ~ 0,
    z_i <= -1.7 ~ 0.0189,
    z_i <= 1.7 ~ 0.1138,
    z_i <= 5 ~ 0.0190,
    TRUE ~ 0
  )
  zeta
}

#' Get the piecewise approximation coefficient \emph{\rho} based on the
#' value of an input random variable \eqn{z_i}. Used in calculations for
#' \eqn{\mu^*}.
#'
#' @param z_i The value of the random variable \eqn{z_i} for the
#'  \emph{i}-th sample
#' @returns The value of \emph{rho}, the piecewise approximation
#'  coefficient based on the value of \eqn{z_i}.
#' @importFrom dplyr case_when
#' @noRd
get_rho <- function(z_i) {
  rho <- case_when(
    z_i <= -5 ~ 0,
    z_i <= -1.7 ~ 0.1696,
    z_i <= 1.7 ~ 0.5,
    z_i <= 5 ~ 0.8303,
    TRUE ~ 1
  )
  rho
}
