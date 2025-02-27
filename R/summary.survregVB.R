#' Summary for Variational Bayes log-logistic AFT models.
#'
#' Produces a summary of a fitted Variational Bayes Parametric Survival
#' Regression Model for a Log-Logistic AFT Model
#'
#' @name summary.survregVB
#' @aliases print.summary.survregVB
#'
#' @param object The result of a \code{survregVB} fit.
#' @param ci The significance level for the credible intervals. (Default:0.95).
#' @param \dots For future arguments.
#'
#' @returns An object of class \code{summary.survregVB} with components:
#' \itemize{
#'   \item \code{ELBO}: The final value of the Evidence Lower Bound (ELBO)
#'    at the last iteration.
#'   \item \code{alpha}: The shape parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#'   \item \code{omega}: The scale parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#'   \item \code{mu}: Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}, a vector
#'    of means.
#'   \item \code{Sigma}: Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}, a
#'    covariance matrix.
#'   \item \code{na.action}: A missing-data filter function, applied to the
#'    \code{model.frame}, after any subset argument has been used.
#'   \item \code{iterations}: The number of iterations performed by the VB
#'    algorithm: before converging or reaching \code{max_iteration}.
#'   \item \code{n}: The number of observations.
#'   \item \code{call}: The function call used to invoke the \code{survregVB}
#'    method.
#'   \item \code{not_converged}: A boolean indicating if the algorithm
#'    converged.
#'   \itemize{
#'     \item \code{TRUE}: If the algorithm did not converge prior to
#'      achieving `max_iteration`.
#'     \item \code{NULL}: If the algorithm converged successfully.
#'   }
#'   \item \code{posteriors}: A matrix with one row for each coefficient,
#'    and one row for the scale parameter, and columns containing:
#'   \itemize{
#'     \item \code{Value}: The estimated value
#'     \item \code{Lower CI}: The lower bound of the credible interval.
#'     \item \code{Upper CI}: The upper bound of the credible interval.
#'   }
#' }
#'
#' @method summary survregVB
#' @export
#' @seealso \code{\link{survregVB}}
summary.survregVB <- function(object, ci = 0.95, ...) {
  if (!is.null(object$not_converged)) {
    warning("SurvregVB did not converge.\n")
  }

  mu <- object$mu

  # column names for posteriors
  ci_label <- paste0(100 * ci, "% CI")
  col_names <- c("Value", "SD", paste(ci_label, c("Upper", "Lower")))

  # regression coefficients (beta)
  coefficients <- matrix(rep(mu, 4), ncol = 4)
  # Std. Dev
  coefficients[, 2] <- sqrt(diag(object$Sigma))

  # Cred. Int
  beta_ci <- beta_ci(mu, object$Sigma, ci)
  coefficients[, 2] <- beta_ci[, 1]
  coefficients[, 3] <- beta_ci[, 2]

  # scale parameter (b)
  alpha <- object$alpha
  omega <- object$omega
  scale_sd <- sqrt((omega^2) / ((alpha - 1)^2 * (alpha - 2)))
  scale <- c(omega / (alpha - 1), scale_sd, b_ci(alpha, omega, ci))

  if (!is.null(object$clustered)) {
    # for models with shared frailty
    # random intercept (sigma_gamma squared)
    lambda <- object$lambda
    eta <- object$eta
    intercept_sd <- sqrt((eta^2) / ((lambda - 1)^2 * (eta - 2)))
    intercept <- c(eta / (lambda - 1), intercept_sd, b_ci(lambda, eta, ci))

    posteriors <- rbind(coefficients, scale, intercept)
    row_names <- c(names(mu), "scale", "intercept")

    dimnames(posteriors) <- list(row_names, col_names)

    x <- object[match(
      c(
        "clustered", "call", "ELBO", "alpha", "omega",
        "mu", "Sigma", "tau", "sigma", "lambda", "eta",
        "na.action", "iterations", "n"
      ),
      names(object),
      nomatch = 0
    )]
  } else {
    # for models without shared frailty
    posteriors <- rbind(coefficients, scale)
    row_names <- c(names(mu), "scale")

    dimnames(posteriors) <- list(row_names, col_names)

    x <- object[match(
      c(
        "call", "ELBO", "alpha", "omega", "mu", "Sigma",
        "na.action", "iterations", "n"
      ),
      names(object),
      nomatch = 0
    )]
  }

  x <- c(x, list(posteriors = posteriors))
  class(x) <- "summary.survregVB"

  x
}
