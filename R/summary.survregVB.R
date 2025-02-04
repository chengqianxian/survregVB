#' Summary for Variational Bayes log-logistic AFT models.
#'
#' Produces a summary of a fitted Variational Bayes Parametric Survival
#' Regression Model for a Log-Logistic AFT Model
#'
#' @name summary.survregVB
#' @aliases print.summary.survregVB
#'
#' @param object The result of a \code{survregVB} fit.
#' @param ci The significance level for the credible intervals.
#'  (Default:0.95).
#' @param \dots For future arguments.
#'
#' @returns An object of class \code{summary.survregVB} with components:
#' \itemize{
#'   \item \code{call}: The function call used to invoke the
#'    \code{survregVB} method.
#'   \item \code{ELBO}: The final value of the Evidence Lower Bound
#'    (ELBO) at the last iteration.
#'   \item \code{alpha}: The updated parameter \eqn{\alpha} of the
#'    approximate posterior distribution of the scale parameter.
#'   \item \code{omega}: The updated parameter \eqn{\omega} of the
#'    approximate posterior distribution of the scale parameter.
#'   \item \code{mu}: The updated vector of means \eqn{\mu} of the
#'    approximate posterior distribution of the regression coefficients.
#'   \item \code{Sigma}: The updated covariance matrix \eqn{\Sigma} of
#'    the approximate posterior distribution of regression coefficients.
#'   \item \code{iterations}: The number of iterations performed by the
#'    VB algorithm: before converging or reaching \code{max_iteration}.
#'   \item \code{coefficients}: A matrix with one row for each
#'    coefficient, and columns containing:
#'   \itemize{
#'     \item \code{Value}: The estimated value of the coefficient.
#'     \item \code{Lower CI}: The lower bound of the credible interval
#'      associated with the coefficient.
#'     \item \code{Upper CI}: The upper bound of the credible interval
#'      associated with the coefficient.
#'   }
#'   \item \code{scale}: A vector containing:
#'   \itemize{
#'     \item \code{Value}: The estimated value of the scale parameter.
#'     \item \code{Lower CI}: The lower bound of the credible interval
#'      associated with the scale parameter.
#'     \item \code{Upper CI}: The upper bound of the credible interval
#'      associated with the scale parameter.
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
  cname <- names(object$mu)

  coefficients <- matrix(rep(mu, 3), ncol = 3)
  dimnames(coefficients) <- list(cname, c("Value", "Lower CI",
                                          "Upper CI"))
  beta_ci <- beta_ci(mu, object$Sigma, ci)
  coefficients[, 2] <-  beta_ci[, 1]
  coefficients[, 3] <-  beta_ci[, 2]

  alpha <- object$alpha
  omega <- object$omega
  scale <- c(expectation_b(alpha, omega), b_ci(alpha, omega))
  names(scale) <- c("Value", "Lower CI", "Upper CI")

  x <- object[match(c('call', 'ELBO', 'alpha', 'omega', 'mu', 'Sigma',
                      'iterations'),
                    names(object), nomatch=0)]
  x <- c(x, list(coefficients=coefficients, scale=scale))

  class(x) <- 'summary.survregVB'
  x
}
