#' Produces a summary table containing the value and confidence intervals
#' for an Inverse-Gamma posterior distribution.
#'
#' @param shape The shape parameter of the distribution.
#' @param scale The scale parameter of the distribution.
summary_table <- function(shape, scale) {
  summary_table <- c(expectation_b(shape, scale), b_ci(shape, scale))
  names(summary_table) <- c("Value", "Lower CI", "Upper CI")
  summary_table
}

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
#'    the approximate posterior distribution of the regression
#'    coefficients.
#'   \item \code{tau}: The updated parameter \eqn{\tau} of the approximate
#'    posterior distribution of the random intercept.
#'   \item \code{sigma}: The updated parameter \eqn{\sigma^2} of the approximate
#'    posterior distribution of the random intercept.
#'   \item \code{lambda}: The updated parameter \eqn{\lambda} of the approximate
#'    posterior distribution of the variance of the random intercept.
#'   \item \code{eta}: The updated parameter \eqn{\eta} of the approximate
#'    posterior distribution of the variance of the random intercept.
#'   \item \code{na.action}: A missing-data filter function, applied to
#'    the model.frame, after any subset argument has been used.
#'   \item \code{iterations}: The number of iterations performed by the
#'    VB algorithm: before converging or reaching \code{max_iteration}.
#'   \item \code{n}: The number of observations.
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
  scale <- summary_table(object$alpha, object$omega)

  # for models with shared frailty
  if (!is.null(object$clustered)) {
    intercept <- summary_table(object$lambda, object$eta)
    x <- object[match(c('clustered', 'call', 'ELBO', 'alpha', 'omega',
                        'mu', 'Sigma', 'tau', 'sigma', 'lambda', 'eta',
                        'na.action', 'iterations', 'n'),
                      names(object), nomatch=0)]
    x <- c(x, list(coefficients=coefficients, scale=scale, intercept=intercept))
  }
  else {
    x <- object[match(c('call', 'ELBO', 'alpha', 'omega', 'mu', 'Sigma',
                        'na.action', 'iterations', 'n'), names(object), nomatch=0)]
    x <- c(x, list(coefficients=coefficients, scale=scale))
  }

  class(x) <- 'summary.survregVB'
  x
}
