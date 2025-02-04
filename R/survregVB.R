#' Variational Bayesian Analysis of Survival Data Using a Log-Logistic
#' Accelerated Failure Time Model
#'
#' Applies a mean-field Variational Bayes (VB) algorithm to infer the
#' parameters of an accelerated failure time (AFT) survival model with
#' right-censored  survival times following a log-logistic distribution.
#'
#' @name survregVB
#' @aliases print.survregVB
#' @aliases survregVB.object
#'
#' @param formula A formula object, with the response on the left of a
#'  `~` operator, and the covariates on the right. The response must be
#'  a survival object of type `right`, as returned by the \code{Surv}
#'  function.
#' @param data A `data.frame` in which to interpret the variables named
#'  in the formula`.
#' @param alpha_0 A numeric scalar specifying the shape hyperparameter
#'  \eqn{\alpha_0} of the prior Inverse-Gamma distribution for \emph{b}.
#' @param omega_0 A numeric scalar specifying the scale hyperparameter
#'  \eqn{\omega_0} of the prior Inverse-Gamma distribution for \emph{b}.
#' @param mu_0 A numeric vector containing the mean hyperparameters
#'  \eqn{\mu_0} for the prior multivariate normal distributions of the
#'  intercept (\eqn{\beta_0}) and the coefficients (\eqn{\beta_i})
#'  corresponding to the covariates.
#' @param v_0 A numeric scalar specifying the precision (inverse
#'  variance) hyperparameter \eqn{v_0} of the prior multivariate normal
#'  prior distribution for \emph{β}.
#' @param max_iteration The maximum number of iterations for the
#'  variational inference optimization. If reached, iteration stops.
#'  (Default:100)
#' @param threshold The convergence threshold for the evidence based
#'  lower bound (ELBO) optimization. If the difference between the
#'  current and previous ELBO's is smaller than this threshold,
#'  iteration stops. (Default:0.0001)
#'
#' @returns An object of class \code{survregVB}.
#' Objects of this class have methods for the functions \code{print} and
#' \code{summary}. The components of this class are:
#' \itemize{
#'   \item \code{ELBO}: The final value of the Evidence Lower Bound
#'    (ELBO) at the last iteration.
#'   \item \code{alpha}: The updated parameter \eqn{\alpha} of the
#'    approximate posterior distribution of \emph{b}.
#'   \item \code{omega}: The updated parameter \eqn{\omega} of the
#'    approximate posterior distribution of \emph{b}.
#'   \item \code{mu}: The updated vector of means \eqn{\mu} of the
#'    approximate posterior distribution of \emph{β}.
#'   \item \code{Sigma}: The updated covariance matrix \eqn{\Sigma} of
#'    the approximate posterior distribution of \emph{β}.
#'   \item \code{iterations}: The number of iterations performed by the
#'    VB algorithm: before converging or reaching \code{max_iteration}.
#'   \item \code{call}: The function call used to invoke the
#'    \code{survregVB} method.
#'   \item \code{not_converged}: A boolean indicating if the algorithm
#'    converged.
#'   \itemize{
#'     \item \code{TRUE}: If the algorithm did not converge prior to
#'      achieving `max_iteration`.
#'     \item \code{NULL}: If the algorithm converged successfully.
#'   }
#' }
#'
#' @details
#' The log-logistic AFT model is specified with the parameters:
#' \emph{β}, the vector of coefficients for the fixed effects, and
#' \emph{b}, a scale parameter.
#' The goal is to maximize the evidence lower bound (ELBO) in order to
#' approximate posterior distributions of the model parameters.
#'
#' The prior distributions are:
#' - \eqn{\beta\sim\text{MVN}(\mu_{0},\sigma_{0}^2I_{p*p})} where
#' precision, \eqn{v_{0}=1/\sigma^2}, and
#' - \eqn{b\sim\text{Inverse-Gamma}(\alpha_0,\omega_0)}.
#'
#' The approximate posterior distributions are:
#' - \eqn{\beta\sim N_p(\mu,\Sigma)}, and
#' - \eqn{q^*(b)\sim\text{Inverse-Gamma}(\alpha,\omega)}.
#'
#' @examples

#' # Data frame containing survival data
#' example_data <- data.frame(
#'   time = c(100, 120, 90, 35, 140),
#'   status = c(1, 0, 1, 0, 1),
#'   age = c(60, 55, 45, 70, 50),
#'   gender = factor(c("M", "F", "M", "M", "F"))
#' )
#' # Formula for the survival model
#' example_formula <- survival::Surv(time, status) ~ age + gender
#' # Call the survregVB function
#' result <- survregVB(formula = example_formula,
#'                     data = example_data,
#'                     alpha_0 = 10,
#'                     omega_0 = 11,
#'                     mu_0 = c(5,0,0),
#'                     v_0 = 0.2)
#' # View results
#' summary(result)
#' @import stats
#'
#' @export
survregVB <- function (formula, data, alpha_0, omega_0, mu_0, v_0,
                       max_iteration = 100, threshold = 0.0001) {
  Call <- match.call()    # save a copy of the call

  if (missing(formula)) stop("a formula argument is required")
  defined <- ls(pattern = '.*_0$')
  passed <- names(as.list(Call)[-1])
  if (any(!defined %in% passed))
    stop(paste("Missing value(s) for", paste(setdiff(defined, passed),
                                             collapse =
                                               ", ")))

  if (is.list(formula))
    stop("formula argument cannot be a list")
  else
    Terms <- if (missing(data))
      terms(formula)
  else
    terms(formula, data = data)

  indx <- match(c("formula", "data"), names(Call), nomatch=0)
  if (indx[1] ==0) stop("A formula argument is required")
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame)   # change the function called

  temp$formula <- if (missing(data))
    terms(formula)
  else
    terms(formula, data = data)
  m <- eval(temp, parent.frame())

  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")

  if (!inherits(Y, "Surv"))
    stop("response must be a survival object")
  type <- attr(Y, "type")
  if (type != "right")
    stop("only Survival objects of type right are supported")

  X <- model.matrix(Terms, m)
  if (!all(is.finite(X)) || !all(is.finite(Y)))
    stop("data contains an infinite predictor")

  result <- survregVB.fit(Y, X, alpha_0, omega_0, mu_0, v_0,
                          max_iteration, threshold)
  result$call <- Call

  class(result) <- 'survregVB'
  result
}
