#' Variational Bayesian Analysis of Survival Data Using a Log-Logistic
#' Accelerated Failure Time Model
#'
#' Applies a mean-field Variational Bayes (VB) algorithm to infer the
#' parameters of an accelerated failure time (AFT) survival model with
#' right-censored  survival times following a log-logistic distribution.
#'
#' @name survregVB
#'
#' @param formula A formula object, with the response on the left of a `~`
#'  operator, and the covariates on the right. The response must be a survival
#'  object of type `right`, as returned by the \code{Surv} function.
#' @param data A `data.frame` in which to interpret the variables named in
#'  the formula`.
#' @param alpha_0 The shape hyperparameter \eqn{\alpha_0} of the prior
#'  distribution of \emph{b}.
#' @param omega_0 The shape hyperparameter \eqn{\omega_0} of the prior
#'  distribution of \emph{b}.
#' @param mu_0 Hyperparameter \eqn{\mu_0}, a vector of means, of the prior
#'  distribution of \emph{β}.
#' @param v_0 The precision (inverse variance) hyperparameter \eqn{v_0},
#'  of the prior distribution of \emph{β}.
#' @param lambda_0 The shape hyperparameter \eqn{\lambda_0} of the prior
#'  distribution of \eqn{\sigma_\gamma^2}.
#' @param eta_0 The scale hyperparameter \eqn{\eta_0} of the prior distribution
#'  of \eqn{\sigma_\gamma^2}.
#' @param na.action A missing-data filter function, applied to the
#'  \code{model.frame}, after any subset argument has been used.
#'  (Default:\code{options()$na.action}).
#' @param cluster An optional variable which clusters the observations to
#'  introduce shared frailty for correlated survival data.
#' @param max_iteration The maximum number of iterations for the variational
#'  inference optimization. If reached, iteration stops. (Default:100)
#' @param threshold The convergence threshold for the evidence based lower
#'  bound (ELBO) optimization. If the difference between the current and
#'  previous ELBO's is smaller than this threshold, iteration stops.
#'  (Default:0.0001)
#'
#' @returns An object of class \code{survregVB}.
#'
#' @details
#' The goal of \code{survregVB} is to maximize the evidence lower bound
#' (ELBO) to approximate posterior distributions of the model parameters.
#'
#' The log-logistic AFT model without shared frailty is specified with the
#' parameters:
#' - \emph{β}, the vector of coefficients for the fixed effects, and
#' - \emph{b}, a scale parameter.
#'
#' We assume prior distributions:
#' - \eqn{\beta\sim\text{MVN}(\mu_{0},\sigma_{0}^2I_{p*p})} with precision
#'   \eqn{v_{0}=1/\sigma^2}, and
#' - \eqn{b\sim\text{Inverse-Gamma}(\alpha_0,\omega_0)},
#'
#' and obtain approximate posterior distributions:
#' - \eqn{q^*(\beta)}, a \eqn{N_p(\mu^*,\Sigma^*)} density function, and
#' - \eqn{q^*(b)}, an \eqn{\text{Inverse-Gamma}(\alpha^*,\omega^*)}
#'    density function.
#'
#' With shared frailty, a model with \eqn{i=1,...,K} clusters is specified
#' with additional parameters:
#' - \eqn{\sigma^2_\gamma}, the random intercept, and
#' - \eqn{\gamma_i|\sigma^2_\gamma} the random effects.
#'
#' We additionally assume prior distributions:
#' - \eqn{\sigma^2_\gamma\sim\text{Inverse-Gamma}(\lambda_0,\eta_0)}, and
#' - \eqn{\gamma_i|\sigma^2_\gamma\mathop{\sim}\limits^{\mathrm{iid}}
#'    N(0,\sigma^2_\gamma)}.
#'
#' and obtain posterior distributions:
#' - \eqn{q^*(\sigma^2_\gamma)}, an \eqn{\text{Inverse-Gamma}(\lambda^*,\eta^*)}
#'   density function, and
#' - \eqn{q^*(\gamma_i)}, a \eqn{N_l(\tau^*_i,\sigma^{2*}_i))} density function.
#'
#' @examples
#' # Data frame containing survival data
#' example_data <- data.frame(
#'   time = c(100, 120, 90, 35, 140),
#'   status = c(1, 0, 1, 0, 1),
#'   age = c(60, 55, 45, 70, 50),
#'   gender = c("M", "F", "M", "M", "F"),
#'   group = c(1, 1, 2, 3, 3)
#' )
#' # Formula for the survival model
#' example_formula <- survival::Surv(time, status) ~ age + gender
#' # Call the survregVB function
#' result1 <- survregVB(
#'   formula = example_formula,
#'   data = example_data,
#'   alpha_0 = 10,
#'   omega_0 = 11,
#'   mu_0 = c(5, 0, 0),
#'   v_0 = 0.2
#' )
#' # View results
#' summary(result1)
#'
#' # Call the survregVB function with shared frailty
#' result2 <- survregVB(
#'   formula = example_formula,
#'   data = example_data,
#'   alpha_0 = 10,
#'   omega_0 = 11,
#'   mu_0 = c(5, 0, 0),
#'   v_0 = 0.2,
#'   lambda_0 = 3,
#'   eta_0 = 2,
#'   cluster = group
#' )
#' @import stats
#'
#' @export
#' @seealso \code{\link{survregVB.object}}
survregVB <- function(formula, data, alpha_0, omega_0, mu_0, v_0,
                      lambda_0, eta_0, na.action, cluster,
                      max_iteration = 100, threshold = 0.0001) {
  Call <- match.call() # save a copy of the call

  if (missing(formula)) stop("a formula argument is required")
  if (is.list(formula)) {
    stop("formula argument cannot be a list")
  }

  Terms <- if (missing(data)) {
    terms(formula)
  } else {
    terms(formula, data = data)
  }

  defined <- if (missing(cluster)) {
    c("alpha_0", "omega_0", "mu_0", "v_0")
  } else {
    c("alpha_0", "omega_0", "mu_0", "v_0", "lambda_0", "eta_0")
  }
  passed <- match(defined, names(Call), nomatch = 0)
  missing <- which(passed == 0)
  if (length(missing) > 0) {
    stop(paste("missing value(s) for", defined[missing]))
  }

  indx <- match(c("formula", "data", "cluster", "na.action"),
    names(Call),
    nomatch = 0
  )
  if (indx[1] == 0) stop("A formula argument is required")

  temp <- Call[c(1, indx)] # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame) # change the function called

  temp$formula <- if (missing(data)) {
    terms(formula)
  } else {
    terms(formula, data = data)
  }

  m <- eval(temp, parent.frame())

  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")

  if (!inherits(Y, "Surv")) {
    stop("response must be a survival object")
  }
  type <- attr(Y, "type")
  if (type != "right") {
    stop("only Survival objects of type right are supported")
  }

  X <- model.matrix(Terms, m)
  if (!all(is.finite(X)) || !all(is.finite(Y))) {
    stop("data contains an infinite predictor")
  }

  cluster <- model.extract(m, "cluster")

  if (length(cluster)) {
    result <- survregVB.frailty.fit(
      Y, X, alpha_0, omega_0, mu_0, v_0,
      lambda_0, eta_0, cluster,
      max_iteration, threshold
    )
  } else {
    result <- survregVB.fit(
      Y, X, alpha_0, omega_0, mu_0, v_0, max_iteration,
      threshold
    )
  }

  na.action <- attr(m, "na.action")
  if (length(na.action)) result$na.action <- na.action
  result$call <- Call
  class(result) <- "survregVB"

  result
}

#' Variational Bayes Accelererated Failure Time Survival Model Object
#'
#' This class of objects is returned by the survregVB function to represent
#' a fitted parametric log-logistic accelerated failure time (AFT) survival
#' model. Objects of this class have methods for the functions \code{print}
#' and \code{summary}.
#'
#' @name survregVB.object
#' @aliases print.survregVB
#'
#' @details
#' For approximate posterior distributions:
#' - \eqn{q^*(\beta)}, a \eqn{N_p(\mu^*,\Sigma^*)} density function, and
#' - \eqn{q^*(b)}, an \eqn{\text{Inverse-Gamma}(\alpha^*,\omega^*)}
#'    density function,
#'
#' the components of this class are:
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
#' }
#'
#' If \code{survregVB} was called with shared frailty (with the `cluster`
#' argument), for approximate posterior distributions:
#' - \eqn{q^*(\sigma^2_\gamma)}, an \eqn{\text{Inverse-Gamma}(\lambda^*,\eta^*)}
#'   density function,
#' - \eqn{q^*(\gamma_i)}, a \eqn{N_l(\tau^*_i,\sigma^{2*}_i))} density function,
#'   for \eqn{i=1,...,K} clusters, and
#'
#' the additional components are present:
#'
#' \itemize{
#'   \item \code{lambda}: The shape parameter \eqn{\lambda^*} of
#'    \eqn{q^*(\sigma^2_\gamma)}.
#'   \item \code{eta}: The scale parameter \eqn{\eta^*} of
#'    \eqn{q^*(\sigma^2_\gamma)}.
#'   \item \code{tau}: Parameter \eqn{\tau^*} of \eqn{q^*(\gamma_i)}, a
#'    vector of means.
#'   \item \code{sigma}: Parameter \eqn{\sigma^{2*}_i} of \eqn{q^*(\gamma_i)},
#'    a vector of variance.
#'  }
#'
NULL
