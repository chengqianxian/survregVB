#' @method print survregVB
#' @export
print.survregVB <- function(x, digits = max(options()$digits - 4, 3),
                            signif.stars=FALSE, ...) {
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  if (!is.null(x$not_converged)) {
    cat("\nThe VB algorithm did not converge.\n")
  }

  cat("\nELBO: ", round(x$ELBO, digits), "\n")
  cat("\nNumber of iterations: ", x$iterations, "\n")

  cat("\nPosterior distributions of Beta (regression coefficients):")
  cat("\nmu:\n")
  print(x$mu, digits = digits)
  cat("\nSigma:\n")
  printCoefmat(x$Sigma, digits = digits, signif.stars=signif.stars,
               P.values=TRUE, has.Pvalue=TRUE)

  cat("\nPosterior distribution of b (scale parameter):")
  cat("\nalpha: ", x$alpha, "  omega: ", x$omega, "\n")

  if (!is.null(x$clustered)) {
    cat("\nPosterior distribution of gamma (random intercept):")
    cat("\ntau:\n")
    print(x$tau, digits = digits)
    cat("\nsigma squared:\n")
    print(x$sigma, digits = digits)

    cat("\nPosterior distribution of sigma_gamma squared (variance of the random intercept):")
    cat("\nlamda: ", x$lambda, "  eta: ", x$eta, "\n")
  }

  omit <- x$na.action
  if (length(omit))
    cat("\nn=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("\nn=", x$n, "\n")

  invisible(x)
}
