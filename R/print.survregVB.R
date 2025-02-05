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
  cat("\nMu:\n")
  print(x$mu, digits = digits)
  cat("\nSigma:\n")
  printCoefmat(x$Sigma, digits = digits, signif.stars=signif.stars,
               P.values=TRUE, has.Pvalue=TRUE)

  cat("\nPosterior distribution of b (scale parameter):")
  cat("\nAlpha: ", x$alpha, "  Omega: ", x$omega, "\n")

  omit <- x$na.action
  if (length(omit))
    cat("\nn=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("\nn=", x$n, "\n")

  invisible(x)
}
