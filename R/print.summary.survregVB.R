#' @method print summary.survregVB
#' @export
print.summary.survregVB <- function(x, digits =
                                      max(options()$digits - 4, 3),
                                    signif.stars=FALSE, ...) {
  if(is.null(digits))
    digits <- options()$digits
  cat("Call:\n")
  dput(x$call)

  cat("\nELBO: ", round(x$ELBO, digits), "\n")
  cat("\nNumber of iterations: ", x$iterations, "\n")

  cat("\nRegression Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits,
               signif.stars=signif.stars, P.values=TRUE,
               has.Pvalue=TRUE)

  cat("\nScale Parameter:\n")
  print(x$scale, digits = digits)

  if (!is.null(x$clustered)) {
    cat("\nVariance of the Random Intercept:\n")
    print(x$intercept, digits = digits)
  }

  omit <- x$na.action
  if (length(omit))
    cat("\nn=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("\nn=", x$n, "\n")

  invisible(NULL)
}
