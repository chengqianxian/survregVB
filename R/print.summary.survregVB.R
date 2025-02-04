#' @method print summary.survregVB
#' @export
print.summary.survregVB <- function(x, digits =
                                      max(options()$digits - 4, 3),
                                    signif.stars=FALSE, ...) {
  if(is.null(digits))
    digits <- options()$digits
  cat("Call:\n")
  dput(x$call)

  cat("\nELBO: ", round(x$ELBO, 3), "\n")
  cat("\nNumber of iterations: ", x$iterations, "\n")

  cat("\nRegression Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits,
               signif.stars=signif.stars, P.values=TRUE,
               has.Pvalue=TRUE)

  cat("\nScale Parameter:\n")
  print(x$scale, digits = digits)
  invisible(NULL)
}
