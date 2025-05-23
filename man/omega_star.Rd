% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{omega_star}
\alias{omega_star}
\title{Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the
evidence based lower bound (ELBO) in \code{survregVB.fit}.}
\usage{
omega_star(y, X, delta, omega_0, mu, expectation_b)
}
\arguments{
\item{y}{A vector of observed log-transformed survival times.}

\item{X}{A design matrix including covariates with first column of ones
to represent the intercept.}

\item{delta}{A binary vector indicating right censoring.}

\item{omega_0}{The shape hyperparameter \eqn{\omega_0} of the prior
distribution of the scale parameter, \emph{b}.}

\item{mu}{Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}, a vector of means.}

\item{expectation_b}{The expected value of \emph{b}.}
}
\value{
Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
}
\description{
Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the
evidence based lower bound (ELBO) in \code{survregVB.fit}.
}
\seealso{
\code{\link{survregVB.fit}}
}
