% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simulation_frailty}
\alias{simulation_frailty}
\title{Simulated data incorporating shared frailty effects to model clustered
time-to-event data.}
\format{
A dataframe with 75 observations grouped into 15 clusters, each
with 5 individuals.
\describe{
\item{x1}{Continuous covariate from \code{N(1, 0.2^2)}}
\item{x2}{Binary covariate from \code{Bernoulli(0.5)}}
\item{Time}{True survival time}
\item{Time.15}{Observed survival time accounting for uniformly distributed
right censoring time from \code{uniform(0,u)}}
\item{delta}{Event indicator for uncensored data (always 1 in this
simulation.)}
\item{delta.15}{Event indicator after censoring (1 = event, 0 =
censored).}
\item{cluster}{Cluster ID (1–15), indicating group-level frailty}
. }
@references Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
& Tian, R. (2024). Fast variational bayesian inference for correlated
survival data: An application to invasive mechanical ventilation
duration analysis. https://doi.org/10.48550/ARXIV.2408.00177
}
\usage{
simulation_frailty
}
\description{
Simulated data incorporating shared frailty effects to model clustered
time-to-event data.
}
\keyword{datasets}
