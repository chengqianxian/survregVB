#' Subset of `rhDNase` from the `survival` package
#'
#' The `dnase` dataset is a subset of the `rhDNase` dataset from the
#' `survival` package.
#' It is included in this package under the LGPL (≥ 2) license.
#'
#' @format ## `dnase`
#' A data frame with 767 observations on the following variables:
#' \describe{
#'   \item{trt}{treatment arm: 0=placebo, 1= rhDNase}
#'   \item{fev}{forced expriatory volume at enrollment, a measure of lung
#'    capacity}
#'   \item{infect}{an infection that required the use of intravenous antibiotics}
#'   \item{time}{difference between the date of entry into the study and
#'    the date of last follow-up capped at 169 days}
#' }
#'
#' @source `survival` package, https://cran.r-project.org/package=survival
"dnase"

#' Simulated Survival Data Without Frailty
#'
#' This dataset is a simulated survival dataset generated to demonstrate
#' Bayesian inference in accelerated failure time (AFT) models.
#'
#' @format A data frame with 300 rows and 8 variables:
#' \describe{
#'   \item{x1}{Continuous covariate from `N(1, 0.2)`}
#'   \item{x2}{Binary covariate from `Bernoulli(0.5)`}
#'   \item{T}{True survival time}
#'   \item{T.10}{Observed time under `cen.time.10`}
#'   \item{T.30}{Observed time under `cen.time.30`}
#'   \item{delta}{Event indicator (1 = event, 0 = censored)}
#'   \item{delta.10}{Censoring indicator under `cen.time.10`}
#'   \item{delta.30}{Censoring indicator under `cen.time.30`}
#' }
#' @source Simulated data
"simulation_nofrailty"

#' Simulated Survival Data with Frailty
#'
#' This dataset is a simulated survival dataset incorporating frailty effects
#' to account for cluster-level uncertainty. It is designed for demonstrating
#' Bayesian inference in accelerated failure time (AFT) models.
#'
#' @format A data frame with 75 rows and 7 variables:
#' \describe{
#'   \item{x1}{Continuous covariate from `N(1, 0.2)`}
#'   \item{x2}{Binary covariate from `Bernoulli(0.5)`}
#'   \item{T}{True survival time}
#'   \item{T.15}{Observed survival time under `cen.time.15`}
#'   \item{delta}{Event indicator (1 = event, 0 = censored)}
#'   \item{delta.15}{Censoring indicator under `cen.time.15`}
#'   \item{cluster}{Cluster ID (1–15), indicating group-level frailty}
#' }
#' @source Simulated data
"simulation_frailty"
