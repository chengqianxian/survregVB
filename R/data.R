## Datasets ============================================================

#' Subset of `rhDNase` from the `survival` package
#'
#' The `dnase` dataset is a subset of the `rhDNase` dataset from the
#' `survival` package.
#' It is included in this package under the LGPL (≥ 2) license.
#'
#' @format A data frame with 767 observations on the following variables:
#' \describe{
#'   \item{trt}{treatment arm: 0=placebo, 1= rhDNase}
#'   \item{fev}{forced expriatory volume at enrollment, a measure of lung
#'    capacity}
#'   \item{infect}{an infection that required the use of intravenous
#'    antibiotics}
#'   \item{time}{difference between the date of entry into the study and
#'    the date of last follow-up capped at 169 days}
#' }
#'
#' @source `survival` package.
#'  \url{https://cran.r-project.org/package=survival}
"dnase"

#' Simulated Survival Data Without Frailty
#'
#' This dataset is a simulated survival dataset generated to demonstrate
#' Bayesian inference in accelerated failure time (AFT) models.
#'
#' @format A data frame with 300 rows and the following variables:
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
#' @format A data frame with 75 rows and the following variables:
#' \describe{
#'   \item{x1}{Continuous covariate from `N(1, 0.2^2)`}
#'   \item{x2}{Binary covariate from `Bernoulli(0.5)`}
#'   \item{T}{True survival time}
#'   \item{T.15}{Observed survival time under `cen.time.15`}
#'   \item{delta}{Event indicator (1 = event, 0 = censored)}
#'   \item{delta.15}{Censoring indicator under `cen.time.15`}
#'   \item{cluster}{Cluster ID (1–15), indicating group-level frailty}
#' }
#' @source Simulated data
"simulation_frailty"

#' Subset of GSE102287: African American (AA) Patients
#'
#' This dataset is a subset of the GSE102287 dataset that includes only
#' characteristics of patients who are identified as African American (AA).
#'
#' @format A data frame with 60 observations on selected patient
#'  characteristics:
#' \describe{
#'   \item{patient}{Patient identification number.}
#'   \item{age}{Patient age.}
#'   \item{Stage}{Lung cancer stage (I, II, III).}
#'   \item{time}{Survival time in days.}
#'   \item{gender}{Gender of the patient.}
#'   \item{smoking}{0 = Never smoked, 1 = Has smoked.}
#'   \item{status}{0 = Alive, 1 = Death due to lung cancer.}
#' }
#'
#' @source  Gene Expression Omnibus (GEO), Accession: GSE102287.
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102287}
#'
#' @references Mitchell, K. A., Zingone, A., Toulabi, L., Boeckelman, J.,
#' & Ryan, B. M. (2017). Comparative Transcriptome Profiling Reveals Coding
#' and Noncoding RNA Differences in NSCLC from African Americans and European
#' Americans. Clinical cancer research: an official journal of the American
#' Association for Cancer Research, 23(23), 7412–7425.
#' doi:10.1158/1078-0432.CCR-17-0527.
"lung_cancer"
