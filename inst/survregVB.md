---
title: "survregVB: Variational Bayesian Analysis of Survival Data Using a Log-Logistic Accelerated Failure Time Model"
author: Alison Zhang
output:
  html_document:
    citation_package: biblatex
bibliography: references.bib
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/apa.csl
link-citations: true
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

## Background

A commonly used survival regression model is the accelerated failure
time (AFT) model, which assumes an accelerative effect of the covariates
directly on survival time [@webber2022]. In particular, the log-logistic
distribution is suitable for modelling a wide variety of survival data
[@rivas-lópez2022]. In cases where correlated survival data arises from
clusters of individuals with shared environmental factors, a shared
frailty AFT model can be used to account for correlations among survival
data [@hougaard1995; @hanagal2011; @gorfine2023].

Bayesian inference is a technique used to derive the posterior
distribution of model parameters based on Bayes’ theorem. Markov Chain
Monte Carlo (MCMC) algorithms are typically used to approximate the
posterior [@geman1984]. Due to high computational costs, research has
been done to explore alternative methods such as variational inference
(VI) developed from machine learning [@jordan1999]. VI uses optimization
to approximate the parameters of Bayesian models, providing similar
estimations to MCMC techniques at a much lower computational cost
[@blei2017].

Mean-field variational Bayes (VB) is a case of mean-field VI that arises
from minimizing Kullback-Leibler (KL) divergence, the dissimilarity
between the approximated and exact posterior densities [@bishop2006].
Minimizing KL divergence can also be thought of as maximizing the
evidence lower-bound (ELBO) [@jordan1999; @blei2017].

## The survregVB Package

The **survregVB** package is an R package for the implementation of
Bayesian inference for accelerated failure time (AFT) models used in
survival analysis. As an alternative to Markov chain Monte Carlo (MCMC)
methods, the **survregVB** package implements a mean-field variational
Bayes (VB) algorithm to infer the parameters of a log-logistic AFT
model. A piecewise approximation technique is embedded into the VB
algorithm to achieve conjugacy. Compared to MCMC, the VB algorithm has
similar performance and significantly reduced computation cost, with an
average speed-up of 300 times [@xian2024]. The survregVB package
includes two methods, `survregVB` and `survregVB_frailty` for analyzing
AFT survival data under a log-logistic distribution with and without
frailty respectively.

The goal of the `survregVB` function is to fit it a log-logistic AFT
regression model of the following form for the $i^{th}$ subject in the
sample, $i=1,...,n$ , $T_i$:

$$
log(T_i):=Y=X_i^T\beta+bz_i
$$

where $X_i$ is column vector of $p-1$ covariates and a constant one
(i.e. $X_i=(1,x_i1,...,x_i(p-1))^T$), $\beta$ is a vector of
coefficients, $z_i$ is a random variable following a standard logistic
distribution and $b$ is a scale parameter. In particular, `survregVB`
uses a Bayesian framework to obtain the optimal variational densities of
parameters $\beta$ and $b$ by maximizing the evidence based lower bound
(ELBO) [@jordan1999; @blei2017]. In order to do so, we assume prior
distributions:

-   $\beta\sim\text{MVN}(\mu_0,\sigma_0^2I_{p*p})$, with precision
    $v_0=1/\sigma^2$, and

-   $b\sim\text{Inverse-Gamma}(\alpha_0,\omega_0)$

where $\mu_0,v_0,\alpha_0$ and $\omega_0$ are known hyperparameters
[@nedelman2005; @faes2011].

At the end of the model fitting process, the `survregVB` function
obtains the approximated posterior distributions given as follows:

-   $q^*(\beta)$ is a $N_p(\mu,\Sigma)$ density function, and

-   $q^*(b)$ is an $\text{Inverse-Gamma}(\alpha,\omega)$ density
    function,

where the parameters $\mu,\Sigma,\alpha$ and $\omega$ are obtained via
the VB algorithm [@xian2024].

We can also use the `survregVB` function is to fit it a shared frailty
log-logistic AFT regression model that accounts for intra-cluster
correlation through a cluster-specific random intercept. For time
$T_{ij}$ of the $j_{th}$ subject from the $i_{th}$ cluster in the
sample, in a sample with $i=1,...,K$ clusters and $j=1,...,n_i$
subjects:

$$
\log(T_{ij})=\gamma_i+X_{ij}^T\beta+b\epsilon_{ij}
$$

where $X_{ij}$ is column vector of $p-1$ covariates and a constant one
(i.e. $X_{ij}=(1,x_{ij1},...,x_{ij(p-1)})^T$), $\beta$ is a vector of
coefficients, $\gamma_i$ is a random intercept for the $i^{th}$ cluster,
$\epsilon_{ij}$ is a variable following a standard logistic distribution
and $b$ is a scale parameter [@xian2024a].

To obtain the optimal variational densities of parameters
$\beta,b,\gamma$ and $\sigma^2_\gamma$, we assume prior distributions:

-   $\beta\sim\text{MVN}(\mu_0,\sigma_0^2I_{p*p})$, with precision
    $v_0=1/\sigma^2$,

-   $\gamma_i|\sigma^2_\gamma\mathop{\sim}\limits^{\mathrm{iid}}N(0,\sigma^2_\gamma)$,

-   $b\sim\text{Inverse-Gamma}(\alpha_0,\omega_0)$, and

-   $\sigma^2\sim\text{Inverse-Gamma}(\lambda_0,\eta_0)$,

where $\mu_0,v_0,\alpha_0$, $\omega_0$, $\lambda_0$ and $\eta_0$ are
known hyperparameters [@xian2024a].

At the end of the model fitting process, the `survregVB` function
obtains the approximated posterior distributions given as follows:

-   $q^*(\beta)$ is a $N_p(\mu,\Sigma)$ density function,

-   $q^*(b)$ is a $\text{Inverse-Gamma}(\alpha,\omega)$ density
    function,

-   $q^*(\gamma_i)$ is a $N_l(\tau^*_i,\sigma^{2*}_i))$ density
    function, and

-   $q^*(\sigma^2_\gamma)$ is a $\text{Inverse-Gamma}(\lambda^*,\eta^*)$
    density function,

where the parameters
$\mu,\Sigma,\alpha,\omega,\tau_i, \sigma_i, \lambda$ and $\eta$ are
obtained via the VB algorithm [@xian2024a].

# Statement of Need

Several R packages, such as the likelihood-based survival regression
`survreg` function from the **survival** package [@survival] and the
MCMC-based Hamiltonian Monte Carlo (HMC) algorithm in the **rstan**
package [@rstan], provide methods for AFT model estimation. However,
there is currently no integrated software solution for implementing
variational inference (VI) for AFT models.

In comparison to MCMC techniques, VI has several advantages. These
include computational efficiency [@blei2017], the ability to make use of
prior information from similar studies, and the ability to conduct
inference for small sample sizes as it does not rely on asymptotics
[@ibrahim2001]. Research has been done to compare estimation results for
a log-logistic AFT model with right censoring from the VB algorithm to
those from the HMC algorithm, and `survreg`. The three methods provide
similar estimates of the posterior model parameters. However, the VB
method offers a significantly reduced computational cost compared to
MCMC, with an average speedup of 300 times. [@xian2024].

Bayesian variational methods for accelerated failure time (AFT) models
are still relatively new, and to date, there is no integrated software
solution for their practical implementation. Developing a practical tool
to implement the complex calculations involved in the algorithm would
enable researchers to more effectively apply AFT models in a variety of
real-world scenarios, particularly those involving large datasets or
complex survival analysis problems. In addition, The **survregVB**
package will perform additional steps such as removing missing values,
one-hot-encoding, and supplying summary statistics, to create a smoother
workflow.

# Example

The following example demonstrates how to use the `survregVB` function
to fit data from the `rhDNase` dataset from the **survival** package.
The `rhDNase` dataset contains the results of a randomized trial of
rhDNase for the treatment of cystic fibrosis. The survival time, *T*, is
defined as the time until the first pulmonary exacerbation. *T* is
calculated as the difference between the date of entry into the study
(`entry.dt`) and date of last follow-up (`end.dt`), with the follow-up
period capped at 169 days [@survival). The covariates of interest are:

-   Treatment (`trt`: 0=placebo, 1=rhDNase), and
-   Forced expiratory volume (`fev`) [@survival]

Our goal is to fit it a log-logistic AFT regression model of the form:

$${\log(T):=Y=\beta_0+\beta_1x_1+\beta_2x_2+bz}$$ where:

-   $x_1$ is `trt`,

-   $x_2$ is `fev`, and

-   $z$ is a random variable following a standard logistic distribution
    with scale parameter *b* [@xian2024].

The event of interest is whether or not the subject experienced
infection, `infect`.

## Preparing the Data

First, we prepare the data for modeling so the dataset is ready with the
survival time (`time`) and event indicator (`infect`) via the example
code provided in the survival package [@survival]:

``` r
library(survival)
# Extract the first row for each subject
first <- subset(rhDNase, !duplicated(id)) 
dnase <- tmerge(first, first, id=id, tstop=as.numeric(end.dt -entry.dt))
# Subject whose fu ended during the 6 day window after finishing their 
# antibiotics are not counted as new infections
temp.end <- with(rhDNase, pmin(ivstop+6, end.dt-entry.dt))
# Create the event indicator
dnase <- tmerge(dnase, rhDNase, id=id,
                infect=event(ivstart),
                end=  event(temp.end))
# Toss out the non-at-risk intervals, and extra variables
dnase <- subset(dnase, (infect==1 | end==0), c(id:trt, fev:infect))
# Remove duplicated subjects
dnase <- subset(dnase, !duplicated(id))
# Set the survival time 
dnase$time <- dnase$tstop - dnase$tstart
```

## Fitting the model

### Setting prior distributions

Next we define the prior distributions of the model parameters, *β* and
*b*.

We set hyperparameters $\alpha_0=501$ and $\omega_0=500$ to achieve a
mean scale of one.

The prior means for $\mu_0$ are chosen based on historical data and
similar analyses on this type of data:

-   for the intercept $\beta_0$, we use $\log(169/2)\approx4.4$ (half
    the follow-up period),

-   for $\beta_1$ (`trt`), we use $\log(1.28)\approx0.25$ [@shah1996]

-   for $\beta_2$ (`fev`) we use $\log(1.04)\approx0.04$ [@block2006]

For the precision hyperparameter, we use a low precision $v_0=1$ to
obtain a flat prior [@xian2024).

### Applying the `survregVB` function

Now that we have set up the prior distributions, we can fit the
`rhDNase` dataset to the log-logistic AFT model using the `survregVB`
function. At the end of the model fitting process, the VB algorithm
obtains the approximated posterior distributions of *b* and *β:*

``` r
fit <- survregVB(formula = Surv(time, infect) ~ trt + fev, 
                 data = dnase, 
                 alpha_0 = 501,
                 omega_0 = 500,
                 mu_0 = c(4.4, 0.25, 0.04), 
                 v_0 = 1, 
                 max_iteration = 10000, 
                 threshold = 0.0005)
```

We can view summary statistics for the results of the fit via the
`summary` function. We can adjust the significance level for the
credible intervals by altering the `ci` argument:

``` r
summary(fit, ci=0.9)
```

# Code Availability

The **survregVB** is available at the following GitHub repository:
<https://github.com/chengqianxian/survregVB>, along with tutorials as
vignettes in R Markdown notebooks.

# References
