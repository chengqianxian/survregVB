---
title: "Functional Documentation"
editor_options: 
  markdown: 
    wrap: 72
---

#### `alpha_star`:

Calculates parameter $\alpha^*$ of the approximate posterior
distribution of *b* in order to optimize the evidence based lower bound
(ELBO) in `survregVB.fit`.

@param alpha_0 Hyperparameter $\alpha_0$ of the prior distribution of
*b*.

@param delta A binary vector indicating right censoring.

@returns Parameter $\alpha^*$ of the approximate posterior distribution
of *b*.

#### `omega_star`:

Calculates parameter $\omega^*$ of the approximate posterior
distribution of *b* in order to optimize the evidence based lower bound
(ELBO) in `survregVB.fit`.

@param y A vector of observed log-transformed survival times.

@param X A matrix of predictors (covariates), including an intercept.

@param delta A binary vector indicating right censoring.

@param omega_0 Hyperparameter $\omega_0$ of the prior distribution of
*b*.

@param mu Parameter $\mu^*$ of the approximate posterior distribution of
*β*.

@param expectation_b The expected value of b.

@returns Parameter $\omega^*$ of the approximate posterior distribution
of *b*.

#### `mu_star`:

Calculates parameter $\mu^*$ of the approximate posterior \#'
distribution of *β* in order to optimize the evidence based lower bound
(ELBO) in `survregVB.fit`.

@param y A vector of observed log-transformed survival times.

@param X A matrix of predictors (covariates), including an intercept.

@param delta A binary vector indicating right censoring.

@param alpha Parameter $\alpha^*$ of the approximate posterior \#'
distribution of *b*.

@param omega Parameter $\omega^*$ of the approximate posterior \#'
distribution of *b*.

@param curr_mu The current value of the parameter $\mu^*$ of the
approximate posterior distribution of *β*.

@param Sigma Parameter $\Sigma^*$ of the approximate posterior
distribution of *β*.

@param expectation_b The expected value of b.

@returns Parameter $\mu^*$ of the approximate posterior distribution of
*β*.

#### `Sigma_star`:

Calculates parameter $\Sigma^*$ of the approximate posterior \#'
distribution of *β* in order to optimize the evidence based lower bound
(ELBO) in `survregVB.fit`.

@param y A vector of observed log-transformed survival times.

@param X A matrix of predictors (covariates), including an intercept.

@param delta A binary vector indicating right censoring.

@param v_0 Hyperparameter $v_0$ of the prior distribution of *β*.

@param alpha Parameter $\alpha^*$ of the approximate posterior
distribution of *b*.

@param omega Parameter $\omega^*$ of the approximate posterior
distribution of *b*.

@param curr_mu The current value of the parameter $\mu^*$ of the
approximate posterior distribution of *β*.

@param expectation_b The expected value of b.

@returns Parameter $\Sigma^*$ of the approximate posterior \#'
distribution of *β*.

#### `expectation_log_likelihood:`

Calculates the approximated expectation over the density of observed
data *D\$ given parameters* β\* and *b*, $\log(p(D|\beta,b))$.

@param y A vector of observed log-transformed survival times.

@param X A matrix of predictors (covariates), including an intercept.

@param delta A binary vector indicating right censoring.

@param alpha Parameter $\alpha^*$ of the approximate posterior
distribution of *b*.

@param omega Parameter $\omega^*$ of the approximate posterior
distribution of *b*.

@param curr_mu The current value of the parameter $\mu^*$ of the
approximate posterior distribution of *β*.

@param expectation_b The expected value of b.

@returns The approximated log-likelihood $\log(p(D|\beta,b))$.

#### `diff_beta`

Calculates the approximated expectation over the density of observed
data *D* given parameters\* β\* and *b*, $\log(p(D|\beta,b))$.

@param v_0 Precision hyperparameter $v_0$ of the prior distribution of
*β*.

@param mu_0 Hyperparameter $\mu_0$ of the prior distribution of *β*.

@param curr_mu The current value of the parameter $\mu^*$ of the
approximate posterior distribution of *β*.

@param Sigma Parameter $\Sigma^*$ of the approximate posterior
distribution of *β*.

@returns The difference between the expectations of $\log(p(\beta))$ and
$\log(q(\beta))$.

#### `diff_b`:

Calculates the difference between the expectations of $\log(p(b))$ and
$\log(q(b))$.

@param alpha_0 Hyperparameter $\alpha_0$ of the prior distribution of
*b*.

@param omega_0 Hyperparameter $\omega_0$ of the prior distribution of
*b*.

@param alpha Parameter $\alpha^*$ of the approximate posterior
distribution of *b*.

@param omega Parameter $\omega^*$ of the approximate posterior
distribution of *b*.

@returns The difference between the expectations of \log(p(b)) and
\log(q(b)).

#### `elbo`:

Calculates the variational Bayes convergence criteria, evidence lower
bound (ELBO), optimized in `survregVB.fit`.

@param y A vector of observed log-transformed survival times.

@param X A matrix of predictors (covariates), including an intercept.

@param delta A binary vector indicating right censoring.

@param alpha_0 Hyperparameter $\alpha_0$ of the prior distribution of
*b*.

@param omega_0 Hyperparameter $\omega_0$ of the prior distribution of
*b*.

@param v_0 Precision hyperparameter $v_0$ of the prior distribution of
*β*.

@param mu_0 Hyperparameter $\mu_0$ of the prior distribution of *β*.

@param alpha Parameter $\alpha^*$ of the approximate posterior
distribution of *b*.

@param omega Parameter $\omega^*$ of the approximate posterior
distribution of *b*.

@param curr_mu The current value of the parameter $\mu^*$ of the
approximate posterior distribution of *β*.

@param Sigma Parameter $\Sigma^*$ of the approximate posterior
distribution of*β*.

@param expectation_b The expected value of b.

@returns The evidence lower bound (ELBO).

#### `beta_ci:`

Calculates the credible interval for the regression coefficient *β*.

@param mu The vector of means $\mu$ of the posterior distribution of
*β*.

@param Sigma The covariance matrix $\Sigma$ of the posterior
distribution *β*.

@param ci The significance level. (Default:0.95).

@returns Matrix containing the credible intervals for *β*.

#### `b_ci`:

Calculates the credible interval for the scale parameter *b* using the
highest density intervals (HDI) for b.

@param alpha The parameter $\alpha$ of the posterior distribution of
*b*.

@param omega The parameter $\omega$ of the posterior distribution *b*.

@param ci The significance level. (Default:0.95).

@returns The credible interval for *b*.

#### `survregVB.fit`:

Variational Bayesian Analysis of Survival Data Using a Log-Logistic
Accelerated Failure Time Model

Called by `survregVB` to do the actual parameter and ELBO computations.

@param Y A `Surv` object containing 2 columns: time and event.

@param X A matrix of predictors (covariates), including an intercept.

@param alpha_0 Hyperparameter $\alpha_0$ of the prior distribution of
*b*.

@param omega_0 Hyperparameter $\omega_0$ of the prior distribution of
*b*.

@param v_0 Precision hyperparameter $v_0$ of the prior distribution of
*β*.

@param mu_0 Hyperparameter $\mu_0$ of the prior distribution of *β*.

@param max_iteration The maximum number of iterations for the
variational inference optimization. If reached, iteration stops.
(Default:100)

@param threshold The convergence threshold for the evidence based lower
bound (ELBO) optimization. If the difference between the current and
previous ELBO's is smaller than this threshold, iteration stops.
(Default:0.0001)

@details This routine does no checking that the arguments are the proper
length or type.

@returns A list containing results of the fit.

#### `survregVB`:

Variational Bayesian Analysis of Survival Data Using a Log-Logistic
Accelerated Failure Time Model

Applies a mean-field Variational Bayes (VB) algorithm to infer the
parameters of an accelerated failure time (AFT) survival model with
right-censored survival times following a log-logistic distribution.

@param formula A formula object, with the response on the left of a `~`
operator, and the covariates on the right. The response must be a
survival object of type `right`, as returned by the `Surv` function.

@param data A `data.frame` in which to interpret the variables named in
the formula\`.

@param alpha_0 A numeric scalar specifying the shape hyperparameter
$\alpha_0$ of the prior Inverse-Gamma distribution for *b*.

@param omega_0 A numeric scalar specifying the scale hyperparameter
$\omega_0$ of the prior Inverse-Gamma distribution for *b*.

@param mu_0 A numeric vector containing the mean hyperparameters $\mu_0$
for the prior multivariate normal distributions of the intercept
($\beta_0$) and the coefficients ($\beta_i$) corresponding to the
covariates.

@param v_0 A numeric scalar specifying the precision (inverse variance)
hyperparameter $v_0$ of the prior multivariate normal prior distribution
for *β*.

@param na.action A missing-data filter function, applied to the
model.frame, after any subset argument has been used.
(Default:\\code{options()$na.action$).

@param max_iteration The maximum number of iterations for the
variational inference optimization. If reached, iteration stops.
(Default:100)

@param threshold The convergence threshold for the evidence based lower
bound (ELBO) optimization. If the difference between the current and
previous ELBO's is smaller than this threshold, iteration stops.
(Default:0.0001)

@returns An object of class `survregVB`. Objects of this class have
methods for the functions `print` and `summary`. The components of this
class are: - `ELBO`: The final value of the Evidence Lower Bound (ELBO)
at the last iteration.

-   `alpha`: The updated parameter $\alpha$ of the approximate posterior
    distribution of *b*.

-   `omega`: The updated parameter $\omega$ of the approximate posterior
    distribution of *b*.

-   `mu`: The updated vector of means $\mu$ of the approximate posterior
    distribution of *β*.

-   `Sigma`: The updated covariance matrix $\Sigma$ of the approximate
    posterior distribution of *β*.

-   `na.action`: A missing-data filter function, applied to the
    model.frame, after any subset argument has been used.

-   `iterations`: The number of iterations performed by the VB
    algorithm: before converging or reaching `max_iteration`.

-   `n`: The number of observations.

-   `call`: The function call used to invoke the `survregVB` method.

-   `not_converged`: A boolean indicating if the algorithm converged.

@details The log-logistic AFT model is specified with the parameters:
*β*, the vector of coefficients for the fixed effects, and *b*, a scale
parameter. The goal is to maximize the evidence lower bound (ELBO) in
order to approximate posterior distributions of the model parameters.

The prior distributions are:

-   $\beta\sim\text{MVN}(\mu_0,\sigma_0^2I_{p*p})$, with precision
    $v_0=1/\sigma^2$, and

-   $b\sim\text{Inverse-Gamma}(\alpha_0,\omega_0)$

The approximate posterior distributions are:

-   $q^*(\beta)$ is a $N_p(\mu,\Sigma)$ density function, and

-   $q^*(b)$ is an $\text{Inverse-Gamma}(\alpha,\omega)$ density
    function
