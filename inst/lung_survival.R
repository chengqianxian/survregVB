library(rstan)
# data cleaning
library(survival)
# data preparation
lung <- na.omit(lung)
Y <- lung$time
n <- nrow(lung)
delta <- lung$status - 1

# Convert categorical variables to factors
lung$sex <- factor(lung$sex)
lung$ph.ecog <- factor(lung$ph.ecog)

# Create the model matrix (dummy encoding for Stage)
X <- model.matrix(~ age + sex + ph.ecog, data = lung)
N <- nrow(X)
M <- ncol(X)

real_data = list(y = Y, event = delta, x = X, N = N, M = M)

stancode_llaft <- "
functions {
  // Log survival function
  vector log_S(vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -log(1 + (t[i] / scale[i])^shape);
    }
    return log_S;
  }

  // Log hazard function
  vector log_h(vector t, real shape, vector scale) {
    vector[num_elements(t)] log_h;
    vector[num_elements(t)] ls = log_S(t, shape, scale);

    for (i in 1:num_elements(t)) {
      log_h[i] = log(shape) - shape * log(scale[i]) + (shape - 1) * log(t[i])
                 - 2 * log(1 + (t[i] / scale[i])^shape) - ls[i];
    }
    return log_h;
  }

  // Log likelihood function for right-censored data
  real surv_llogist_lpdf(vector t, vector d, real shape, vector scale) {
    vector[num_elements(t)] log_lik = d .* log_h(t, shape, scale) + log_S(t, shape, scale);
    return sum(log_lik);
  }
}

// Data block
data {
  int<lower=1> N;  // Number of observations
  vector<lower=0>[N] y;  // Observation times
  vector<lower=0, upper=1>[N] event;  // Censoring indicator (1 = observed, 0 = censored)
  int<lower=1> M;  // Number of covariates
  matrix[N, M] x;  // Covariate matrix
}

// Parameters block
parameters {
  vector[M] beta;  // Coefficients for linear predictor
  real<lower=0> sigma;  // Scale parameter (sigma = 1 / shape)
}

// Transformed parameters block
transformed parameters {
  vector[N] linpred = x * beta;
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = exp(linpred[i]);
  }
}

// Model block
model {
  sigma ~ inv_gamma(11, 10);  // Prior for sigma, we use IG(11, 10)
  beta ~ multi_normal(to_vector({6.7, 0, 0, 0, 0, 0}),  // Mean vector for beta
    diag_matrix(rep_vector(1, M))  // Identity covariance matrix
  );  // Multivariate normal prior for beta
  y ~ surv_llogist(event,1/sigma,mu);  // Likelihood function
}

// Generated quantities block
generated quantities {
  vector[N] y_rep;  // Posterior predictive values
  vector[N] log_lik;  // Log-likelihood values

  for (n in 1:N) {
    log_lik[n] = ((log(1 / sigma) - (1 / sigma) * (x[n,] * beta) + ((1 / sigma) - 1) * log(y[n])
                   - 2 * log(1 + (y[n] / exp(x[n,] * beta))^(1 / sigma)))
                   - (-log(1 + (y[n] / exp(x[n,] * beta))^(1 / sigma)))) * event[n]
                   + (-log(1 + (y[n] / exp(x[n,] * beta))^(1 / sigma)));
  }

  for (n in 1:N) {
    real u = uniform_rng(0,1);
    y_rep[n] = exp(x[n,] * beta) * (((1 - u)^(-1) - 1)^sigma);
  }
}
"
start_time <- Sys.time()
M3 <- stan(model_code = stancode_llaft, data = real_data,
           iter = 2000, chains = 4)
end_time <- Sys.time()
#print(M3, pars = c("beta", "sigma"), digits = 3, probs = c(0.025, 0.975))
time1 <- end_time - start_time

start <- Sys.time()
vb <- survregVB(formula = Surv(time, status == 2) ~ age + factor(sex) +
            factor(ph.ecog), data = lung, alpha_0 = 501, omega_0 = 500,
          mu_0 = c(6.7, 0, 0, 0, 0, 0), v_0 = 1, threshold = 0.01)
end <- Sys.time()
time2 <- end - start

surv <- survreg(formula = Surv(time, status == 2) ~ age + factor(sex) +
          factor(ph.ecog), data = lung, dist = "loglogistic")
