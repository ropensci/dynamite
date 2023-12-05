// modified version of get_code(gaussian_example_fit)
// Instead of normal prior for the random intercepts,
// here we use t-distribution.
data {
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int<lower=0> K; // total number of covariates across all channels
  array[T] matrix[N, K] X; // covariates as an array of N x K matrices
  row_vector[K] X_m; // Means of all covariates at first time point
  int<lower=1> D; // number of B-splines
  matrix[D, T] Bs; // B-spline basis matrix
  int<lower=0> M; // number group-level effects (including intercepts)
  // number of fixed, varying and random coefficients, and related indices
  int<lower=0> K_fixed_y;
  int<lower=0> K_varying_y;
  int<lower=0> K_random_y; // includes the potential random intercept
  int<lower=0> K_y; // K_fixed + K_varying
  array[K_fixed_y] int J_fixed_y;
  array[K_varying_y] int J_varying_y;
  array[K_y] int J_y; // fixed and varying
  array[K_fixed_y] int L_fixed_y;
  array[K_varying_y] int L_varying_y;
  // Parameters of vectorized priors
  matrix[K_fixed_y, 2] beta_prior_pars_y;
  matrix[K_varying_y, 2] delta_prior_pars_y;
  matrix[K_varying_y, 2] tau_prior_pars_y;
  matrix[K_random_y, 2] sigma_nu_prior_pars_y;
  matrix[N, T] y_y;
}
transformed data {
}
parameters {
  // Random group-level effects
  vector<lower=0>[M] sigma_nu; // standard deviations of random effects
  matrix[N, M] nu_raw;
  vector[K_fixed_y] beta_y; // Fixed coefficients
  matrix[K_varying_y, D] omega_y; // Spline coefficients
  vector<lower=0>[K_varying_y] tau_y; // SDs for the random walks
  real a_y; // Mean of the first time point
  row_vector[D - 1] omega_raw_alpha_y; // Coefficients for alpha
  real<lower=0> tau_alpha_y; // SD for the random walk
  real<lower=0> sigma_y; // SD of the normal distribution
  real<lower=2> df;
}
transformed parameters {
  vector[1] sigma_nu_y = sigma_nu[1:1];
  matrix[N, 1] nu_y = diag_post_multiply(nu_raw[, 1:1], sigma_nu_y);
  // Time-varying coefficients
  array[T] vector[K_varying_y] delta_y;
  // Time-varying intercept
  array[T] real alpha_y;
  // Spline coefficients
  real omega_alpha_1_y;
  row_vector[D] omega_alpha_y;
  for (t in 1:T) {
    delta_y[t] = omega_y * Bs[, t];
  }
  // Define the first alpha using mean a_y
  {
    vector[K_y] gamma__y;
    gamma__y[L_fixed_y] = beta_y;
    gamma__y[L_varying_y] = delta_y[1];
    omega_alpha_1_y = a_y - X_m[J_y] * gamma__y;
  }
  omega_alpha_y[1] = omega_alpha_1_y;
  omega_alpha_y[2:D] = omega_raw_alpha_y;
  for (t in 1:T) {
    alpha_y[t] = omega_alpha_y * Bs[, t];
  }
}
model {
  df ~ gamma(2, 0.1);
  to_vector(nu_raw) ~ student_t(df, 0, 1);
  sigma_nu_y ~ normal(sigma_nu_prior_pars_y[, 1], sigma_nu_prior_pars_y[, 2]);
  a_y ~ normal(1.5, 3.1);
  omega_raw_alpha_y[1] ~ normal(omega_alpha_1_y, tau_alpha_y);
  for (i in 2:(D - 1)) {
    omega_raw_alpha_y[i] ~ normal(omega_raw_alpha_y[i - 1], tau_alpha_y);
  }
  tau_alpha_y ~ normal(0, 3.1);
  beta_y ~ normal(beta_prior_pars_y[, 1], beta_prior_pars_y[, 2]);
  omega_y[, 1] ~ normal(delta_prior_pars_y[, 1], delta_prior_pars_y[, 2]);
  for (i in 2:D) {
    omega_y[, i] ~ normal(omega_y[, i - 1], tau_y);
  }
  tau_y ~ normal(tau_prior_pars_y[, 1], tau_prior_pars_y[, 2]);
  sigma_y ~ exponential(0.65);
  {
    real ll = 0.0;
    vector[K_y] gamma__y;
    gamma__y[L_fixed_y] = beta_y;
    for (t in 1:T) {
      vector[N] intercept_y = alpha_y[t] + nu_y[, 1];
      gamma__y[L_varying_y] = delta_y[t];
      ll += normal_id_glm_lupdf(y_y[, t] | X[t][, J_y], intercept_y, gamma__y, sigma_y);
    }
    target += ll;
  }
}
