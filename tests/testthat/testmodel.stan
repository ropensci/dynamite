// Model for testing predictions, see test-recovery.R
data {
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int y[N, T];
}

parameters {
  real beta;
  real alpha;
  vector[N] nu_raw;
  real<lower=0> sigma_nu;
}
transformed parameters {
  vector[N] nu = sigma_nu * nu_raw;
}

model {
  alpha ~ std_normal();
  beta ~ std_normal();
  nu_raw ~ std_normal();
  sigma_nu ~ std_normal();
  for(i in 1:N) {
    for(t in 2:T) {
      y[i, t] ~ bernoulli_logit(nu[i] + alpha + beta * y[i, t - 1]);
    }
  }
}
generated quantities {
  matrix[N, T] y_m;
  matrix[N, T] y_rep;
  vector[T] mean_y_m;
  vector[T] sd_y_m;
  vector[T] mean_y;
  vector[T] sd_y;
  for(i in 1:N) {
    y_m[i, 1] = y[i, 1];
    y_rep[i, 1] = y[i, 1];
    for(t in 2:T) {
      y_m[i, t] = inv_logit(nu[i] + alpha + beta * y_rep[i, t - 1]);
      y_rep[i, t] = bernoulli_rng(y_m[i, t]);
    }
  }
  for(t in 1:T) {
    mean_y[t] = mean(y_rep[, t]);
    sd_y[t] = sd(y_rep[, t]);
    mean_y_m[t] = mean(y_m[, t]);
    sd_y_m[t] = sd(y_m[, t]);
  }
}
