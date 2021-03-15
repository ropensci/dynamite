//
// data {
//   int<lower=1> N; // number of sequences
//   int<lower=1> T; // number of time points
//   int<lower=1> S; // number of symbols
//   int<lower=0> y[T, N]; //observations, 0 is missing value
//   int<lower=0> K; // number of covariates
//   matrix[N, K] x[T];
// }
//
// parameters {
//   matrix[S, T] mu_raw;
//   matrix[K, S-1] beta0;
// //  matrix[K, S] beta;
//   vector<lower=0>[S] sigma_mu;
// }
//
// transformed parameters {
//   matrix[K, S] beta = append_col(beta0, rep_vector(0, S));
//   matrix[S, T] mu;
//   mu[, 1] = 5 * mu_raw[, 1];
//   for(t in 2:T) mu[, t] = mu[, t-1] + sigma_mu .* mu_raw[, t];
// }
// model {
//   to_vector(beta0)  ~ normal(0, 5);
//   //to_vector(beta)  ~ normal(0, 5);
//   sigma_mu ~ normal(0, 5);
//   mu[, 1] ~ normal(0, 5);
//
//   to_vector(mu_raw) ~ std_normal();
//   // for(t in 2:T) {
//   //   mu[, t] ~ normal(mu[, t-1], sigma_mu);
//   // }
//    for(t in 1:T) {
//      // note: categorical_logit_glm_lupmf is not available with current (15/2/2021) version of stanheaders/rstan
//       target += categorical_logit_glm_lupmf(y[t,] | x[t], mu[, t], beta);
//     }
// }
// generated quantities {
//   matrix[S, S] A[T];
//   for(t in 1:T) {
//     for(s in 1:S) {
//       A[t, s, ] = softmax(mu[, t] + beta[, s])';
//     }
//   }
// }
