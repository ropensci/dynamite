// Model with time-constant effects and correlated id-specific noise terms
// Note: We could fit this model easily with brms as well.
// something like:
// bf1 <- bf(y1 ~ y1lag + y2lag + y3lag + (1|p|id), family = categorical)
// bf2 <- bf(y2 ~ y1lag + y2lag + y3lag + (1|p|id), family = categorical)
// bf3 <- bf(y3 ~ y1lag + y2lag + y3lag + (1|p|id), family = gaussian)
// options(contrasts = rep("contr.sum", 2))
// fit <- brm(bf1 + bf2 + bf3,data=df[-1,], chains = 1)

functions {

  /* Create an array of vectors from one long vector
  * alpha[i] contains the intercept terms for the sequence i
  * Each alpha[i] has length S_max but only first S[i] are used
  * (with last one being zero for identifiability purposes)
  */
  vector[] alpha_vector(vector alphav, int Pc, int[] S) {
    // intercept terms for the sequences
    // alphav:
    // 1:(S[1] - 1) for the first series
    // S[1]:(S[1] + S[2] - 2) for the second
    // (S[1] + S[2] - 1):(S[1] + S[2] + S[3] - 3) for the third etc
    int S_max = max(S);
    vector[S_max] alpha[Pc];
    int k = 0;
    for(i in 1:Pc) {
      alpha[i] = rep_vector(0, S_max);
      alpha[i, 1:(S[i]-1)] = alphav[(k + 1):(k + S[i] - 1)];
      k += S[i] - 1;
    }
    return alpha;
  }

  /* Create array of matrices of beta coefficients for sequences  */
  matrix[] beta_matrix(vector betav, int Pc, int[] S, int[] K) {
    // regression coefficients terms for sequences
    // betav:
    // 1:(K[1] * (S[1] - 1)) for the first series
    // (K[1] * (S[1] - 1) + 1):(K[1] * (S[1] - 1) + 1 + K[2] * (S[2] - 1)) for the second
    int S_max = max(S);
    int K_max = max(K);
    matrix[K_max, S_max] beta[Pc];
    int k = 0;
    for(i in 1:Pc) {
      beta[i] = rep_matrix(0, K_max, S_max);
      beta[i, 1:K[i], 1:(S[i] - 1)] = to_matrix(betav[(k + 1):(k + K[i] * (S[i] - 1))], K[i], S[i] - 1);
      k += K[i] * (S[i] - 1);
    }
    return beta;
  }
  /* Array of vectors of regression coefficients for non-sequence data */
  vector[] beta_vector(vector betav, int P, int[] K) {
    int K_max = max(K);
    vector[K_max] beta[P];
    int k = 0;
    for(i in 1:P) {
      beta[i] = rep_vector(0, K_max);
      beta[i, 1:K[i]] = betav[(k + 1):(k + K[i])];
      k += K[i];
    }
    return beta;
  }



  // from brms
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }


}
data {
  int<lower=1> N; // number of sequences
  int<lower=1> T; // number of time points

  // categorical observations
  int<lower=0> Pc; // number of categorical sequences
  int<lower=1> S[Pc]; // number of symbols per categorical series
  int<lower=1> yc[Pc, T, N]; // categorical sequences

  // Gaussian observations
  int<lower=0> Pg;
  real yg[Pg, T, N];


  //int<lower=0> Pp; // number of Poisson series
  //real<lower=0> exposure[Pp, T, N];
  //int<lower=0> yp[Pp, T, N];

  //int<lower=0> Pb; // number of Binomial series
  //real<lower=0> trials[Pb, T, N];
  //int<lower=0> yb[Pb, T, N];

  // covariates, can include variables above,
  // but we need to separate them at least because of int/real discrepancy
  int<lower=0> K; // total number of unique covariates
  matrix[N, K] x[T];
  int<lower=0, upper=K> Kj[Pc + Pg]; //number covariates per series
  // covariate indices for ith series (use only 1:Kj[i] values)
  int<lower=0,upper=K> Xidx[K, Pc + Pg];// Xidx[1:Kj[i], i]
}

transformed data {
  int S_max = 0;
  int S_sum = 0;
  int Kc_max = 0;
  int Kc_sum = 0;

  int Kg_max = 0;
  int Kg_sum = 0;

  if (Pc > 0) {
    S_max = max(S);
    S_sum = sum(S) - Pc;
    Kc_max = max(Kj[1:Pc]);
    Kc_sum = sum(Kj[1:Pc]);
  }
  if (Pg > 0) {
    Kg_max = max(Kj[(Pc + 1):(Pc + Pg)]);
    Kg_sum = sum(Kj[(Pc + 1):(Pc + Pg)]);
  }
}

parameters {
  // intercept terms for sequences
  // 1:(S[1] - 1) for first type
  // S[1]:(S[1] + S[2] - 1 for second etc)
  vector[S_sum] alpha_c0;
  vector[S_sum * Kc_sum] beta_c0;
  vector<lower=0>[Pg] sigma;
  vector[Pg] alpha_g;
  vector[Pg * Kg_sum] beta_g0;

  // correlated zero-mean random effect for each series
  // for categorical we need multiple terms
  // cholesky_factor_corr[Pg + S_sum] L;
  // vector<lower=0>[Pg + S_sum] sigma_re;
}

transformed parameters {
  vector[S_max] alpha_c[Pc] = alpha_vector(alpha_c0, Pc, S);
  matrix[Kc_max, S_max] beta_c[Pc] = beta_matrix(beta_c0, Pc, S, Kj[1:Pc]);
  vector[Kg_max] beta_g[Pg] = beta_vector(beta_g0, Pg, Kj[(Pc + 1):(Pc + Pg)]);
}
model {
  // hard-coded for simplicity here,
  // we can use alpha_c0 ~ normal(alpha_c_mean, alpha_c_sd)
  // where alpha_c_mean and alpha_c_sd are given by the user
  alpha_c0  ~ normal(0, 10);
  beta_c0  ~ normal(0, 10);
  alpha_g ~ normal(0, 10);
  beta_g0 ~ normal(0, 10);
  sigma ~ normal(0, 10);
  for(i in 1:Pc) {
    for(t in 2:T) {
      // note: categorical_logit_glm_lupmf is not available with current (15/2/2021) version of stanheaders/rstan
      // but should make this faster
      matrix[N, S[i]] xbeta = x[t, ,Xidx[1:Kj[i], i]] * beta_c[i, 1:Kj[i], 1:S[i]];
      for(ii in 1:N) {
        yc[i, t, ii] ~ categorical_logit(alpha_c[i, 1:S[i]] + xbeta[ii]');
      }
    }
  }
  for(i in 1:Pg) {
    for(t in 2:T) {
      yg[i, t, ] ~ normal(alpha_g[i] + x[t, , Xidx[1:Kj[i], i]] * beta_g[i, 1:Kj[i]], sigma[i]); //+ random_effect_g[i]
    }
  }
}
