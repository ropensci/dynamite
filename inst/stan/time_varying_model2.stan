// Model with time-varying effects
// creates coefficient matrices in generated quantities
// could be slightly faster than the other version

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

}

// awful names, need to come up with better naming conventions...
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
  int<lower=0, upper=K> Kt[Pc + Pg]; //number covariates with time varying effects per series

  // covariate indices for ith series (use only 1:Kj[i] values)
  int<lower=0,upper=K> Xidx[K, Pc + Pg];// Xidx[1:Kj[i], i]
  // covariate indices for time varying effects for ith series (use only 1:Kt[i] values)
  int<lower=0,upper=K> Xidx_t[K, Pc + Pg];
}

transformed data {
  int S_max = 0;
  int S_sum = 0;
  int Kc_max = 0;
  int Kc_sum = 0;
  int Kc_max_t = 0;
  int Kc_sum_t = 0;

  int Kg_max = 0;
  int Kg_sum = 0;
  int Kg_max_t = 0;
  int Kg_sum_t = 0;

  if (Pc > 0) {
    S_max = max(S);
    S_sum = sum(S) - Pc;
    Kc_max = max(Kj[1:Pc]);
    Kc_sum = sum(Kj[1:Pc]);
    Kc_max_t = max(Kt[1:Pc]);
    Kc_sum_t = sum(Kt[1:Pc]);
  }
  if (Pg > 0) {
    Kg_max = max(Kj[(Pc + 1):(Pc + Pg)]);
    Kg_sum = sum(Kj[(Pc + 1):(Pc + Pg)]);
    Kg_max_t = max(Kt[(Pc + 1):(Pc + Pg)]);
    Kg_sum_t = sum(Kt[(Pc + 1):(Pc + Pg)]);
  }
}

parameters {

  // distribution specific parameters such as sigma
  vector<lower=0>[Pg] sigma;

  // time-invariant coefficients
  //(including intercept, need to figure out good way choose between tv and tiv intercepts)

  vector[S_sum] alpha_c0;
  vector[S_sum * Kc_sum] beta_c0;

  vector[Pg] alpha_g;
  vector[Pg * Kg_sum] beta_g0;

  // time-varying coefficients
  // no intercept yet (unless in X, but some identification issues likely arise with constant + rw)
  // assume at least for now potentially different standard deviations for each coefficient
  vector<lower=0>[S_sum * Kc_sum_t] sigma_beta_c;
  vector<lower=0>[Pg * Kg_sum_t] sigma_beta_g;

  // standard normal noise for each time step
  matrix[S_sum * Kc_sum_t, T] beta_c0_t_raw;
  matrix[Pg * Kg_sum_t, T] beta_g0_t_raw;

}

transformed parameters {
  matrix[S_sum * Kc_sum_t, T] beta_c0_t;
  matrix[Pg * Kg_sum_t, T] beta_g0_t;
  vector[S_max] alpha_c[Pc];

  if(Pc > 0) alpha_c = alpha_vector(alpha_c0, Pc, S);

  if(Kc_max_t > 0) {
    // t = 1
    // should have prior for the first time point
    // i.e. beta = prior_mean + prior_sd * beta_raw
    beta_c0_t[, 1] = beta_c0_t_raw[, 1];
    for(t in 2:T) {
      beta_c0_t[, t] = beta_c0_t[, t - 1] + sigma_beta_c .* beta_c0_t_raw[,t];
    }
  }

  if(Kg_max_t > 0) {
    beta_g0_t[, 1] = beta_g0_t_raw[, 1];
    for(t in 2:T) {
      beta_g0_t[, t] = beta_g0_t[, t - 1] + sigma_beta_g .* beta_g0_t_raw[,t];
    }
  }
}

model {
  int kj = 0;
  int kt = 0;
  // hard-coded for simplicity here,
  // we can use alpha_c0 ~ normal(alpha_c_mean, alpha_c_sd)
  // where alpha_c_mean and alpha_c_sd are given by the user
  alpha_c0  ~ normal(0, 10);
  beta_c0 ~ normal(0, 10);
  alpha_g ~ normal(0, 10);
  beta_g0 ~ normal(0, 10);
  sigma ~ gamma(2, 0.001);

  to_vector(beta_c0_t_raw) ~ std_normal();
  to_vector(beta_g0_t_raw) ~ std_normal();
  sigma_beta_c ~ gamma(2, 0.001);
  sigma_beta_g ~ gamma(2, 0.001);
  for(i in 1:Pc) {
    matrix[N, Kj[i] + Kt[i]] xi;
    matrix[Kj[i] + Kt[i], S[i]] betai;
    if(Kj[i] > 0) {
      betai[1:Kj[i], 1:S[i]] = append_col(to_matrix(beta_c0[(kj + 1):(kj + Kj[i] * (S[i] - 1))], Kj[i], S[i] - 1), rep_vector(0, Kj[i]));
      kj += Kj[i] * (S[i] - 1);
    }
    for(t in 2:T) {
      if(Kj[i] > 0) {
        xi[, 1:Kj[i]] = x[t, ,Xidx[1:Kj[i], i]];
      }
      if(Kt[i] > 0) {
        xi[, (Kj[i] + 1):(Kj[i] + Kt[i])] = x[t, ,Xidx_t[1:Kt[i], i]];
        betai[(Kj[i] + 1):(Kj[i] + Kt[i]), 1:S[i]] =
        append_col(to_matrix(beta_c0_t[(kt + 1):(kt + Kt[i] * (S[i] - 1)),t], Kt[i], S[i] - 1), rep_vector(0, Kt[i]));

      }
      yc[i, t, ] ~ categorical_logit_glm(xi, alpha_c[i, 1:S[i]], betai);
    }
    kt += Kt[i] * (S[i] - 1);
  }
  kj = 0;
  kt = 0;
  for(i in 1:Pg) {
    for(t in 2:T) {
      vector[N] xbeta = rep_vector(0, N);
      if(Kj[Pc + i] > 0) {
        xbeta = x[t, ,Xidx[1:Kj[Pc + i], i]] * beta_g0[(kj + 1):(kj + Kj[Pc + i])];
      }
      if(Kt[Pc + i] > 0) {
        xbeta += x[t, ,Xidx_t[1:Kt[Pc + i], i]] *  beta_g0_t[(kt + 1):(kt + Kt[Pc + i]), t];
      }
      yg[i, t, ] ~ normal(alpha_g[i] + xbeta, sigma[i]);
    }
    kj += Kj[Pc + i];
    kt += Kt[Pc + i];
  }
}

generated quantities {
  //vector[S_max] alpha_c[Pc];
  matrix[Kc_max, S_max] beta_c[Pc];
  vector[Kg_max] beta_g[Pg];

  matrix[Kc_max_t, S_max] beta_c_t[T, Pc];
  vector[Kg_max_t] beta_g_t[T, Pg];
  if(Pc > 0) {
    if(Kc_max > 0) {
      beta_c = beta_matrix(beta_c0, Pc, S, Kj[1:Pc]);
    }
    if(Kc_max_t > 0) {
      for(t in 1:T) {
        beta_c_t[t] = beta_matrix(beta_c0_t[, t], Pc, S, Kt[1:Pc]);
      }
    }
  }
  if(Pg > 0) {
    if(Kg_max > 0) {
      beta_g = beta_vector(beta_g0, Pg, Kj[(Pc + 1):(Pc + Pg)]);
    }
    if(Kg_max_t > 0) {
      for(t in 1:T) {
        beta_g_t[t] = beta_vector(beta_g0_t[, t], Pg, Kt[(Pc + 1):(Pc + Pg)]);
      }
    }
  }
}
