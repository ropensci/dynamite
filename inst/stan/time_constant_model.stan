functions {

  vector[] alpha_vector(vector alphav, int Pc, int[] S, int S_max, int[] S_cumsum) {
    // intercept terms for the sequences
    // alphav:
    // (S_cumsum[1]+1):S_cumsum[2] = 1:(S[1]-1)for the first series
    // (S_cumsum[2]+1):S_cumsum[3] = S[1]:(S[1])  for the second etc
    // (S[1] + S[2]):(S[1] + S[2] + S[3] - 1) for the third
    vector[S_max] alpha[Pc];
    for(i in 1:Pc) {
      alpha[i] = rep_vector(0, S_max);
      alpha[i, 1:(S[i]-1)] = alphav[(S_cumsum[i] + 1):S_cumsum[i+1]];
    }
    return alpha;
  }

  matrix[] beta_matrix(vector betav, int Pc, int[] S, int S_max, int[] K, int K_max, int[] KS_cumsum) {
    // regression coefficients terms for sequences
    // betav:
    // 1:(K[1] * (S[1] - 1)) for the first series
    // (K[1] * (S[1] - 1) + 1):(K[1] * (S[1] - 1) + 1 + K[2] * (S[2] - 1)) for the second
    matrix[K_max, S_max] beta[Pc];
    for(i in 1:Pc) {
      beta[i] = rep_matrix(0, K_max, S_max);
      beta[i, 1:K[i], 1:(S[i] - 1)] = to_matrix(betav[(KS_cumsum[i] + 1):KS_cumsum[i+1]], K[i], S[i] - 1);
    }
    return beta;
  }

  vector[] beta_vector(vector betav, int P, int[] K, int K_max, int[] K_cumsum) {
    vector[K_max] beta[P];
    for(i in 1:P) {
      beta[i] = rep_vector(0, K_max);
      beta[i, 1:K[i]] = betav[(K_cumsum[i] + 1):K_cumsum[i + 1]];
    }
    return beta;
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
  int S_cumsum[Pc + 1] = rep_array(0, Pc + 1);
  int S_sum = 0;
  int Kc_max = 0;
  int Kc_sum = 0;
  int KS_cumsum[Pc + 1] = rep_array(0, Pc + 1);

  int Kg_max = 0;
  int Kg_sum = 0;
  int Kg_cumsum[Pg + 1] = rep_array(0, Pg + 1);

  if (Pc > 0) {
    S_max = max(S);
    S_sum = sum(S) - Pc;
    Kc_max = max(Kj[1:Pc]);
    Kc_sum = sum(Kj[1:Pc]);
    S_cumsum[1] = 0;
    KS_cumsum[1] = 0;
    for(i in 2:(Pc + 1)) {
      S_cumsum[i] = S_cumsum[i-1] + S[i-1] - 1;
      KS_cumsum[i] = KS_cumsum[i-1] + Kj[i-1] * (S[i-1] - 1);
    }
  }
  if (Pg > 0) {
    Kg_max = max(Kj[(Pc + 1):(Pc + Pg)]);
    Kg_sum = sum(Kj[(Pc + 1):(Pc + Pg)]);
    Kg_cumsum[1] = 0;
    for(i in 2:(Pg + 1)) {
      Kg_cumsum[i] = Kg_cumsum[i-1] + Kj[Pc + i - 1];
    }
  }
}

parameters {
  // intercept terms for sequences
  // 1:(S[1] - 1) for first type
  // S[1]:(S[1] + S[2] - 1 for second etc)
  vector[S_sum] alpha_c0;
  vector[S_sum * Kc_sum] beta_c0;
  vector[S_sum] pi0;
  vector<lower=0>[Pg] sigma;
  vector[Pg] alpha_g;
  vector[Pg * Kg_sum] beta_g0;
}

transformed parameters {
  vector[S_max] alpha_c[Pc] = alpha_vector(alpha_c0, Pc, S, S_max, S_cumsum);
  matrix[Kc_max, S_max] beta_c[Pc] = beta_matrix(beta_c0, Pc, S, S_max, Kj[1:Pc], Kc_max, KS_cumsum);
  vector[S_max] pi[Pc] = alpha_vector(pi0, Pc, S, S_max, S_cumsum);
  vector[Kg_max] beta_g[Pg] = beta_vector(beta_g0, Pg, Kj[(Pc + 1):(Pc + Pg)], Kg_max, Kg_cumsum);
}
model {
  alpha_c0  ~ normal(0, 5);
  beta_c0  ~ normal(0, 5);
  pi0 ~ normal(0, 5);
  alpha_g ~ normal(0, 5);
  beta_g0 ~ normal(0, 5);
  sigma ~ normal(0, 5);
  for(i in 1:Pc) {
    for(ii in 1:N) yc[i, 1, ii] ~ categorical_logit(pi[i, 1:S[i]]);
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
    // QR parameterization could be better (see Stan manual and Case study), also for categorical case...
    // assume first time point fixed for simplicity at least for now
    for(t in 2:T) {
      yg[i, t, ] ~ normal(alpha_g[i] + x[t, , Xidx[1:Kj[i], i]] * beta_g[i, 1:Kj[i]], sigma[i]);
    }
  }
}
