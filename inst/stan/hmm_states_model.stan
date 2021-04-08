// Hidden markov model for one-channel sequences where each individuals share the common hidden states

functions {
  // Forward algorithm for model where only transition to next state (or stay put) allowed
  /* y observed sequence
  * K number of hidden states
  * T number of time points
  * logA log-transition probabilities from state i to i+1 (strict left-to-right model)
  * logB log-emission probabilities b_ks for symbol s when at state k
  */
  vector forward_algorithm(int[,,] yc, int Pc, int[] S, real [,,] yg,  int Pg, int T, int N, int[] Kj, int[] Kt,
  int n_states, matrix[] x, int[,] Xidx, int[,] Xidx_t,
  vector[] alpha_c, matrix[] beta_c, vector alpha_g, vector[] beta_g, matrix[,] beta_c_t, vector[,] beta_g_t, matrix logA,
  vector sigma) {


    matrix[n_states, T] log_alpha;
    log_alpha[1, 1] = 0;
    log_alpha[2:n_states, 1] = rep_vector(negative_infinity(), n_states - 1);

    // should compute some of these things outside of the loop, needed later as well
    // probably better to use separate function for this piece
    for(i in 1:Pc) {
      int t = 1;
      int k = 1;
      matrix[N, Kj[i] + Kt[i]] xi;
      matrix[Kj[i] + Kt[i], S[i]] betai;
      // if(Kj[i] > 0) {
      //   betai[1:Kj[i], 1:S[i]] = beta_c[i, 1:Kj[i], 1:S[i]];
      //   xi[, 1:Kj[i]] = x[t, ,Xidx[1:Kj[i], i]];
      // }
      //if(Kt[i] > 0) {
        xi[, (Kj[i] + 1):(Kj[i] + Kt[i])] = x[t, ,Xidx_t[1:Kt[i], i]];
        betai[(Kj[i] + 1):(Kj[i] + Kt[i]), 1:S[i]] = beta_c_t[k, i, 1:Kt[i], 1:S[i]];
      //}

      // not yet available
      //yc[i, t, ] ~ categorical_logit_glm(xi, alpha_c[i, 1:S[i]], betai);
      for(n in 1:N) {
        vector[S[i]] theta = alpha_c[i, 1:S[i]] + (xi[n, ] * betai)';
        log_alpha[1, 1] += categorical_logit_lpmf(yc[i, t, n] | theta);
      }

    }
    // for(i in 1:Pg) {
    //   int t = 1;
    //   int k = 1;
    //   vector[N] xbeta = rep_vector(0, N);
    //   if(Kj[Pc + i] > 0) {
    //     xbeta = x[t, ,Xidx[1:Kj[Pc + i], i]] * beta_g[i, 1:Kj[Pc + i]];
    //   }
    //   if(Kt[Pc + i] > 0) {
    //     xbeta += x[t, ,Xidx_t[1:Kt[Pc + i], i]] * beta_g_t[k, i, 1:Kt[Pc + i]];
    //   }
    //
    //   log_alpha[1, 1] += normal_lpdf(yg[i, t, ] | alpha_g[i] + xbeta, sigma[i]);
    // }


    for (t in 2:T) {
      // stay at the starting state 1
      log_alpha[1, t] = log_alpha[1, t - 1] + logA[1, 1];

      for(i in 1:Pc) {
        matrix[N, Kj[i] + Kt[i]] xi;
        matrix[Kj[i] + Kt[i], S[i]] betai;
        // if(Kj[i] > 0) {
        //   betai[1:Kj[i], 1:S[i]] = beta_c[i, 1:Kj[i], 1:S[i]];
        //   xi[, 1:Kj[i]] = x[t, ,Xidx[1:Kj[i], i]];
        // }
        // if(Kt[i] > 0) {
          xi[, (Kj[i] + 1):(Kj[i] + Kt[i])] = x[t, ,Xidx_t[1:Kt[i], i]];
          betai[(Kj[i] + 1):(Kj[i] + Kt[i]), 1:S[i]] = beta_c_t[1, i, 1:Kt[i], 1:S[i]];
       // }

        // not yet available
        //yc[i, t, ] ~ categorical_logit_glm(xi, alpha_c[i, 1:S[i]], betai);
        for(n in 1:N) {
          vector[S[i]] theta = alpha_c[i, 1:S[i]] + (xi[n, ] * betai)';
          log_alpha[1, t] += categorical_logit_lpmf(yc[i,t, n] | theta);
        }

      }
      // for(i in 1:Pg) {
      //   vector[N] xbeta = rep_vector(0, N);
      //   if(Kj[Pc + i] > 0) {
      //     xbeta = x[t, ,Xidx[1:Kj[Pc + i], i]] * beta_g[i, 1:Kj[Pc + i]];
      //   }
      //   if(Kt[Pc + i] > 0) {
      //     xbeta += x[t, ,Xidx_t[1:Kt[Pc + i], i]] * beta_g_t[1, i, 1:Kt[Pc + i]];
      //   }
      //
      //   log_alpha[1, t] += normal_lpdf(yg[i, t, ] | alpha_g[i] + xbeta, sigma[i]);
      // }

      // move (or stay) on state 2 to n_states
      for (k in 2:n_states) {
        real log_py = 0;
        real move = log_alpha[k - 1, t - 1] + logA[k - 1, 2];
        real stay = log_alpha[k, t - 1]   + logA[k, 1];
        for(i in 1:Pc) {
          matrix[N, Kj[i] + Kt[i]] xi;
          matrix[Kj[i] + Kt[i], S[i]] betai;
          // if(Kj[i] > 0) {
          //   betai[1:Kj[i], 1:S[i]] = beta_c[i, 1:Kj[i], 1:S[i]];
          //   xi[, 1:Kj[i]] = x[t, ,Xidx[1:Kj[i], i]];
          // }
          // if(Kt[i] > 0) {
            xi[, (Kj[i] + 1):(Kj[i] + Kt[i])] = x[t, ,Xidx_t[1:Kt[i], i]];
            betai[(Kj[i] + 1):(Kj[i] + Kt[i]), 1:S[i]] = beta_c_t[k, i, 1:Kt[i], 1:S[i]];
          //}

          // not yet available
          //yc[i, t, ] ~ categorical_logit_glm(xi, alpha_c[i, 1:S[i]], betai);
          for(n in 1:N) {
            vector[S[i]] theta = alpha_c[i, 1:S[i]] + (xi[n, ] * betai)';
            log_py += categorical_logit_lpmf(yc[i,t, n] | theta);
          }

        }
        // for(i in 1:Pg) {
        //   vector[N] xbeta = rep_vector(0, N);
        //   if(Kj[Pc + i] > 0) {
        //     xbeta = x[t, ,Xidx[1:Kj[Pc + i], i]] * beta_g[i, 1:Kj[Pc + i]];
        //   }
        //   if(Kt[Pc + i] > 0) {
        //     xbeta += x[t, ,Xidx_t[1:Kt[Pc + i], i]] * beta_g_t[k, i, 1:Kt[Pc + i]];
        //   }
        //
        //   log_py += normal_lpdf(yg[i, t, ] | alpha_g[i] + xbeta, sigma[i]);
        // }

        move += log_py;
        stay += log_py;
        // Stan gives error if both terms are -Inf
        // need to check if it faster to just add tiny values to all initial state probabilities
        if(!is_inf(move) || !is_inf(stay)) {
          log_alpha[k, t] = log_sum_exp(move, stay);
        } else {
          log_alpha[k, t] = negative_infinity();
        }
      }
    }
    // we need all log_alphas later if we want to extract posterior probabilities of the states etc
    return log_alpha[, T];
  }

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

  int<lower=1> n_states; // number of latent states
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
  // transition probabilities of hidden states
  vector<lower=0,upper=1>[n_states - 1] A;

  // distribution specific parameters such as sigma (could vary in time)
  vector<lower=0>[Pg] sigma;

  // time-invariant coefficients
  //(including intercept, need to figure out good way choose between tv and tiv intercepts)

  vector[S_sum] alpha_c0;
  vector[S_sum * Kc_sum] beta_c0;

  vector[Pg] alpha_g;
  vector[Pg * Kg_sum] beta_g0;

  // time-varying coefficients, one column for each state
  matrix[S_sum * Kc_sum_t, n_states] beta_c0_t;
  matrix[Pg * Kg_sum_t, n_states] beta_g0_t;
}

transformed parameters {

  vector[S_max] alpha_c[Pc];

  matrix[Kc_max, S_max] beta_c[Pc];
  vector[Kg_max] beta_g[Pg];
  matrix[Kc_max_t, S_max] beta_c_t[n_states, Pc];
  vector[Kg_max_t] beta_g_t[n_states, Pg];

  matrix[n_states, 2] logA;
  logA[1:(n_states - 1),] = log(append_col(A, 1 - A));
  logA[n_states, 1] = 0;
  logA[n_states, 2] = negative_infinity();



 // if(Pc > 0) {
    alpha_c = alpha_vector(alpha_c0, Pc, S);
    // if(Kc_max > 0) {
    //   beta_c = beta_matrix(beta_c0, Pc, S, Kj[1:Pc]);
    // }
    // if(Kc_max_t > 0) {
      for(t in 1:n_states) {
        beta_c_t[t] = beta_matrix(beta_c0_t[, t], Pc, S, Kt[1:Pc]);
      }
   // }
  //}
  // if(Pg > 0) {
  //   if(Kg_max > 0) {
  //     beta_g = beta_vector(beta_g0, Pg, Kj[(Pc + 1):(Pc + Pg)]);
  //   }
  //   if(Kg_max_t > 0) {
  //     for(t in 1:n_states) {
  //       beta_g_t[t] = beta_vector(beta_g0_t[, t], Pg, Kt[(Pc + 1):(Pc + Pg)]);
  //     }
  //   }
  // }

}
model {
  // more likely to stay than to move
  // prior probability for staying 0.9, probably should depend on T
  A ~ beta(22.5, 2.5);

  // hard-coded for simplicity here,
  // we can use alpha_c0 ~ normal(alpha_c_mean, alpha_c_sd)
  // where alpha_c_mean and alpha_c_sd are given by the user
  alpha_c0  ~ normal(0, 10);
  beta_c0  ~ normal(0, 10);
  // alpha_g ~ normal(0, 10);
  // beta_g0 ~ normal(0, 10);
  // sigma ~ normal(0, 10);

  to_vector(beta_c0_t) ~ normal(0, 10);
  //to_vector(beta_g0_t) ~ normal(0, 10);
  target += log_sum_exp(forward_algorithm(yc, Pc, S, yg, Pg, T, N, Kj, Kt,
  n_states,x, Xidx, Xidx_t, alpha_c, beta_c, alpha_g, beta_g, beta_c_t, beta_g_t, logA, sigma));
}
