// Double chain model for Gaussian sequences 
// where each individuals have common hidden state process

functions {
  /* Log-likelihood function where we check for -Inf values
  *  We need this because of autodiff issues
  */
  real loglik(matrix log_py, int M, int T, matrix log_A) {
    
    vector[M] log_alpha;
    vector[M] log_alpha_new;
    vector[M] tmp;
    
    log_alpha[1] = log_py[1, 1];
    log_alpha[2:M] = rep_vector(negative_infinity(), M - 1);
    for (t in 2:T) {
      for (k in 1:M) {
        // variable tmp below is constructed as theta + x where x might contain -Inf elements
        int n_finite = 0;
        tmp = log_alpha + log_A[, k] + log_py[k, t];
        // reorder elements of tmp so that all finite values are at the start of the vector
        for(i in 1:M) {
          if(tmp[i] > negative_infinity()) {
            n_finite += 1;
            tmp[n_finite] = tmp[i];
          }
        }
        if(n_finite == 0) {
          log_alpha_new[k] = negative_infinity();
        } else {
          log_alpha_new[k] = log_sum_exp(tmp[1:n_finite]);
        }
      }
      
      log_alpha = log_alpha_new;
    }
    return log_sum_exp(log_alpha);
  }
  
  /* Viterbi algorithm
  * We always assume initial state distribution pi as (1, 0, ..., 0)
  * transition matrix A is unconstrained
  * 
  * log_py: Observational log-densities as M x T matrix
  * M:    Number of hidden states
  * T:    Number of time points
  * log_A: logs of transition probabilities row->column, M x M matrix
  */
  int[] viterbi(matrix log_py, int M, int T, matrix log_A) {
    
    int viterbi_path[T];
    real logp;
    int psi[T, M];
    matrix[M, T] delta;
    // fixed starting state
    delta[1, 1] = log_py[1, 1];
    delta[2:M, 1] = rep_vector(negative_infinity(), M - 1);
    
    for (t in 2:T) {
      for (j in 1:M) {
        // j = current (t)
        delta[j, t] = negative_infinity();
        for (i in 1:M) {
          // i = previous (t-1)
          logp = delta[i, t-1] + log_A[i, j] + log_py[j, t];
          if (logp > delta[j, t]) {
            psi[t, j] = i;
            delta[j, t] = logp;
          }
        }
      }
    }
    
    logp = max(delta[,T]);
    
    for (j in 1:M) {
      if (delta[j, T] == logp) {
        viterbi_path[T] = j;
      }
    }
    
    for (t in 1:(T - 1)) {
      viterbi_path[T - t] = psi[T - t + 1, viterbi_path[T - t + 1]];
    }
    
    return viterbi_path;
  }
  
}

data {
  int<lower=1> M; // number of latent states
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int<lower=1> C; // Number of channels
  matrix[N, C] y[T]; //observations
  int<lower=0> T_fixed_start; // number of time points fixed to first state
  int<lower=0> T_fixed_end; // number of time points fixed to last state
}
transformed data {
  matrix[M, M] inf_mat = rep_matrix(negative_infinity(), M, M);
}
parameters {
  
  vector<lower=0,upper=1>[M-1] A;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_intercept;
  real<lower=0> sigma_sigma;
  // use noncentered paramterisation for sampling efficiency
  matrix[C, C] beta_raw[M];
  vector[C] intercept_raw[M]; 
  vector[C] log_sigma_raw[M];
}

transformed parameters {  
  // log(p(y_t | y_t-1, z_t))
  matrix[M, T] log_py = rep_matrix(0, M, T);
  // regression coefficients
  // col c gives coefficients for response channel c, i.e.
  // y^c_t = beta_{c1} * y^1_{t-1} + ... + beta_{cC} * y^C_{t-1}
  matrix[C, C] beta[M];
  vector[C] intercept[M]; // intercepts
  vector<lower=0>[C] sigma[M];
  matrix[M, M] log_A = inf_mat;
  for(i in 1:(M-1)) {
    log_A[i, i] = log(A[i]);
    log_A[i, i+1] = log(1 - A[i]);
  }
  log_A[M, M] = 0;
  {
    vector[C] log_sigma_tmp = log_sigma_raw[1];
    sigma[1] = exp(log_sigma_tmp);
    beta[1] = beta_raw[1];
    intercept[1] = intercept_raw[1];
    for(m in 2:M) {
      beta[m] = beta[m-1] + sigma_beta * beta_raw[m];
      intercept[m] = intercept[m-1] + sigma_intercept * intercept_raw[m];
      log_sigma_tmp += sigma_sigma * log_sigma_raw[m];
      sigma[m] = exp(log_sigma_tmp);
    }
  }
  
  
  // assume first time point as fixed for simplicity
  
  // assume that at least the first T_fixed time points are from state 1
  // This should help with multimodality together with random walk prior
  log_py[2:M, 1:T_fixed_start] += negative_infinity();
  for(t in 2:T_fixed_start) {
    for(c in 1:C) {
      log_py[1, t] += normal_lpdf(y[t, , c] | intercept[1, c] + y[t-1] * beta[1, , c], sigma[1, c]);
    }
  }
  // similarly assume that we actually end up in the last state at least
  // for last T_fixed_end time points
  for(t in (T_fixed_start + 1):(T - T_fixed_end)) {
    for(m in 1:M) {
      for(c in 1:C) {
        log_py[m, t] += normal_lpdf(y[t, , c] | intercept[m, c] + y[t-1] * beta[m, , c], sigma[m, c]);
      }
    }
  } 
  log_py[1:(M-1), (T - T_fixed_end + 1):T] += negative_infinity();
  for(t in (T - T_fixed_end + 1):T) {
    for(c in 1:C) {
      log_py[M, t] += normal_lpdf(y[t, , c] | intercept[M, c] + y[t-1] * beta[M, , c], sigma[M, c]);
    }
  }
}

model {
  
  A ~ beta(T, M); // prior average holding time (T + M) / M
  sigma_beta ~ std_normal();
  sigma_intercept ~ std_normal();
  sigma_sigma ~ std_normal();
  for(m in 1:M){
    to_vector(beta_raw[m]) ~ std_normal();
    intercept_raw[m] ~ std_normal();
    log_sigma_raw[m] ~ std_normal();
  }
  
  target += loglik(log_py, M, T, log_A);
}

// generated quantities {
//    int<lower=1,upper=M> viterbi_path[T] = viterbi(log_py, M, T, log_A);
// }
