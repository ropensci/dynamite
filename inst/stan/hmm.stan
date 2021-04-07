// Hidden markov model for one-channel sequences where each individual has own hidden paths

functions {
  // Forward algorithm for model where only transition to next state (or stay put) allowed
  /* y observed sequence
  * K number of hidden states
  * T number of time points
  * logA log-transition probabilities from state i to i+1 (strict left-to-right model)
  * logB log-emission probabilities b_ks for symbol s when at state k
  */
  vector forward_algorithm(int[] y, int K, int T, matrix logA, matrix logB) {

    matrix[K, T] log_alpha;

    log_alpha[1, 1] = logB[1, y[1]];
    log_alpha[2:K, 1] = rep_vector(negative_infinity(), K - 1);

    for (t in 2:T) {
      // stay at the starting state 1
      log_alpha[1, t] = log_alpha[1, t - 1] + logA[1, 1] + logB[1, y[t]];
      // move (or stay) on state 2 to K
      for (k in 2:K) {
        real move = log_alpha[k - 1, t - 1] + logA[k - 1, 2] + logB[k, y[t]];
        real stay = log_alpha[k, t - 1]   + logA[k, 1]   + logB[k, y[t]];
        // Stan gives error if both terms are -Inf
        // need to check if it faster to just add tiny values to all initial state probabilities
        if(!is_inf(move) || !is_inf(stay)) {
          log_alpha[k, t] = log_sum_exp(move, stay);
        } else {
          log_alpha[k, t] = negative_infinity();
        }
      }
    }
    return log_alpha[, T];
  }
}

data {
  int<lower=1> N; // number of sequences
  int<lower=1> T; // number of subsequences
  int<lower=1> S; // number of symbols
  int<lower=1> K; // number of latent states
  int<lower=1> y[T, N]; //observations
}
transformed data {
  vector<lower=0>[S] ones_S = rep_vector(1, S);
  vector<lower=0>[K] ones_K = rep_vector(1, K);

}
parameters {
  simplex[S] B[K];
  vector<lower=0,upper=1>[K-1] A;
}

transformed parameters {
  vector[N] loglik;
  matrix[K, 2] logA;
  matrix[K, S] logB;
  logA[1:(K-1),] = log(append_col(A, 1 - A));
  logA[K, 1] = 0;
  logA[K, 2] = negative_infinity();
  for (k in 1:K) {
    logB[k, ] = log(B[k])';
  }
  for(n in 1:N) {
    loglik[n] = log_sum_exp(forward_algorithm(y[, n], K, T, logA, logB));
  }
}
model {
  // more likely to stay than to move
  // prior probability for staying 0.9, probably should depend on T
  A ~ beta(22.5, 2.5);
  for(k in 1:K) {
    B[k] ~ dirichlet(ones_S);
  }
  target += sum(loglik);
}
