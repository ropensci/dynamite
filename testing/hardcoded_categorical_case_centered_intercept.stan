
data {

    int<lower=1> T; // number of time points
    int<lower=1> N; // number of individuals
    int y[T,N];
    int<lower=0> K;
    matrix[N, K] X[T]; // past observations as covariates, see generated quantities!
    int<lower=0> S; // number of observed symbols
    int D;
    matrix[D, T] B;

}
transformed data {
    row_vector[K] zeros_K = rep_row_vector(0, K);
    vector[S] zeros_S = rep_vector(0, S);
}

parameters {
    row_vector[D] a[S-1, K]; // S - 1 as we need to fix one set of betas to zero
    vector<lower=0>[K] tau;
    matrix[S-1, K] beta0;
}

transformed parameters {
    matrix[S, K] beta[T];

    for(s in 1:(S-1)) {
        for(k in 1:K) {
            for(t in 1:T){
                beta[t, s, k] = beta0[s, k] + a[s, k] * B[, t];
            }
        }
    }
    for(t in 1:T){
        beta[t, S, ] = zeros_K;
    }
}

model {
    to_vector(beta0) ~ normal(0, 2);
    tau ~ std_normal();
    for(s in 1:(S-1)) {
        for(k in 1:K) {
            a[s, k, 1] ~ std_normal();
            for (i in 2:D) {
                a[s, k, i] ~ normal(a[s, k, i-1], tau[k]);
            }
        }
    }

    for(t in 1:T) {
        y[t] ~ categorical_logit_glm(X[t], zeros_S, beta[t]');
    }
}
