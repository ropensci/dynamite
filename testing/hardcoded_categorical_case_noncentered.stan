
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
    row_vector[D] a_raw[S-1, K]; // S - 1 as we need to fix one set of betas to zero
    vector<lower=0>[K] tau;
}

transformed parameters {
    matrix[S, K] beta[T];
    row_vector[D] a[S-1, K];

    for(s in 1:(S-1)) {
        for(k in 1:K) {
            a[s, k, 1] = a_raw[s, k, 1];
            for (i in 2:D) {
                a[s, k, i] = a[s, k, i-1] + tau[k] * a_raw[s, k, i];
            }
            for(t in 1:T){
                beta[t, s, k] = a[s, k] * B[, t];
            }
        }
    }
    for(t in 1:T){
        beta[t, S, ] = zeros_K;
    }
}

model {
    tau ~ std_normal();
    for(s in 1:(S-1)) {
        for(k in 1:K) {
            a_raw[s, k] ~ std_normal();
        }
    }
    for(t in 1:T) {
        y[t] ~ categorical_logit_glm(X[t], zeros_S, beta[t]');
    }
}
