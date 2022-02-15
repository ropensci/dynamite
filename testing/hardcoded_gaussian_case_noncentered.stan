
data {

    int<lower=1> T; // number of time points
    int<lower=1> N; // number of individuals
    real y[T,N];
    int<lower=0> K;
    matrix[N, K] X[T];
    int D;
    matrix[D, T] B;

}

parameters {
    row_vector[D] a_raw[K];
    vector<lower=0>[K] tau;
    real<lower=0> sigma;
}

transformed parameters {
    vector[K] beta[T];
    row_vector[D] a[K];

    for(k in 1:K) {
        a[k, 1] = a_raw[k, 1];
        for (i in 2:D) {
            a[k, i] = a[k, i-1] + tau[k] * a_raw[k, i];
        }
        for(t in 1:T){
            beta[t, k] = a[k] * B[, t];
        }
    }
}

model {
    tau ~ std_normal();
    for(k in 1:K) {
        a_raw[k] ~ std_normal();
    }
    sigma ~ std_normal();
    for(t in 1:T) {
        y[t] ~ normal(X[t] * beta[t], sigma);
    }
}


