.loglik_stan <- "
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
                for (i in 1:M) {
                    if (tmp[i] > negative_infinity()) {
                        n_finite += 1;
                        tmp[n_finite] = tmp[i];
                    }
                }
                if (n_finite == 0) {
                    log_alpha_new[k] = negative_infinity();
                } else {
                    log_alpha_new[k] = log_sum_exp(tmp[1:n_finite]);
                }
            }
            log_alpha = log_alpha_new;
        }
        return log_sum_exp(log_alpha);
    }"
