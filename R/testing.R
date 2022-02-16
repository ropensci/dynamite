# Dummy data for stanmodel conversion

T <- 10
N <- 4
TN <- T * N

test_data <- data.frame(y1 = as.factor(sample(2, size = TN, replace = TRUE)),
                        y2 = as.factor(sample(3, size = TN, replace = TRUE)),
                        y3 = as.factor(sample(5, size = TN, replace = TRUE)),
                        x1 = rnorm(TN),
                        x2 = as.factor(sample(4, size = TN, replace = TRUE)),
                        x3 = rnorm(TN),
                        x4 = rnorm(TN),
                        ID = gl(N, T))

test_form <- obs(y1 ~ x1 + x2 + x4, family = categorical()) +
    obs(y2 ~ x1 + x3 + x4, family = categorical()) +
    obs(y3 ~ x1 + x3, family = categorical()) +
    splines() +
    lags()

# Same with manual lags
test_form2 <- obs(y1 ~ x1 + x2 + x4 + lag(y1, 1) + lag(y2, 1) + lag(y3, 1), family = categorical()) +
    obs(y2 ~ x1 + x3 + x4 + lag(y1, 1) + lag(y2, 1) + lag(y3, 1), family = categorical()) +
    obs(y3 ~ x1 + x3 + lag(y1, 1) + lag(y2, 1) + lag(y3, 1), family = categorical()) +
    splines()

# Use these to control what is returned, wrap in argument debug:
#     no_compile = TRUE => stan model is not compiled, implies no_sampling = TRUE
#     no_sampling = TRUE => sampling is not called, stanfit == NULL
#     model_matrix = TRUE => also return the model matrix
#     model_data = TRUE => also return the data object passed to stan
#     model_code = TRUE => also return the model code

# Rstudio ignores .Rbuildignore, so this needs to be commented out when building the package
# test_fit <- btvcm:::btvcmfit(test_form, test_data, ID, debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))
# test_fit2 <- btvcm:::btvcmfit(test_form2, test_data, ID, debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))

# Should give identical model matrices
# all.equal(test_fit$model_matrix, test_fit2$model_matrix)

set.seed(1)
T <- 20
N <- 10
x <- matrix(rnorm(T*N), N, T)
intercept <- t(replicate(N, cumsum(rnorm(T, sd = 0.2))))
y <- matrix(0, N, T+1)
y[, 1] <- rnorm(N)
for(t in 1:T) {
    y[, t+1] <- intercept[, t] + x[, t] + 0.5 * y[, t] + rnorm(N, sd = 0.1)
}
y <- y[, -1]
ts.plot(t(y))
d <- data.frame(y = c(t(y)), x = c(t(x)), ID = gl(N, T))

#fit <- btvcm:::btvcmfit(obs(y ~ x, family = gaussian()) + lags() + splines(df = 5), d, ID, chains = 1,
#    debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))

# print(fit$stanfit, "tau_1")
