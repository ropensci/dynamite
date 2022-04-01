# Dummy data for stanmodel conversion

T <- 10
N <- 4
TN <- T * N

test_data <- data.frame(
    y1 = as.factor(sample(2, size = TN, replace = TRUE)),
    y2 = as.factor(sample(3, size = TN, replace = TRUE)),
    y3 = as.factor(sample(5, size = TN, replace = TRUE)),
    y4 = sample(5, size = TN, replace = TRUE),
    y5 = rnorm(n = TN, mean = 1, sd = 2),
    x1 = rnorm(TN),
    x2 = as.factor(sample(4, size = TN, replace = TRUE)),
    x3 = rnorm(TN),
    x4 = rnorm(TN),
    x5 = sample(letters[1:4], size = TN, replace = TRUE),
    t = sample(50:100, size = TN, replace = TRUE),
    ID = gl(N, T),
    time = 1:T)

test_data$n <- test_data$y4 + sample(10, size = TN, replace = TRUE)

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

test_all <- obs(y4 ~ -1 + x1 + varying(~x2) + offset(log(t)), family = poisson()) +
    obs(y2 ~ -1 + x1 + varying(~x3) + trials(n), family = binomial()) +
    obs(y5 ~ x4 + fixed(~x5) + varying(~-1 + lag(y5)), family = gaussian()) +
    obs(y1 ~ x1 + x2 + x4 + lag(y1, 1) + lag(y2, 1), family = categorical()) +
    splines()
#test_fit <- btvcm:::btvcmfit(test_all, test_data, ID, t, debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))
#test_fit <- btvcm:::btvcmfit(test_all, test_data, ID, t)

test_trials <- obs(y4 ~ x1 + trials(n), family = binomial()) + splines()
#test_fit <- btvcm:::btvcmfit(test_trials, test_data, ID, t, debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))
#test_fit <- btvcm:::btvcmfit(test_trials, test_data, ID, t)

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

#
# set.seed(1)
# T <- 20
# N <- 500
# x <- matrix(rnorm(T*N), N, T)
# intercept <- cumsum(rnorm(T))
# y <- matrix(0, N, T+1)
# y[, 1] <- rnorm(N)
# betax <- cumsum(rnorm(T, sd = 0.1))
# betay <- cumsum(rnorm(T, sd = 0.1))
# for(t in 1:T) {
#     y[, t+1] <- intercept[t] + betax[t] * x[, t] + betay[t] * y[, t] + rnorm(N, sd = 0.1)
# }
# y <- y[, -1]
# ts.plot(t(y))
# d <- data.frame(y = c(t(y)), x = c(t(x)), ID = gl(N, T), month = 1:T)
# #
# fit <- btvcm:::btvcmfit(
#     obs(y ~ x, family = gaussian()) +
#         lags() +
#         splines(knots = 3:19), #t=1 fixed
#     d, ID, month, chains = 1)
#
# print(fit$stanfit, pars = c("sigma_1", "tau_1"))
# b <- apply(rstan::extract(fit$stanfit, "beta_1")[[1]], 2:3, mean)
# ts.plot(cbind(b, intercept[-1], betax[-1], betay[-1]), col=1:3,lty=rep(1:2,each=3))
#
# cf <- coef(fit) %>% group_by(time, variable) %>%
#     summarise(mean = mean(value),
#         lwr = quantile(value, 0.025), upr = quantile(value, 0.975))
# cf$true <- c(rbind(intercept, betay, betax)[, -1])
# cf %>%
#     ggplot(aes(time, mean)) + theme_bw() +
#     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7) +
#     geom_line() + geom_line(aes(y = true), colour = "red") +
#     facet_wrap(~ variable, scales = "free_y")
#
# newdata <- d %>% filter(ID == 1)
# newdata$x[] <- 0
# newdata$y[] <- 2
# fits <- fitted(fit, newdata = newdata)
# sumr <- fits %>% group_by(ID, month) %>% summarise(mean = mean(y), lwr = quantile(y, 0.025), upr = quantile(y, 0.975))
# sumr$y <- d$y[d$ID == 1]
# sumr %>%
#     ggplot(aes(month, mean, group = ID)) +
#     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7, fill = "blue") +
#     geom_line(aes(y = y), colour = "red") +
#     geom_line() +
#     theme_bw()
#
# newdata <- d %>% filter(ID == 1)
# newdata$x[] <- 0
# # the uncertainty is minimal due to the strong correlation of intercept and beta_y, and fixed negative y
# fits <- fitted(fit, newdata = newdata)
# sumr <- fits %>% group_by(ID, month) %>% summarise(mean = mean(y), lwr = quantile(y, 0.025), upr = quantile(y, 0.975))
# sumr$y <- d$y[d$ID == 1]
# sumr %>%
#     ggplot(aes(month, mean, group = ID)) +
#     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7, fill = "blue") +
#     geom_line(aes(y = y), colour = "red") +
#     geom_line() +
#     theme_bw()
#
# # simulate also y:
# newdata <- d %>% filter(ID == 1)
# newdata$x[] <- 0
# newdata$y[-1] <- NA
# pred <- predict(fit, newdata = newdata)
# sumr <- pred %>% group_by(ID, month) %>% summarise(mean = mean(y), lwr = quantile(y, 0.025), upr = quantile(y, 0.975))
# sumr$y <- d$y[d$ID == 1]
# sumr %>%
#     ggplot(aes(month, mean, group = ID)) +
#     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7, fill = "blue") +
#     geom_line(aes(y = y), colour = "red") +
#     geom_line() +
#     theme_bw()



# newdata <- d
# newdata$y[newdata$month > 1] <- NA
# system.time(pred <- predict(fit, newdata = newdata))
# sumr <- pred %>% group_by(ID, month) %>% summarise(mean = mean(y), lwr = quantile(y, 0.05), upr = quantile(y, 0.95))
#
# sumr %>% ggplot(aes(month, mean, group = ID)) + geom_line()
# print(fit$stanfit, "tau_1")
# b <- apply(rstan::extract(fit$stanfit, "beta_1")[[1]], 2:3, mean)
# ts.plot(cbind(b, intercept[-1], betax[-1], betay[-1]), col=1:3,lty=rep(1:2,each=3))
#
# set.seed(1)
# T <- 20
# N <- 500
# x <- matrix(rnorm(T*N), N, T)
# intercept <- 1
# y <- matrix(0, N, T+1)
# y[, 1] <- rnorm(N)
# for(t in 1:T) {
#     y[, t+1] <- intercept +  x[, t] + 0.5 * y[, t] + rnorm(N, sd = 0.1)
# }
# y <- y[, -1]
# ts.plot(t(y))
# d <- data.frame(y = c(t(y)), x = c(t(x)), ID = gl(N, T))

# fit <- btvcm:::btvcmfit(
#     obs(y ~ x, family = gaussian()) +
#         lags() +
#         splines(df = 10),
#     d, ID, chains = 1, cores = 1)
#      #, debug = list(no_compile = TRUE, model_matrix = TRUE, model_data = TRUE, model_code = TRUE))
# print(fit$stanfit, "beta_1")
# print(fit$stanfit, "tau_1")
#
