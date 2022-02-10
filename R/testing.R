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

# Rstudio ignores .Rbuildignore, so this needs to be commented out when building the package
# test_fit <- btvcm:::btvcmfit(test_form, test_data, ID)

