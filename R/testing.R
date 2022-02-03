# Dummy data for stanmodel conversion

T <- 25
N <- 4
TN <- T * N

test_data <- data.frame(y1 = rnorm(TN),
                        y2 = rnorm(TN),
                        y3 = sample(5, size = 100, replace = TRUE),
                        x1 = rnorm(TN),
                        x2 = rnorm(TN),
                        x3 = rnorm(TN),
                        x4 = rnorm(TN),
                        ID = gl(N, T))

test_form <- obs(y1 ~ x1 + x2 + x4, "gaussian") +
    obs(y2 ~ x1 + x3 + x4, "gaussian") +
    obs(y3 ~ x1 + x3, "categorical") +
    states(n = 2)

test_fit <- because(test_form, test_data)
