set.seed(1)
n_id <- 500
n_time <- 50
d <- data.frame(y = sample(factor(c("a", "b", "c")), size = n_id, replace = TRUE),
                x = sample(factor(c("A", "B", "C")), size = n_id, replace = TRUE),
                time = 1, id = 1:n_id)
d <- dplyr::right_join(d,
                       data.frame(
                         time = rep(1:n_time, each = n_id),
                         id = 1:n_id))

d$z <- rnorm(nrow(d))
f <-
  obs(x ~ z + lag(x) + lag(y), family = categorical()) +
  obs(y ~ z + lag(x) + lag(y), family = categorical())

get_priors(f, d, "id", "time")
init <- list(beta_x = matrix(c(c(2, 0.8, 0.2, 0, 0,
                                 1, 0.5, 2, 0.2, 0.1)), 5, 2),
             a_x = c(-0.1, 0.2),
             beta_y = matrix(c(c(0, 1, 0.8, 0.3, 0.5,
                                 1, 0.2, 0, 0.3, -0.5)), 5, 2),
             a_y = c(0.1, -0.5))

fit <- dynamite(f,  d, "id", "time",
                chains = 1, iter = 1,
                algorithm = "Fixed_param", init = list(init))


p <- predict(fit, type = "response")
#
# p |> dplyr::select(id, time, x_new) |>
#   tidyr::pivot_wider(names_from = time, values_from = x_new) |>
#   TraMineR::seqdef(1:n_time + 1) |>
#   TraMineR::seqdplot()
#
#
# p |> dplyr::select(id, time, y_new) |>
#   tidyr::pivot_wider(names_from = time, values_from = y_new) |>
#   TraMineR::seqdef(1:n_time + 1) |>
#   TraMineR::seqdplot()

categorical_example <- p |> dplyr::mutate(x = x_new, y = y_new) |>
  dplyr::select(id, time, x, y, z)

usethis::use_data(categorical_example, overwrite = TRUE)
