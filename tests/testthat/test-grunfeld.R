#' @srrstats {G5.0, G5.4a, G5.4b, G5.4c} The parameter estimates of
#'   model fitted to a classic Grunfeld data match with brms and plm.
#'
run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "1")

test_that("parameters of the Grunfield model are recovered", {

  skip_if_not(run_extended_tests)
  library("plm")


  data(Grunfeld, package = "plm")
  fit_plm <- plm(inv ~ value + capital, data = Grunfeld,
    index = c("firm", "year"), effect = "individual", model = "within")

  set.seed(1)
  # set very vague priors
  p <- get_priors(obs(inv ~ value + capital,
    family = "gaussian", random_intercept = TRUE),
    Grunfeld, "firm", "year")
  p$prior[] <- rep("normal(0, 1000)", nrow(p))
  fit <- dynamite(obs(inv ~ value + capital,
    family = "gaussian", random_intercept = TRUE),
    Grunfeld, "firm", "year", refresh = 0,
    chains = 1, cores = 1, iter = 2000)

  sumr <- summary(fit)

  expect_equal(coef(fit_plm), coef(fit)$mean[3:2],
    tolerance = 0.01, ignore_attr = TRUE)
  library("brms")
  fit_brm <- brm(inv ~ value + capital + (1 | firm), Grunfeld,
    prior = prior(normal(0, 1000), class = "b") +
      prior(normal(0, 1000), class = "Intercept") +
      prior(normal(0, 1000), class = "sd") +
      prior(normal(0, 1000), class = "sigma"),
    refresh = 0,
    chains = 1, cores = 1, iter = 2000)
  expect_equal(coef(fit)$mean[c(1,3,2)],
    fixef(fit_brm)[, "Estimate"], tolerance = 0.1, ignore_attr = TRUE)
  coef(fit,"nu")$mean-ranef(fit_brm)$firm[, "Estimate", ]
  #WIP
})
