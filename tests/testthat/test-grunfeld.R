#' @srrstats {G5.0, G5.4, G5.4a, G5.4b, G5.4c} The parameter estimates of
#'   model fitted to a classic Grunfeld data match with brms and plm.

test_that("parameters of the Grunfield model are recovered", {

  skip_if_not(run_extended_tests)
  library("plm")

  data(Grunfeld, package = "plm")
  fit_plm <- plm(inv ~ value + capital, data = Grunfeld,
    index = c("firm", "year"), effect = "individual", model = "within")

  set.seed(1)
  # dynamite defines prior for the intercept based on the mean at the first time
  # point, which differs from the brms, so use dummy intercept instead in both
  Grunfeld$intercept <- 1
  p <- get_priors(obs(inv ~ -1 + intercept + value + capital,
    family = "gaussian", random_intercept = TRUE),
    Grunfeld, "firm", "year")
  # set very vague priors
  p$prior[] <- rep("normal(0, 1000)", nrow(p))
  fit <- dynamite(obs(inv ~ value + capital,
    family = "gaussian", random_intercept = TRUE),
    Grunfeld, "firm", "year", refresh = 0,
    chains = 2, cores = 2, iter = 10000)

  expect_equal(coef(fit_plm), coef(fit)$mean[2:3],
    tolerance = 0.01, ignore_attr = TRUE)

  # Not run, values are stored
  # library(brms)
  # fit_brm <- brm(inv ~ 0 + Intercept + value + capital + (1 | firm),
  #   Grunfeld,
  #   prior = prior(normal(0, 1000), class = "b") +
  #     prior(normal(0, 1000), class = "sd") +
  #     prior(normal(0, 1000), class = "sigma"),
  #     refresh = 0, seed = 1,
  #   chains = 2, cores = 2, iter = 50000, warmup = 5000)
  # brms_est <- round(unname(posterior_summary(fit_brm)[, "Estimate"]), 6)[1:15]
  brms_est <- c(-57.837641, 0.109821, 0.308666, 101.016325, 53.070699,
    -10.084975, 158.146755, -173.565297, 30.013812, -55.072784, 34.483764,
    -8.106796, 0.644707, -28.183244, 50.348913)
  sumr <- as_draws(fit) |> posterior::summarise_draws(
    default_mcse_measures(),
    default_summary_measures())
  for(i in 1:15) {
    expect_equal(sumr$mean[i], brms_est[i],
      tolerance = 10 * sumr$mcse_mean[i], label = sumr$variable[i])
  }
})
