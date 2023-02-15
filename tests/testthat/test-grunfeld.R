#' @srrstats {G5.0, G5.4, G5.4a, G5.4b, G5.4c} The parameter estimates of
#'   model fitted to a classic Grunfeld data match with brms and plm.
run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("parameters of the Grunfield model are recovered", {
  skip_if_not(run_extended_tests)
  if (require("plm", quietly = TRUE)) {

    data(Grunfeld, package = "plm")
    fit_plm <- plm(inv ~ value + capital,
      data = Grunfeld,
      index = c("firm", "year"), effect = "individual", model = "within"
    )

    set.seed(1)
    # dynamite defines prior for the intercept based on the mean at the first time
    # point, which differs from the brms, so use dummy intercept instead in both
    Grunfeld$intercept <- 1
    f <- obs(inv ~ -1 + intercept + value + capital + random(~-1 + intercept),
      family = "gaussian") + random_spec(noncentered = FALSE)
    p <- get_priors(f,
      Grunfeld, time = "year", group = "firm"
    )
    # set very vague priors
    p$prior[] <- rep("normal(0, 1000)", nrow(p))
    fit <- dynamite(f,
      Grunfeld, time = "year", group = "firm",
      refresh = 0, seed = 1,
      chains = 2, cores = 2, iter = 20000, warmup = 1000
    )

    expect_equal(coef(fit_plm), coef(fit)$mean[2:3],
      tolerance = 0.01, ignore_attr = TRUE
    )

    # Not run, values are stored
    # library(brms)
    # fit_brm <- brm(inv ~ 0 + Intercept + value + capital + (1 | firm),
    #   Grunfeld,
    #   prior = prior(normal(0, 1000), class = "b") +
    #     prior(normal(0, 1000), class = "sd") +
    #     prior(normal(0, 1000), class = "sigma"),
    #     refresh = 0, seed = 1,
    #   chains = 2, cores = 2, iter = 150000, warmup = 5000,
    #   control = list(adapt_delta = 0.9))
    # brms_est <- round(unname(posterior_summary(fit_brm)[, "Estimate"]), 6)[1:15]
    brms_est <- c(
      -57.917705, 0.109846, 0.308349, 100.941716, 53.086506, -9.896139,
      158.194407, -173.491487, 29.99692, -54.875308, 34.459838, -7.931004,
      0.646235, -28.227079, 50.50187
    )
    # reorder parameters to match dynamite
    brms_est <- brms_est[c(3, 1, 2, 6:15, 5, 4)]
    sumr <- as_draws(fit) |> posterior::summarise_draws(
      posterior::default_mcse_measures(),
      posterior::default_summary_measures()
    )
    for (i in 1:15) {
      expect_equal(sumr$mean[i], brms_est[i],
        tolerance = 100 * sumr$mcse_mean[i], label = sumr$variable[i]
      )
    }
  }
})
