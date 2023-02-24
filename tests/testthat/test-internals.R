test_that("formula model topology is correct", {
  f1 <- obs(c(y1, y2) ~ 1 | x, family = "mvgaussian") +
    obs(x ~ 1 + z, family = "gaussian")
  f2 <- obs(z ~ w1, family = "gaussian") +
    obs(c(w1, w2, w3) ~ 1 | y | y, family = "mvgaussian") +
    obs(y ~ 1, family = "gaussian")
  expect_identical(attr(f1, "model_topology"), 2:1)
  expect_identical(attr(f2, "model_topology"), 3:1)
  expect_identical(attr(f2 + f1, "model_topology"), c(3:1, 5:4))
  expect_identical(attr(f1 + f2, "model_topology"), 5:1)
})

test_that("formula parts are correct", {
  f1 <- y ~ x + z * h
  f2 <- ~ x + z * h
  expect_identical(formula_rhs(f1), formula_rhs(f2))
  expect_identical(formula_lhs(f1), quote(y))
  expect_identical(formula_lhs(f2), NULL)
})

test_that("formula incrementation logic is correct", {
  obs_list <- list(
    obs(y ~ 1, family = "gaussian")[[1]],
    obs(y ~ -1 + varying(~1), family = "gaussian")[[1]],
    obs(y ~ x, family = "gaussian")[[1]],
    obs(y ~ -1 + varying(~x), family = "gaussian")[[1]],
    obs(y ~ -1 + x + varying(~z), family = "gaussian")[[1]],
    obs(y ~ x + varying(~ -1 + z), family = "gaussian")[[1]],
    obs(y ~ -1 + varying(~ x + z), family = "gaussian")[[1]],
    obs(y ~ -1 + random(~1), family = "gaussian")[[1]],
    obs(y ~ -1 + random(~x), family = "gaussian")[[1]],
    obs(y ~ -1 + x + random(~z), family = "gaussian")[[1]],
    obs(y ~ x + random(~ -1 + z), family = "gaussian")[[1]],
    obs(y ~ -1 + random(~ x + z), family = "gaussian")[[1]],
    obs(y ~ 1 + varying(~ -1 + x) + random(~ -1 + z), family = "gaussian")[[1]]
  )
  out_list <- list(
    list(
      y ~ 1 + w,
      y ~ 1 + varying(~ -1 + w),
      y ~ 1 + random(~ -1 + w)
    ),
    list(
      y ~ -1 + w + varying(~1),
      y ~ -1 + varying(~ 1 + w),
      y ~ -1 + varying(~1) + random(~ -1 + w)
    ),
    list(
      y ~ 1 + x + w,
      y ~ 1 + x + varying(~ -1 + w),
      y ~ 1 + x + random(~ -1 + w)
    ),
    list(
      y ~ -1 + w + varying(~ 1 + x),
      y ~ -1 + varying(~ 1 + x + w),
      y ~ -1 + varying(~ 1 + x) + random(~ -1 + w)
    ),
    list(
      y ~ -1 + x + w + varying(~ 1 + z),
      y ~ -1 + x + varying(~ 1 + z + w),
      y ~ -1 + x + varying(~ 1 + z) + random(~ -1 + w)
    ),
    list(
      y ~ 1 + x + w + varying(~ -1 + z),
      y ~ 1 + x + varying(~ -1 + z + w),
      y ~ 1 + x + varying(~ -1 + z) + random(~ -1 + w)
    ),
    list(
      y ~ -1 + w + varying(~ 1 + x + z),
      y ~ -1 + varying(~ 1 + x + z + w),
      y ~ -1 + varying(~ 1 + x + z) + random(~ -1 + w)
    ),
    list(
      y ~ -1 + w + random(~1),
      y ~ -1 + varying(~ -1 + w) + random(~1),
      y ~ -1 + random(~ 1 + w)
    ),
    list(
      y ~ -1 + w + random(~ 1 + x),
      y ~ -1 + varying(~ -1 + w) + random(~ 1 + x),
      y ~ -1 + random(~ 1 + x + w)
    ),
    list(
      y ~ -1 + x + w + random(~ 1 + z),
      y ~ -1 + x + varying(~ -1 + w) + random(~ 1 + z),
      y ~ -1 + x + random(~ 1 + z + w)
    ),
    list(
      y ~ 1 + x + w + random(~ -1 + z),
      y ~ 1 + x + varying(~ -1 + w) + random(~ -1 + z),
      y ~ 1 + x + random(~ -1 + z + w)
    ),
    list(
      y ~ -1 + w + random(~ 1 + x + z),
      y ~ -1 + varying(~ -1 + w) + random(~ 1 + x + z),
      y ~ -1 + random(~ 1 + x + z + w)
    ),
    list(
      y ~ 1 + w + varying(~ -1 + x) + random(~ -1 + z),
      y ~ 1 + varying(~ -1 + x + w) + random(~ -1 + z),
      y ~ 1 + varying(~ -1 + x) + random(~ -1 + z + w)
    )
  )
  for (i in seq_along(obs_list)) {
    o <- obs_list[[i]]
    form_list <- list(
      increment_formula(
        formula = o$formula,
        x = "w",
        type = "fixed",
        varying_idx = o$varying,
        fixed_idx = o$fixed,
        random_idx = o$random,
        varying_icpt = o$has_varying_intercept,
        fixed_icpt = o$has_fixed_intercept,
        random_icpt = o$has_random_intercept
      ),
      increment_formula(
        formula = o$formula,
        x = "w",
        type = "varying",
        varying_idx = o$varying,
        fixed_idx = o$fixed,
        random_idx = o$random,
        varying_icpt = o$has_varying_intercept,
        fixed_icpt = o$has_fixed_intercept,
        random_icpt = o$has_random_intercept
      ),
      increment_formula(
        formula = o$formula,
        x = "w",
        type = "random",
        varying_idx = o$varying,
        fixed_idx = o$fixed,
        random_idx = o$random,
        varying_icpt = o$has_varying_intercept,
        fixed_icpt = o$has_fixed_intercept,
        random_icpt = o$has_random_intercept
      )
    )
    expect_equal(form_list[[1L]], out_list[[!!i]][[1L]], ignore_attr = TRUE)
    expect_equal(form_list[[2L]], out_list[[!!i]][[2L]], ignore_attr = TRUE)
    expect_equal(form_list[[3L]], out_list[[!!i]][[3L]], ignore_attr = TRUE)
  }
})

test_that("non-lag term extraction from language objects is correct", {
  expect_identical(find_nonlags(quote(x + z)), c("x", "z"))
  expect_identical(find_nonlags(quote(x + fun(y))), c("x", "y"))
  expect_identical(find_nonlags(quote(x + lag(y))), c("x"))
  expect_identical(find_nonlags(quote(x + fun(y, lag(z)))), c("x", "y"))
})

test_that("internally unsupported families fail", {
  expect_error(
    dynamitefamily("new family"),
    '"new family" is not a supported family\\.'
  )
})

test_that("is.dynamitefamily works", {
  a <- dynamitefamily("gaussian")
  b <- stats::gaussian()
  expect_identical(is.dynamitefamily(a), TRUE)
  expect_identical(is.dynamitefamily(b), FALSE)
})

test_that("paste_rows works correctly", {
  out <- "a\nb\nc\nd\ne\nf"
  expect_identical(
    paste_rows(c("a", "b", "c"), c("d", "e"), "f", .parse = TRUE),
    out
  )
  expect_identical(
    paste_rows(c("a", "b", "c"), c("d", "e"), "f", .parse = FALSE),
    out
  )
})

test_that("message_ outputs", {
  expect_message(message_("This is a message"), "This is a message")
})
