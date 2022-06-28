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
    obs(y ~ -1 + varying(~x + z), family = "gaussian")[[1]]
  )
  out_list <- list(
    list(
      y ~ 1 + w,
      y ~ 1 + varying(~-1 + w)
    ),
    list(
      y ~ -1 + w + varying(~1),
      y ~ -1 + varying(~1 + w)
    ),
    list(
      y ~ x + w,
      y ~ x + varying(~-1 + w)
    ),
    list(
      y ~ -1 + w + varying(~1 + x),
      y ~ -1 + varying(~1 + x + w)
    ),
    list(
      y ~ x - 1 + w + varying(~1 + z),
      y ~ x - 1 + varying(~1 + z + w)
    ),
    list(
      y ~ x + w + varying(~ -1 + z),
      y ~ x + varying(~ -1 + z + w)
    ),
    list(
      y ~ -1 + w + varying(~1 + x + z),
      y ~ -1 + varying(~1 + x + z + w)
    )
  )
  for (i in seq_along(obs_list)) {
    o <- obs_list[[i]]
    form_list <- list(
      increment_formula(o$formula, "w", type = "fixed", o$varying,
                        o$has_varying_intercept, o$has_fixed_intercept),
      increment_formula(o$formula, "w", type = "varying", o$varying,
                        o$has_varying_intercept, o$has_fixed_intercept)
    )
    expect_equal(form_list[[1]], out_list[[i]][[1]], ignore_attr = TRUE)
    expect_equal(form_list[[2]], out_list[[i]][[2]], ignore_attr = TRUE)
  }
  obs_aux <- aux(numeric(d) ~ 1 + x)
  obs_aux_inc <- aux(numeric(d) ~ 1 + x + lag(y))
  expect_equal(
    increment_formula_deterministic(obs_aux[[1]]$formula, "lag(y)"),
    obs_aux_inc[[1]]$formula,
    ignore_attr = TRUE)
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
