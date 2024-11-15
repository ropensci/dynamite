# Working directory needs to be /vignettes

rebuild_vignettes <- function() {
  files <- c(
    "dynamite_custom",
    "dynamite_simulation"
  )
  for (i in seq_along(files)) {
    e <- new.env()
    knitr::knit(
      input = paste0(files[i], ".qmd.orig"),
      output = paste0(files[i], ".qmd"),
      envir = e
    )
  }
}
