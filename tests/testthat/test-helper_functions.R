test_that("add_standardized_indicators() adds norm_ and rank_ columns", {
  dat <- tibble::tibble(id = 1:4, raw_a = c(10, 20, 30, 40))
  out <- add_standardized_indicators(dat, prefix = "raw_a")

  expect_true(all(c("norm_a", "rank_a") %in% names(out)))
  expect_equal(out$norm_a, c(0, 1 / 3, 2 / 3, 1))
  expect_equal(out$rank_a, c(25, 50, 75, 100))
  expect_equal(out$raw_a, dat$raw_a)
})

test_that("add_standardized_indicators(invert = TRUE) flips the direction", {
  dat <- tibble::tibble(raw_a = c(10, 20, 30, 40))
  out <- add_standardized_indicators(dat, prefix = "raw_a", invert = TRUE)

  expect_equal(out$norm_a, c(1, 2 / 3, 1 / 3, 0))
  expect_equal(out$rank_a, c(75, 50, 25, 0))
})

test_that("add_standardized_indicators() gives NA norm for an all-NA column", {
  dat <- tibble::tibble(raw_a = c(NA_real_, NA_real_))
  out <- add_standardized_indicators(dat, prefix = "raw_a")

  expect_true(all(is.na(out$norm_a)))
})
