alpha <- 311
alpha.10 <- 279
alpha.30 <- 220

omega <- 450.90591
omega.10 <- 413.84352
omega.30 <- 353.663491

# expected values
expectation_b <- 1.4545352
expectation_b.10 <- 1.48864576
expectation_b.30 <- 1.6149018

expectation_inverse_b <- 0.6897226
expectation_inverse_b.10 <- 0.67416786
expectation_inverse_b.30 <- 0.62206025

expectation_inverse_b_2 <- 0.47724691
expectation_inverse_b_2.10 <- 0.45613134
expectation_inverse_b_2.30 <- 0.38871786

expectation_log_b <- 0.37307436
expectation_log_b.10 <- 0.39606933
expectation_log_b.30 <- 0.47699277

test_that("expectation_b", {
  expect_equal(expectation_b(alpha, omega), expectation_b)
  expect_equal(expectation_b(alpha.10, omega.10), expectation_b.10)
  expect_equal(expectation_b(alpha.30, omega.30), expectation_b.30)
})

test_that("expectation_inverse_b", {
  expect_equal(
    expectation_inverse_b(alpha,omega), expectation_inverse_b
  )
  expect_equal(
    expectation_inverse_b(alpha.10,omega.10), expectation_inverse_b.10
  )
  expect_equal(
    expectation_inverse_b(alpha.30,omega.30), expectation_inverse_b.30
  )
})

test_that("expectation_inverse_b_2", {
  expect_equal(
    expectation_inverse_b_2(alpha, omega), expectation_inverse_b_2
  )
  expect_equal(
    expectation_inverse_b_2(alpha.10, omega.10), expectation_inverse_b_2.10
  )
  expect_equal(
    expectation_inverse_b_2(alpha.30, omega.30), expectation_inverse_b_2.30
  )
})

test_that("expectation_log_b", {
  expect_equal(
    expectation_log_b(alpha, omega), expectation_log_b
  )
  expect_equal(
    expectation_log_b(alpha.10, omega.10), expectation_log_b.10
  )
  expect_equal(
    expectation_log_b(alpha.30, omega.30), expectation_log_b.30
  )
})
