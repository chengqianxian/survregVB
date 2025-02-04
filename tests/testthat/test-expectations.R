# Test data
alpha <- 3
omega <- 2

# Expected values
expected_b <- omega / (alpha - 1)
expected_inv_b <- alpha / omega
expected_inv_b_2 <- (alpha + alpha^2) / omega^2
expected_log_b <- log(omega) - digamma(alpha)

test_that("Test expectation_b", {
  expect_equal(expectation_b(alpha, omega), expected_b)
})

test_that("Test expectation_inverse_b", {
  expect_equal(expectation_inverse_b(alpha, omega), expected_inv_b)
})

test_that("Test expectation_inverse_b_2", {
  expect_equal(expectation_inverse_b_2(alpha, omega), expected_inv_b_2)
})

test_that("Test expectation_log_b", {
  expect_equal(expectation_log_b(alpha, omega), expected_log_b)
})
