test_that("test get_phi", {
  expect_equal(get_phi(-5), 0)
  expect_equal(get_phi(-1.701), 0.0426)
  expect_equal(get_phi(0), 0.3052)
  expect_equal(get_phi(1.702), 0.6950)
  expect_equal(get_phi(5), 0.9574)
  expect_equal(get_phi(6), 1)
})

test_that("test get_zeta", {
  expect_equal(get_zeta(-5), 0)
  expect_equal(get_zeta(-1.7), 0.0189)
  expect_equal(get_zeta(1.7), 0.1138)
  expect_equal(get_zeta(5), 0.0190)
  expect_equal(get_zeta(6), 0)
})

test_that("test get_rho", {
  expect_equal(get_rho(-5), 0)
  expect_equal(get_rho(-1.7), 0.1696)
  expect_equal(get_rho(1.7), 0.5)
  expect_equal(get_rho(5), 0.8303)
  expect_equal(get_rho(6), 1)
})
