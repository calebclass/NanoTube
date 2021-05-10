test_that("gm_mean works", {
  expect_equal(gm_mean(c(2,5,7,11)), 5.26772,
               tolerance = 1e-4)
})
