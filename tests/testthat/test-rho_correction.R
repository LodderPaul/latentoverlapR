test_that("rho_correction returns a result", {
  mod <- reliability_model(relA = .80, relB = .85, kA = 6, kB = 6)

  out <- rho_correction(
    r_obs = .60,
    type = "continuous",
    measurement_model = mod,
    cA = .30,
    cB = .30,
    k_overlap = 1,
    ci = "none"
  )

  expect_s3_class(out, "latentoverlap_result")
  expect_true(is.numeric(out$rho_hat))
})
