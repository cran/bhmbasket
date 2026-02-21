# Tests for getMuVar -----------------------------------------------------------

# ------------------------------------------------------------------
# Test: correct mu_var for simple reference case
# Input:
#   - response_rate = 0.5, tau_scale = 1, n_worth = 2
# Behaviour:
#   - getMuVar() implements: mu_var = (n_worth * p * (1 - p))^-1 - tau_scale^2
# Expectations:
#   - Result equals manual calculation, is a finite double.
# Why:
#   - Sanity-checks the core formula and numeric type for a simple case.
# ------------------------------------------------------------------
test_that("mu_var calculation is correct for simple case", {
  
  # Manual calculation: (n_worth * p * (1 - p))^-1 - tau_scale^2
  expected <- (2 * 0.5 * (1 - 0.5))^-1 - 1^2
  
  result <- getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 2)
  
  expect_equal(result, expected)
  
  expect_type(result, "double")
  
  expect_true(is.finite(result))
  
})

# ------------------------------------------------------------------
# Test: response_rate must lie strictly between 0 and 1
# Input:
#   - response_rate = -0.1 and 1.5, tau_scale = 1
# Behaviour:
#   - Function should reject probabilities outside (0,1).
# Expectations:
#   - Errors for each invalid response_rate.
# Why:
#   - The variance formula assumes a valid Bernoulli probability.
# ------------------------------------------------------------------
test_that("error if response_rate out of bounds", {
  
  expect_error(getMuVar(response_rate = -0.1, tau_scale = 1))
  
  expect_error(getMuVar(response_rate = 1.5, tau_scale = 1))
  
})

# ------------------------------------------------------------------
# Test: tau_scale must be non-negative
# Input:
#   - response_rate = 0.5, tau_scale = -1
# Behaviour:
#   - Negative scale parameter is invalid and must trigger an error.
# Expectations:
#   - Error for tau_scale < 0.
# Why:
#   - tau_scale represents a scale for a prior; negative values are not meaningful.
# ------------------------------------------------------------------
test_that("error if tau_scale negative", {
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = -1))
  
})

# ------------------------------------------------------------------
# Test: n_worth must be a positive integer
# Input:
#   - n_worth = 0 and 1.5
# Behaviour:
#   - Non-positive or non-integer n_worth should fail validation.
# Expectations:
#   - Errors in both cases.
# Why:
#   - n_worth describes an “effective sample size”; must be >= 1 and integer.
# ------------------------------------------------------------------
test_that("error if n_worth not positive integer", {
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 0))
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 1.5))
  
})


# Tests for getPriorParametersBerry --------------------------------------------

# ------------------------------------------------------------------
# Test: getPriorParametersBerry returns correctly structured object
# Input:
#   - target_rates = c(0.2, 0.8), tau_scale = 1, n_worth = 2
# Behaviour:
#   - Function builds a prior_parameters_list for "berry" with fields
#     mu_mean, mu_sd, tau_scale.
# Expectations:
#   - Class "prior_parameters_list", name "berry", inner list with
#     mu_mean, mu_sd, tau_scale.
# Why:
#   - Verifies the returned container type and its internal naming.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersBerry(target_rates = c(0.2, 0.8), tau_scale = 1, n_worth = 2)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
  
  # Check inner structure
  expect_true(is.list(result$berry))
  expect_named(result$berry, c("mu_mean", "mu_sd", "tau_scale"))
})

# ------------------------------------------------------------------
# Test: mu_sd is computed from mu_var as defined by max-variance target rate
# Input:
#   - target_rates = c(0.2, 0.8), tau_scale = 1, n_worth = 2
# Behaviour:
#   - Pick target_rate furthest from 0.5, compute
#     mu_var = (n_worth * p * (1 - p))^-1 - tau_scale^2, then mu_sd = sqrt(mu_var).
# Expectations:
#   - result$berry$mu_sd equals the manually computed value.
# Why:
#   - Confirms the link between target_rates, mu_var and mu_sd.
# ------------------------------------------------------------------
test_that("mu_sd is computed correctly for simple case", {
  target_rates <- c(0.2, 0.8)
  tau_scale <- 1
  n_worth <- 2
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_var <- (n_worth * target_rate_max_var * (1 - target_rate_max_var))^-1 - tau_scale^2
  expected_mu_sd <- sqrt(expected_mu_var)
  
  result <- getPriorParametersBerry(target_rates, tau_scale, n_worth)
  expect_equal(result$berry$mu_sd, expected_mu_sd)
})

# ------------------------------------------------------------------
# Test: getPriorParametersBerry errors if mu_var <= 0
# Input:
#   - target_rates = 0.5, tau_scale = 100, n_worth = 1
# Behaviour:
#   - Very large tau_scale leads to mu_var <= 0, which is invalid.
# Expectations:
#   - Error is thrown.
# Why:
#   - A non-positive mu_var cannot produce a real standard deviation.
# ------------------------------------------------------------------
test_that("error if mu_var <= 0", {
  # Large tau_scale will make mu_var negative
  expect_error(
    getPriorParametersBerry(target_rates = c(0.5), tau_scale = 100, n_worth = 1)
  )
})

# ------------------------------------------------------------------
# Test: single target_rate works and returns a berry prior
# Input:
#   - target_rates = 0.5, tau_scale = 1, n_worth = 1
# Behaviour:
#   - Function should still construct a valid prior_parameters_list.
# Expectations:
#   - Class "prior_parameters_list", name "berry".
# Why:
#   - Ensures single-cohort use is supported.
# ------------------------------------------------------------------
test_that("single target_rate works", {
  result <- getPriorParametersBerry(target_rates = c(0.5), tau_scale = 1, n_worth = 1)
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
})

# Tests for setPriorParametersBerry --------------------------------------------

# ------------------------------------------------------------------
# Test: setPriorParametersBerry creates correct structure
# Input:
#   - mu_mean = 0.1, mu_sd = 0.5, tau_scale = 1
# Behaviour:
#   - Function wraps values into a "berry" prior_parameters_list.
# Expectations:
#   - Class "prior_parameters_list", name "berry",
#     berry list contains mu_mean, mu_sd, tau_scale.
# Why:
#   - Tests manual prior specification path for Berry.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersBerry(
    mu_mean = 0.1,
    mu_sd = 0.5,
    tau_scale = 1
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
  expect_true(is.list(result$berry))
  expect_named(result$berry, c("mu_mean", "mu_sd", "tau_scale"))
})

# ------------------------------------------------------------------
# Test: mu_sd must be strictly positive
# Input:
#   - mu_sd = 0
# Behaviour:
#   - Function should reject non-positive standard deviations.
# Expectations:
#   - Error mentioning "mu_sd".
# Why:
#   - Standard deviations must be > 0.
# ------------------------------------------------------------------
test_that("error if mu_sd is non-positive", {
  expect_error(
    setPriorParametersBerry(mu_mean = 0.1, mu_sd = 0, tau_scale = 1),
    "mu_sd"
  )
})

# ------------------------------------------------------------------
# Test: tau_scale must be strictly positive
# Input:
#   - tau_scale = 0
# Behaviour:
#   - Function should reject non-positive scale values.
# Expectations:
#   - Error mentioning "tau_scale".
# Why:
#   - Scale parameter for tau must be > 0.
# ------------------------------------------------------------------
test_that("error if tau_scale is non-positive", {
  expect_error(
    setPriorParametersBerry(mu_mean = 0.1, mu_sd = 0.5, tau_scale = 0),
    "tau_scale"
  )
})

# ------------------------------------------------------------------
# Test: mu_mean must be numeric
# Input:
#   - mu_mean = "a" (character)
# Behaviour:
#   - Non-numeric mu_mean should be rejected.
# Expectations:
#   - Error mentioning "mu_mean".
# Why:
#   - Prior mean is a numeric parameter on the link scale.
# ------------------------------------------------------------------
test_that("error if mu_mean is not numeric", {
  expect_error(
    setPriorParametersBerry(mu_mean = "a", mu_sd = 0.5, tau_scale = 1),
    "mu_mean"
  )
})

# Tests for getPriorParametersExNeX --------------------------------------------

# ------------------------------------------------------------------
# Test: getPriorParametersExNex returns correctly structured object
# Input:
#   - target_rates = c(0.3, 0.9), tau_scale = 1, n_worth = 2, w_j = 0.5
# Behaviour:
#   - Builds an "exnex" prior with fields mu_mean, mu_sd, tau_scale,
#     mu_j, tau_j, w_j.
# Expectations:
#   - Class "prior_parameters_list", name "exnex",
#     exnex list has all expected fields.
# Why:
#   - Checks the container structure for the ExNex prior definition.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersExNex(target_rates = c(0.3, 0.9), tau_scale = 1, n_worth = 2, w_j = 0.5)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex")
  
  # Check inner structure
  expect_true(is.list(result$exnex))
  expect_named(result$exnex, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})

# ------------------------------------------------------------------
# Test: ExNex prior uses logit transform for mu_mean and mu_j
# Input:
#   - target_rates = 0.8, tau_scale = 1, n_worth = 2, w_j = 0.5
# Behaviour:
#   - mu_mean = logit(p_max_var), mu_j = logit(target_rates).
# Expectations:
#   - result$exnex$mu_mean and mu_j match manual logit calculations.
# Why:
#   - Confirms link-scale parameterisation is correctly applied.
# ------------------------------------------------------------------
test_that("mu_mean and mu_j use logit transform", {
  target_rates <- 0.8
  tau_scale <- 1
  n_worth <- 2
  w_j <- 0.5
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_mean <- logit(target_rate_max_var)
  expected_mu_j <- logit(target_rates)
  result <- getPriorParametersExNex(target_rates, tau_scale, n_worth, w_j)
  
  expect_equal(result$exnex$mu_mean, expected_mu_mean)
  expect_equal(result$exnex$mu_j, expected_mu_j)
})

# ------------------------------------------------------------------
# Test: getPriorParametersExNex errors if mu_var <= 0
# Input:
#   - target_rates = 0.5, tau_scale = 100, n_worth = 1
# Behaviour:
#   - Large tau_scale makes mu_var non-positive.
# Expectations:
#   - Error thrown.
# Why:
#   - Same reason as Berry: must have positive variance for mu.
# ------------------------------------------------------------------
test_that("error if mu_var <= 0", {
  # Large tau_scale will make mu_var negative
  expect_error(
    getPriorParametersExNex(target_rates = c(0.5), tau_scale = 100, n_worth = 1)
  )
})

# ------------------------------------------------------------------
# Test: mu_sd and tau_j computations follow theoretical formulas
# Input:
#   - target_rates = 0.8, tau_scale = 1, n_worth = 2
# Behaviour:
#   - mu_sd from max-variance target, tau_j from per-cohort variance formula.
# Expectations:
#   - result$exnex$mu_sd and tau_j match manual calculations.
# Why:
#   - Ensures ExNex prior scale parameters are derived correctly.
# ------------------------------------------------------------------
test_that("mu_sd and tau_j are computed correctly", {
  target_rates <- 0.8
  tau_scale <- 1
  n_worth <- 2
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_var <- (n_worth * target_rate_max_var * (1 - target_rate_max_var))^-1 - tau_scale^2
  expected_mu_sd <- sqrt(expected_mu_var)
  
  expected_tau_j <- sqrt((n_worth * target_rates * (1 - target_rates))^-1)
  
  result <- getPriorParametersExNex(target_rates, tau_scale, n_worth, w_j = 0.5)
  
  expect_equal(result$exnex$mu_sd, expected_mu_sd)
  expect_equal(result$exnex$tau_j, expected_tau_j)
})

# Tests for setPriorParametersExNeX --------------------------------------------

# ------------------------------------------------------------------
# Test: setPriorParametersExNeX creates full exnex structure
# Input:
#   - Vectors for mu_mean, mu_sd, tau_scale, mu_j, tau_j, w_j.
# Behaviour:
#   - Returns a prior_parameters_list named "exnex" with all components.
# Expectations:
#   - Class "prior_parameters_list", name "exnex",
#     exnex element has mu_mean, mu_sd, tau_scale, mu_j, tau_j, w_j.
# Why:
#   - Validates the manual ExNex prior construction path.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersExNex(
    mu_mean = c(0.1, 0.2),
    mu_sd = c(0.5, 0.6),
    tau_scale = 1,
    mu_j = c(0.1, 0.2),
    tau_j = c(0.3, 0.4),
    w_j = c(0.3, 0.3, 0.4)
  )

  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex")
  expect_true(is.list(result$exnex))
  expect_named(result$exnex, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})

# ------------------------------------------------------------------
# Test: w_j length rules when mu_mean length = 1
# Input:
#   - mu_mean length 1, w_j length 3
# Behaviour:
#   - For a single mu_mean, w_j length 3 is not allowed.
# Expectations:
#   - Error.
# Why:
#   - Tests internal consistency checks between mixing weights and mean length.
# ------------------------------------------------------------------
test_that("error if w_j length is invalid when mu_mean length = 1", {
  expect_error(
    setPriorParametersExNex(
      mu_mean = c(0.1),
      mu_sd = c(0.5),
      tau_scale = 1,
      mu_j = c(0.1),
      tau_j = c(0.3),
      w_j = c(0.3, 0.3, 0.4)  # length 3 is invalid for mu_mean length 1
    )
  )
})

# Tests for getPriorParametersExNeXAdj -----------------------------------------

# ------------------------------------------------------------------
# Test: getPriorParametersExNeXAdj returns exnex_adj with centered mu's
# Input:
#   - target_rates = c(0.2, 0.3), tau_scale = 1, n_worth = 2, w_j = 0.5
# Behaviour:
#   - Returns "exnex_adj" prior with mu_mean and mu_j shifted to 0.
# Expectations:
#   - Class "prior_parameters_list", name "exnex_adj",
#     mu_mean = 0, mu_j = 0-vector.
# Why:
#   - Confirms adjusted ExNex prior is properly centered.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersExNexAdj(
    target_rates = c(0.2, 0.3),
    tau_scale = 1,
    n_worth = 2,
    w_j = 0.5
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex_adj")
  expect_true(is.list(result$exnex_adj))
  expect_named(result$exnex_adj, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
  
  # Check adjusted values
  expect_equal(result$exnex_adj$mu_mean, 0)
  expect_equal(result$exnex_adj$mu_j, rep(0, length(c(0.2, 0.3))))
})

# Tests for setPriorParametersExNeXAdj -----------------------------------------

# ------------------------------------------------------------------
# Test: setPriorParametersExNeXAdj creates correct exnex_adj structure
# Input:
#   - Vectors for mu_mean, mu_sd, tau_scale, mu_j, tau_j, w_j.
# Behaviour:
#   - Returns "exnex_adj" prior_parameters_list with all fields.
# Expectations:
#   - Class "prior_parameters_list", name "exnex_adj",
#     exnex_adj has mu_mean, mu_sd, tau_scale, mu_j, tau_j, w_j.
# Why:
#   - Checks manual specification of adjusted ExNex prior.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersExNexAdj(
    mu_mean = c(0.1, 0.2),
    mu_sd = c(0.5, 0.6),
    tau_scale = 1,
    mu_j = c(0.1, 0.2),
    tau_j = c(0.3, 0.4),
    w_j = c(0.3, 0.3, 0.4)
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex_adj")
  expect_true(is.list(result$exnex_adj))
  expect_named(result$exnex_adj, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})

# ------------------------------------------------------------------
# Test: acceptable w_j lengths when mu_mean length = 1
# Input:
#   - mu_mean length 1, w_j length 1 (valid).
# Behaviour:
#   - Function should accept such configuration silently.
# Expectations:
#   - No error raised.
# Why:
#   - Sanity-checks "short" w_j is allowed in simple case.
# ------------------------------------------------------------------
test_that("'w_j' length rule works when mu_mean length = 1", {
  # Valid cases: w_j length 1 or 2
  expect_silent(
    setPriorParametersExNexAdj(
      mu_mean = 0.1,
      mu_sd = 0.5,
      tau_scale = 1,
      mu_j = c(0.1),
      tau_j = c(0.3),
      w_j = c(0.5)  # length 1
    )
  )
})

# Tests for getPriorParametersPooled --------------------------------------------

# ------------------------------------------------------------------
# Test: getPriorParametersPooled returns pooled Beta parameters
# Input:
#   - target_rates = c(0.2, 0.4, 0.6), n_worth = 2
# Behaviour:
#   - Selects a target_rate (internally defined rule) and computes
#     a = p * n_worth, b = (1 - p) * n_worth.
# Expectations:
#   - Class "prior_parameters_list", name "pooled", with positive a,b.
# Why:
#   - Validates the pooled Beta( a, b ) construction.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersPooled(
    target_rates = c(0.2, 0.4, 0.6),
    n_worth = 2
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "pooled")
  expect_true(is.list(result$pooled))
  expect_named(result$pooled, c("a", "b"))
  expect_true(result$pooled$a > 0)
  expect_true(result$pooled$b > 0)
})

# ------------------------------------------------------------------
# Test: pooled prior requires valid target_rates in (0,1)
# Input:
#   - target_rates including 0 or 1.
# Behaviour:
#   - Must error on invalid probabilities.
# Expectations:
#   - Errors mentioning "target_rates".
# Why:
#   - Beta prior requires positive support away from 0 and 1.
# ------------------------------------------------------------------
test_that("error if target_rates contains invalid values", {
  expect_error(getPriorParametersPooled(target_rates = c(0, 0.5), n_worth = 1), "target_rates")
  expect_error(getPriorParametersPooled(target_rates = c(1, 0.5), n_worth = 1), "target_rates")
})

# ------------------------------------------------------------------
# Test: n_worth must be positive for pooled prior
# Input:
#   - n_worth = 0
# Behaviour:
#   - Function should reject non-positive n_worth.
# Expectations:
#   - Error mentioning "n_worth".
# Why:
#   - n_worth again represents effective sample size.
# ------------------------------------------------------------------
test_that("error if n_worth is non-positive", {
  expect_error(getPriorParametersPooled(target_rates = c(0.3, 0.4), n_worth = 0), "n_worth")
})

# ------------------------------------------------------------------
# Test: pooled prior selects target_rate closest to 0.5
# Input:
#   - target_rates = c(0.1, 0.4, 0.8), n_worth = 3
# Behaviour:
#   - Among candidates, p = 0.4 is closest to 0.5, so a = 0.4*3, b = 0.6*3.
# Expectations:
#   - result$pooled$a == 0.4*3, result$pooled$b == 0.6*3.
# Why:
#   - Checks selection logic for the central rate.
# ------------------------------------------------------------------
test_that("selects target_rate closest to 0.5", {
  result <- getPriorParametersPooled(target_rates = c(0.1, 0.4, 0.8), n_worth = 3)
  expect_equal(result$pooled$a, 0.4 * 3)
  expect_equal(result$pooled$b, (1 - 0.4) * 3)
})

# Tests for setPriorParametersPooled -------------------------------------------

# ------------------------------------------------------------------
# Test: setPriorParametersPooled builds pooled Beta structure
# Input:
#   - a = 2, b = 3
# Behaviour:
#   - Wraps into "pooled" prior_parameters_list.
# Expectations:
#   - Class "prior_parameters_list", name "pooled", a and b preserved.
# Why:
#   - Tests manual specification path for pooled prior.
# ------------------------------------------------------------------
test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersPooled(a = 2, b = 3)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "pooled")
  expect_true(is.list(result$pooled))
  expect_named(result$pooled, c("a", "b"))
  expect_equal(result$pooled$a, 2)
  expect_equal(result$pooled$b, 3)
})

# Tests for getPriorParametersStratified ---------------------------------------

# ------------------------------------------------------------------
# Test: stratified prior computes a_j and b_j per target rate
# Input:
#   - target_rates = c(0.2, 0.4, 0.6), n_worth = 2
# Behaviour:
#   - a_j = p * n_worth, b_j = (1 - p) * n_worth elementwise.
# Expectations:
#   - result$stratified$a_j and b_j match manual multiplication.
# Why:
#   - Verifies Beta-prior parameters for each stratum.
# ------------------------------------------------------------------
test_that("computes a_j and b_j correctly for multiple target rates", {
  target_rates <- c(0.2, 0.4, 0.6)
  n_worth <- 2
  result <- getPriorParametersStratified(target_rates, n_worth)
  
  expect_equal(result$stratified$a_j, target_rates * n_worth)
  expect_equal(result$stratified$b_j, (1 - target_rates) * n_worth)
})

# ------------------------------------------------------------------
# Test: stratified prior handles a single target rate
# Input:
#   - target_rates = 0.3, n_worth = 5
# Behaviour:
#   - Computes scalar a_j and b_j.
# Expectations:
#   - a_j = 0.3*5, b_j = 0.7*5.
# Why:
#   - Ensures single-stratum usage is supported.
# ------------------------------------------------------------------
test_that("handles single target rate correctly", {
  result <- getPriorParametersStratified(target_rates = 0.3, n_worth = 5)
  expect_equal(result$stratified$a_j, 0.3 * 5)
  expect_equal(result$stratified$b_j, (1 - 0.3) * 5)
})

# ------------------------------------------------------------------
# Test: stratified result has stable structure
# Input:
#   - target_rates = c(0.25, 0.75), n_worth = 3
# Behaviour:
#   - Returns "stratified" prior with a_j and b_j fields.
# Expectations:
#   - Class "prior_parameters_list", name "stratified", fields a_j,b_j.
# Why:
#   - Confirms naming and type are consistent.
# ------------------------------------------------------------------
test_that("result structure is consistent", {
  result <- getPriorParametersStratified(c(0.25, 0.75), 3)
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "stratified")
  expect_named(result$stratified, c("a_j", "b_j"))
})

# Tests for setPriorParametersStratified ---------------------------------------

# ------------------------------------------------------------------
# Test: setPriorParametersStratified preserves a_j,b_j and structure
# Input:
#   - a_j = c(1,2,3), b_j = c(4,5,6)
# Behaviour:
#   - Returns "stratified" prior with these vectors.
# Expectations:
#   - a_j,b_j unchanged; class prior_parameters_list.
# Why:
#   - Checks manual specification path for stratified Beta priors.
# ------------------------------------------------------------------
test_that("returns correct structure and preserves input values", {
  a_j <- c(1, 2, 3)
  b_j <- c(4, 5, 6)
  result <- setPriorParametersStratified(a_j, b_j)
  
  expect_equal(result$stratified$a_j, a_j)
  expect_equal(result$stratified$b_j, b_j)
  expect_s3_class(result, "prior_parameters_list")
})

# ------------------------------------------------------------------
# Test: setPriorParametersStratified works with single-element vectors
# Input:
#   - a_j = 10, b_j = 20
# Behaviour:
#   - Returns scalar a_j,b_j under "stratified".
# Expectations:
#   - Values preserved.
# Why:
#   - Confirms scalar case is supported smoothly.
# ------------------------------------------------------------------
test_that("works with single-element vectors", {
  result <- setPriorParametersStratified(a_j = 10, b_j = 20)
  expect_equal(result$stratified$a_j, 10)
  expect_equal(result$stratified$b_j, 20)
})

# ------------------------------------------------------------------
# Test: length consistency for a_j and b_j
# Input:
#   - a_j length 2, b_j length 1
# Behaviour:
#   - Function should error on mismatched lengths.
# Expectations:
#   - Error mentioning "a_j and b_j".
# Why:
#   - Each stratum must have both a and b.
# ------------------------------------------------------------------
test_that("length consistency check works", {
  expect_error(setPriorParametersStratified(a_j = c(1, 2), b_j = c(3)), "a_j and b_j")
})

# Tests for getPriorParameters -------------------------------------------------

# ------------------------------------------------------------------
# Test: getPriorParameters returns combined prior_parameters_list
# Input:
#   - method_names = c("berry", "pooled", "stratified", "exnex", "exnex_adj"),
#     plus reasonable target_rates, n_worth, tau_scale, w_j.
# Behaviour:
#   - Constructs all requested method priors and merges them.
# Expectations:
#   - Class "prior_parameters_list", 5 elements, names sorted as in implementation.
# Why:
#   - Validates the high-level prior factory for multiple methods.
# ------------------------------------------------------------------
test_that("valid input returns a prior_parameters_list with correct names and class", {
  result <- getPriorParameters(
    method_names = c("berry", "pooled", "stratified", "exnex", "exnex_adj"),
    target_rates = c(0.2, 0.1, 0.3, 0.2, 0.4, 0.6),
    n_worth = 2,
    tau_scale = 1,
    w_j = 0.5
  )
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, c("berry", "exnex", "exnex_adj", "pooled", "stratified"))
  expect_type(result, "list")
  expect_length(result, 5)

})

# ------------------------------------------------------------------
# Test: getPriorParameters errors on invalid method_names
# Input:
#   - method_names = "invalid"
# Behaviour:
#   - Function must reject unknown methods.
# Expectations:
#   - Error.
# Why:
#   - Guards against typos or unsupported priors.
# ------------------------------------------------------------------
test_that("invalid method_names throws error", {
  expect_error(getPriorParameters(
    method_names = "invalid",
    target_rates = c(0.2, 0.3)
  ))
})

# ------------------------------------------------------------------
# Test: getPriorParameters errors on invalid target_rates
# Input:
#   - target_rates outside [0,1]
# Behaviour:
#   - Should reject probabilities <=0 or >=1.
# Expectations:
#   - Error.
# Why:
#   - All underlying constructions assume valid rates.
# ------------------------------------------------------------------
test_that("invalid target_rates throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(-0.1, 1.2)
  ))
})

# ------------------------------------------------------------------
# Test: getPriorParameters errors on invalid tau_scale
# Input:
#   - tau_scale = -1
# Behaviour:
#   - Negative tau_scale should fail.
# Expectations:
#   - Error.
# Why:
#   - Same constraint as lower-level prior builders.
# ------------------------------------------------------------------
test_that("invalid tau_scale throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    tau_scale = -1
  ))
})

# ------------------------------------------------------------------
# Test: getPriorParameters errors on invalid n_worth
# Input:
#   - n_worth = 0
# Behaviour:
#   - Must reject non-positive effective sample size.
# Expectations:
#   - Error.
# Why:
#   - Consistency with other n_worth checks.
# ------------------------------------------------------------------
test_that("invalid n_worth throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    n_worth = 0
  ))
})

# ------------------------------------------------------------------
# Test: getPriorParameters errors on invalid w_j
# Input:
#   - w_j = 2 (outside [0,1] or invalid)
# Behaviour:
#   - Weight must be in acceptable range.
# Expectations:
#   - Error.
# Why:
#   - Mixture weights need to be interpretable probabilities.
# ------------------------------------------------------------------
test_that("invalid w_j throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    w_j = 2
  ))
})

# Tests for combinePriorParameters ---------------------------------------------

# ------------------------------------------------------------------
# Test: combinePriorParameters merges multiple prior_parameter objects
# Input:
#   - prior_parameters_stratified + prior_parameters_berry
# Behaviour:
#   - Flattens and merges, returning a single prior_parameters_list.
# Expectations:
#   - Class "prior_parameters_list", length 2, names are union of input names,
#     each element is a list.
# Why:
#   - Confirms high-level combination for multi-method priors.
# ------------------------------------------------------------------
test_that("combinePriorParameters returns correct structure with real objects", {
  prior_parameters_stratified <- setPriorParametersStratified(c(1, 2), c(3, 4))
  prior_parameters_berry      <- setPriorParametersBerry(1, 1, 2)
  
  result <- combinePriorParameters(list(prior_parameters_berry, prior_parameters_stratified))
  
  # Check class and type
  expect_s3_class(result, "prior_parameters_list")
  expect_type(result, "list")
  
  # Names should match method names from input
  expect_named(result, sort(c(names(prior_parameters_berry), names(prior_parameters_stratified))))
  
  # Length should equal number of input elements
  expect_length(result, 2)
  
  # Each element should be a list (inner structure)
  expect_true(all(vapply(result, is.list, logical(1))))
})

# ------------------------------------------------------------------
# Test: combinePriorParameters sorts method names alphabetically
# Input:
#   - [stratified, berry] in that order
# Behaviour:
#   - Resulting list names should be sorted.
# Expectations:
#   - names(result) == sort of input method names.
# Why:
#   - Ensures deterministic ordering in combined prior.
# ------------------------------------------------------------------
test_that("combinePriorParameters sorts names alphabetically", {
  prior_parameters_stratified <- setPriorParametersStratified(c(1, 2), c(3, 4))
  prior_parameters_berry      <- setPriorParametersBerry(0, 1, 2)
  
  result <- combinePriorParameters(list(prior_parameters_stratified, prior_parameters_berry))
  
  expect_equal(names(result), sort(c(names(prior_parameters_stratified), names(prior_parameters_berry))))
})

# ------------------------------------------------------------------
# Test: combinePriorParameters errors on duplicate method names
# Input:
#   - Two separate prior_parameters_berry objects.
# Behaviour:
#   - Duplicate "berry" name should be rejected.
# Expectations:
#   - Error.
# Why:
#   - Each method may appear only once in the combined object.
# ------------------------------------------------------------------
test_that("combinePriorParameters errors on duplicate method names", {
  prior_parameters_berry1 <- setPriorParametersBerry(0, 1, 2)
  prior_parameters_berry2 <- setPriorParametersBerry(0, 1, 2)
  
  expect_error(
    combinePriorParameters(list(prior_parameters_berry1, prior_parameters_berry2))
  )
})

# ------------------------------------------------------------------
# Test: combinePriorParameters errors when input is not a list
# Input:
#   - "not_a_list"
# Behaviour:
#   - Type check should fail for non-list input.
# Expectations:
#   - Error.
# Why:
#   - Enforces input contract: a list of prior_parameters_list objects.
# ------------------------------------------------------------------
test_that("combinePriorParameters errors if input is not a list", {
  expect_error(combinePriorParameters("not_a_list"))
})

# ------------------------------------------------------------------
# Test: combinePriorParameters errors if any element has wrong class
# Input:
#   - One prior_parameters_berry, one plain list(dummy = TRUE)
# Behaviour:
#   - Mixed-class list must be rejected.
# Expectations:
#   - Error.
# Why:
#   - Every element must be a proper prior_parameters_list.
# ------------------------------------------------------------------
test_that("combinePriorParameters errors if any element is not prior_parameters_list", {
  prior_parameters_berry <- setPriorParametersBerry(0, 1, 2)
  bad_input <- list(prior_parameters_berry, list(dummy = TRUE))
  
  expect_error(combinePriorParameters(bad_input))
})
