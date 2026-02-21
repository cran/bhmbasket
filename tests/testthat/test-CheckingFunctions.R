# Tests for check.evidence.levels ----------------------------------------------

# ------------------------------------------------------------------
# Test: check.evidence.levels passes with valid numeric input
# Input:
#   - evidence_levels = c(0.8, 0.9) inside [0,1]
#   - cohort_names = c("A","B") with same length
#   - analyses_list providing quantiles (not used for numeric check)
# Behaviour:
#   - Function validates that numeric inputs are in range and lengths match.
# Expectations:
#   - No error or warning (expect_silent).
# Why:
#   - Confirms the functioning for simple numeric evidence levels.
# ------------------------------------------------------------------
test_that("check.evidence.levels passes with valid numeric input", {
  analyses_list <- list(list(analysis_parameters = list(quantiles = c(0.1, 0.2, 0.3))))
  expect_silent(check.evidence.levels(
    evidence_levels = c(0.8, 0.9),
    cohort_names = c("A", "B"),
    analyses_list = analyses_list,
    error_evidence_levels = "Invalid evidence levels"
  ))
})

# ------------------------------------------------------------------
# Test: check.evidence.levels fails when lengths differ
# Input:
#   - evidence_levels length = 1
#   - cohort_names length = 2
# Behaviour:
#   - Function checks that each cohort has a corresponding evidence level.
# Expectations:
#   - Error containing "must have the same length".
# Why:
#   - Ensures consistent one-to-one mapping between cohorts and evidence.
# ------------------------------------------------------------------
test_that("check.evidence.levels fails when lengths differ", {
  analyses_list <- list(list(analysis_parameters = list(quantiles = c(0.1, 0.2))))
  expect_error(check.evidence.levels(
    evidence_levels = c(0.8),
    cohort_names = c("A", "B"),
    analyses_list = analyses_list,
    error_evidence_levels = "Invalid evidence levels"
  ), "must have the same length")
})

# ------------------------------------------------------------------
# Test: check.evidence.levels fails when character input has invalid strings
# Input:
#   - evidence_levels = c("mean", "wrong")
# Behaviour:
#   - Character inputs must match allowed keywords (e.g. "mean" only).
# Expectations:
#   - Error "The only string allowed".
# Why:
#   - Guards against unsupported descriptive evidence level strings.
# ------------------------------------------------------------------
test_that("check.evidence.levels fails when character input has invalid strings", {
  analyses_list <- list(list(analysis_parameters = list(quantiles = c(0.1, 0.2))))
  expect_error(check.evidence.levels(
    evidence_levels = c("mean", "wrong"),
    cohort_names = c("A", "B"),
    analyses_list = analyses_list,
    error_evidence_levels = "Invalid evidence levels"
  ), "The only string allowed")
})

# ------------------------------------------------------------------
# Test: check.evidence.levels fails when numeric values out of [0,1]
# Input:
#   - evidence_levels = c(-0.1, 1.2)
# Behaviour:
#   - Numeric evidence levels must be valid probabilities.
# Expectations:
#   - Error containing "Invalid evidence levels".
# Why:
#   - Rejects out-of-bound probability inputs.
# ------------------------------------------------------------------
test_that("check.evidence.levels fails when numeric values out of [0,1]", {
  analyses_list <- list(list(analysis_parameters = list(quantiles = c(0.1, 0.2))))
  expect_error(check.evidence.levels(
    evidence_levels = c(-0.1, 1.2),
    cohort_names = c("A", "B"),
    analyses_list = analyses_list,
    error_evidence_levels = "Invalid evidence levels"
  ), "Invalid evidence levels")
})

# ------------------------------------------------------------------
# Test: check.evidence.levels fails when quantiles do not match
# Input:
#   - evidence_levels interpreted as quantiles (0.5, 0.6)
#   - analysis quantiles = c(0.1, 0.2)
# Behaviour:
#   - Evidence levels used as quantiles must match analysis quantiles exactly.
# Expectations:
#   - Error "must have matches".
# Why:
#   - Prevents mismatched quantile specifications for analysis.
# ------------------------------------------------------------------
test_that("check.evidence.levels fails when quantiles do not match", {
  analyses_list <- list(list(analysis_parameters = list(quantiles = c(0.1, 0.2))))
  expect_error(check.evidence.levels(
    evidence_levels = c(0.5, 0.6),
    cohort_names = c("A", "B"),
    analyses_list = analyses_list,
    error_evidence_levels = "Invalid evidence levels"
  ), "must have matches")
})

# Tests for checkForParallelBackend --------------------------------------------

# ------------------------------------------------------------------
# Test: checkForParallelBackend is silent when backend registered
# Input:
#   - doFuture backend registered + sequential plan set
# Behaviour:
#   - Function checks availability of a foreach-compatible future backend.
# Expectations:
#   - No error (expect_silent).
# Why:
#   - Confirms that with a valid backend registered the function accepts it.
# ------------------------------------------------------------------
test_that("checkForParallelBackend is silent when backend registered", {
  # Register a dummy backend
  doFuture::registerDoFuture()
  future::plan(future::sequential)
  expect_silent(checkForParallelBackend())
})

