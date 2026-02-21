# Tests for getEstimates -------------------------------------------------------

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30)),
  response_rates_list = list(c(0.1, 0.2, 3)),  # 3 triggers historic-rate branch
  n_trials            = 10
)

analyses_list <- performAnalyses(
  scenario_list       = scenarios_list,
  target_rates        = c(0.5, 0.5, 0.5),
  calc_differences    = matrix(c(3, 2, 2, 1), ncol = 2),  # ensures diff cohorts exist
  n_mcmc_iterations   = 100,
  method_names        = "berry"
)

outcome <- createTrial(
  n_subjects   = c(10, 20, 30),
  n_responders = c( 1,  2,  3)
)

outcome_analysis <- performAnalyses(
  scenario_list     = outcome,
  target_rates      = c(0.5, 0.5, 0.5),
  n_mcmc_iterations = 100
)


# -------------------------------------------------------------------
# Test: getEstimates basic structure and contents for simulated scenarios
# Input:
#   - analyses_list: analysis_list with 10 simulated trials, 3 cohorts, diff-cohorts.
# Behaviour:
#   - getEstimates() should return list-per-method of matrices with posterior summaries,
#     Bias and MSE for response-rate parameters and diff parameters.
# Expectations:
#   - Output is a non-empty list.
#   - First inner element is a matrix with columns Mean, SD, 2.5%, 50%, 97.5%, Bias, MSE.
#   - Row names contain p_* for cohorts and "diff" rows for differences.
# Why:
#   - This is the main “happy path”; other tests rely on this shape.
# -------------------------------------------------------------------
test_that("getEstimates works for simulated scenarios and returns sensible structure", {
  res <- getEstimates(analyses_list)
  
  expect_type(res, "list")
  expect_true(length(res) > 0)
  
  first_obj <- res[[1]]
  if (is.list(first_obj)) {
    first_mat <- first_obj[[1]]
  } else {
    first_mat <- first_obj
  }
  
  expect_true(is.matrix(first_mat))
  expect_true(all(c("Mean", "SD", "2.5%", "50%", "97.5%", "Bias", "MSE") %in% colnames(first_mat)))
  expect_true(any(grepl("^p_",   rownames(first_mat))))
  expect_true(any(grepl("diff", rownames(first_mat))))
})


# -------------------------------------------------------------------
# Test: getEstimates with add_parameters and historic-rate handling
# Input:
#   - Same analyses_list as above.
#   - add_parameters = c("mu", "tau", "w_1", "w_2", "w_3").
# Behaviour:
#   - Additional model parameters are appended as extra rows.
#   - For historic cohorts (rr <= 0 or >= 1) true_rr is derived as responders/n_subjects.
#   - Bias and MSE are computed only for p_* rows, not for mu/tau/w_*.
# Expectations:
#   - p_3 row exists and its Bias + true_rr_used equals the median ("50%") estimate.
#   - Column structure identical with and without add_parameters.
#   - Extra rows correspond to mu/tau/w_* and have NA for Bias/MSE.
# Why:
#   - Ensures correct handling of historic cohorts and extra parameters.
# -------------------------------------------------------------------
test_that("additional parameters are added and have NA bias/MSE", {
  res_base <- getEstimates(analyses_list)
  base_obj <- res_base[[1]]
  base_mat <- if (is.list(base_obj)) base_obj[[1]] else base_obj
  
  res_add <- getEstimates(
    analyses_list   = analyses_list,
    add_parameters  = c("mu", "tau", "w_1", "w_2", "w_3")
  )
  add_obj <- res_add[[1]]
  add_mat <- if (is.list(add_obj)) add_obj[[1]] else add_obj
  
  p3_row <- grep("^p_3$", rownames(add_mat))
  expect_length(p3_row, 1)
  
  bias_p3   <- add_mat[p3_row, "Bias"]
  median_p3 <- add_mat[p3_row, "50%"]
  
  true_rr_raw       <- analyses_list[[1]]$scenario_data$response_rates
  n_subj            <- analyses_list[[1]]$scenario_data$n_subjects[1, ]
  expected_true_rr3 <- true_rr_raw[3] / n_subj[[3]]
  
  point_estimate3 <- bias_p3 + expected_true_rr3
  expect_equal(point_estimate3, median_p3)
  
  expect_identical(colnames(base_mat), colnames(add_mat))
  
  extra_rows <- setdiff(rownames(add_mat), rownames(base_mat))
  expect_true(length(extra_rows) > 0)
  expect_true(all(grepl("mu|tau|w_", extra_rows)))
  
  expect_true(all(is.na(add_mat[extra_rows, "Bias"])))
  expect_true(all(is.na(add_mat[extra_rows, "MSE"])))
})


# -------------------------------------------------------------------
# Test: getEstimates for single-trial outcome
# Input:
#   - outcome_analysis: analysis_list from createTrial() (one trial only).
# Behaviour:
#   - With n_trials = 1, getEstimates() should *not* compute Bias/MSE,
#     only posterior summaries.
# Expectations:
#   - Returned inner matrix has columns exactly: Mean, SD, 2.5%, 50%, 97.5%.
# Why:
#   - Distinguishes between simulation-based metrics and outcome-only analyses.
# -------------------------------------------------------------------
test_that("single-trial outcome returns only posterior summaries (no bias/MSE)", {
  res_single <- getEstimates(outcome_analysis)
  single_obj <- res_single[[1]]
  single_mat <- if (is.list(single_obj)) single_obj[[1]] else single_obj
  
  expect_true(is.matrix(single_mat))
  expect_identical(
    colnames(single_mat),
    c("Mean", "SD", "2.5%", "50%", "97.5%")
  )
})


# -------------------------------------------------------------------
# Test: validation of alpha_level and add_parameters in getEstimates
# Input:
#   - analyses_list from the main fixture.
# Behaviour:
#   - alpha_level must be in (0, 1) and must correspond to stored quantiles.
#   - add_parameters must exist in at least one method’s posterior draws.
# Expectations:
#   - alpha_level = 0.07 (not among stored quantiles) -> error.
#   - add_parameters = "totally_unknown_param" -> error with appropriate message.
# Why:
#   - Ensures defensive checks for user-supplied quantiles and parameter names.
# -------------------------------------------------------------------
test_that("alpha_level and add_parameters are validated correctly", {
  expect_error(
    getEstimates(analyses_list, alpha_level = 0.07)
  )
  
  expect_error(
    getEstimates(analyses_list, add_parameters = c("totally_unknown_param"))
  )
})


# -------------------------------------------------------------------
# Test: point_estimator behaviour and basic type checks in getEstimates
# Input:
#   - analyses_list, point_estimator = "median" and "mean".
# Behaviour:
#   - Both estimators should work and yield same matrix shape.
#   - An invalid point_estimator should raise an error.
#   - Non-analysis_list input should raise an error.
# Expectations:
#   - dim() and colnames() identical for median vs mean.
#   - point_estimator = "mode" -> error.
#   - analyses_list = list(a = 1) -> error on class.
# Why:
#   - Confirms configurability of the estimator and robust input validation.
# -------------------------------------------------------------------
test_that("point_estimator argument is respected and input type is validated", {
  res_median <- getEstimates(analyses_list, point_estimator = "median")
  res_mean   <- getEstimates(analyses_list, point_estimator = "mean")
  
  med_obj  <- res_median[[1]]
  mean_obj <- res_mean[[1]]
  
  med_mat  <- if (is.list(med_obj))  med_obj[[1]]  else med_obj
  mean_mat <- if (is.list(mean_obj)) mean_obj[[1]] else mean_obj
  
  expect_identical(dim(med_mat),   dim(mean_mat))
  expect_identical(colnames(med_mat), colnames(mean_mat))
  
  expect_error(
    getEstimates(analyses_list, point_estimator = "mode")
  )
  
  expect_error(
    getEstimates(list(a = 1))
  )
})


# -------------------------------------------------------------------
# Additional explicit validation tests for getEstimates:
# - alpha_level outside (0,1)
# - alpha_level with unavailable quantiles
# - add_parameters not found
# - basic structure for single-trial
# - duplicate explicit checks for invalid point_estimator and analyses_list
# -------------------------------------------------------------------

test_that("throws error for invalid alpha_level (outside 0,1)", {
  expect_error(getEstimates(analyses_list, alpha_level = 1.5), "alpha_level")
})

test_that("throws error if alpha_level quantiles not available", {
  expect_error(
    getEstimates(analyses_list, alpha_level = 0.123),
    "must be among the stored quantiles"
  )
})

test_that("throws error if add_parameters not found in any method", {
  expect_error(
    getEstimates(analyses_list, add_parameters = c("nonexistent")),
    "do not occur"
  )
})

test_that("works for single trial outcome (basic check)", {
  result <- getEstimates(outcome_analysis)
  expect_type(result, "list")
})

test_that("invalid point_estimator throws error explicitly", {
  expect_error(
    getEstimates(analyses_list, point_estimator = "mode"),
    "point_estimator"
  )
})

test_that("invalid analyses_list class throws error explicitly", {
  expect_error(
    getEstimates(list()),
    "analyses_list"
  )
})


# Tests for getGoDecisions / getGoProbabilities / print.decision_list ----------

set.seed(123)

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30)),
  response_rates_list = list(c(0.3, 0.5, 0.7)),
  n_trials            = 10
)

analyses_list <- performAnalyses(
  scenario_list     = scenarios_list,
  target_rates      = c(0.3, 0.3, 0.3),
  n_mcmc_iterations = 100
)

default_cohorts <- c("p_1", "p_2", "p_3")


# -------------------------------------------------------------------
# Test: getGoDecisions – analyses_list must have correct class
# Input:
#   - analyses_list = list() (wrong type)
# Behaviour:
#   - checkmate::assertClass enforces class 'analysis_list'.
# Expectations:
#   - Error mentioning 'analysis_list' is thrown.
# Why:
#   - Prevents silent misuse with arbitrary lists.
# -------------------------------------------------------------------
test_that("getGoDecisions: errors if analyses_list is not of class 'analysis_list'", {
  expect_error(
    getGoDecisions(
      analyses_list   = list(),
      cohort_names    = "p_1",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE))
    ),
    "analysis_list"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – evidence_levels and boundary_rules required
# Input:
#   - Missing evidence_levels or missing boundary_rules.
# Behaviour:
#   - Custom error messages if either argument is missing.
# Expectations:
#   - Missing evidence_levels -> error message refers to evidence_levels.
#   - Missing boundary_rules  -> error message refers to boundary_rules.
# Why:
#   - These are essential pieces of the decision rule definition.
# -------------------------------------------------------------------
test_that("getGoDecisions: errors when evidence_levels or boundary_rules are missing", {
  expect_error(
    getGoDecisions(
      analyses_list  = analyses_list,
      cohort_names   = "p_1",
      boundary_rules = quote(c(TRUE))
    ),
    "evidence_levels"
  )
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "p_1",
      evidence_levels = 0.5
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – invalid cohort_names
# Input:
#   - cohort_names = "invalid" (not present in posterior columns).
# Behaviour:
#   - assertSubset ensures cohort_names correspond to posterior parameters.
# Expectations:
#   - Error referencing 'cohort_names'.
# Why:
#   - Prevents mis-specification of cohort parameter names.
# -------------------------------------------------------------------
test_that("getGoDecisions: errors if cohort_names are not valid posterior parameters", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "invalid",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE))
    ),
    "cohort_names"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – overall_min_gos must be positive integer
# Input:
#   - overall_min_gos = 0L.
# Behaviour:
#   - assertInt(overall_min_gos, lower = 1) rejects non-positive values.
# Expectations:
#   - Error referencing 'overall_min_gos'.
# Why:
#   - The minimum number of cohort-wise Go decisions must be >= 1.
# -------------------------------------------------------------------
test_that("getGoDecisions: errors if overall_min_gos is not a positive integer", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "p_1",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE)),
      overall_min_gos = 0L
    ),
    "overall_min_gos"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – numeric evidence_levels in (0,1)
# Input:
#   - evidence_levels = c(-0.1, 0.5, 1.3).
# Behaviour:
#   - check.evidence.levels enforces all levels are in (0,1).
# Expectations:
#   - Error is thrown for invalid vector.
# Why:
#   - Evidence thresholds must be valid probabilities.
# -------------------------------------------------------------------
test_that("getGoDecisions: numeric evidence_levels outside (0,1) cause an error", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(-0.1, 0.5, 1.3),
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    )
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – list evidence_levels with valid entries
# Input:
#   - evidence_levels as list of numeric vectors in (0,1).
# Behaviour:
#   - check.evidence.levels is applied to each list element.
# Expectations:
#   - No error is thrown for valid list structure.
# Why:
#   - Supports method-specific evidence levels.
# -------------------------------------------------------------------
test_that("getGoDecisions: list evidence_levels all elements valid", {
  ev_list_ok <- list(
    rep(0.5, length(default_cohorts)),
    rep(0.5, length(default_cohorts))
  )
  
  expect_silent(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = ev_list_ok,
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    )
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – non-list boundary_rules must be a language object
# Input:
#   - boundary_rules = 1 (numeric, not a call).
# Behaviour:
#   - check_boundary_rules rejects non-language boundary rules.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Decision rules must be expressed as quoted R expressions.
# -------------------------------------------------------------------
test_that("getGoDecisions: non-list boundary_rules must be a language object", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = 1
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – non-list boundary_rules must use c(...)
# Input:
#   - boundary_rules = quote(list(TRUE, TRUE, TRUE)).
# Behaviour:
#   - The call head must be `c`, so list(...) is rejected.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - The internal logic expects a c(...) vector of logical decisions.
# -------------------------------------------------------------------
test_that("getGoDecisions: non-list boundary_rules must start with c()", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = quote(list(TRUE, TRUE, TRUE))
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – non-list boundary_rules length must match #cohorts
# Input:
#   - 3 cohorts but boundary_rules = quote(c(TRUE, TRUE)).
# Behaviour:
#   - Length mismatch triggers error in check_boundary_rules.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Each cohort must have one logical decision rule.
# -------------------------------------------------------------------
test_that("getGoDecisions: non-list boundary_rules must have one entry per cohort", {
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = quote(c(TRUE, TRUE))
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – list boundary_rules elements must be language
# Input:
#   - boundary_rules = list(123).
# Behaviour:
#   - Each element of the list must be a call; a numeric element is invalid.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Supports multiple methods with method-specific decision rules.
# -------------------------------------------------------------------
test_that("getGoDecisions: list boundary_rules each element must be a language object", {
  br <- list(123)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – list boundary_rules head must be c()
# Input:
#   - boundary_rules = list(quote(list(TRUE, TRUE, TRUE))).
# Behaviour:
#   - Each element must start with c(...); list(...) is rejected.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Same reasoning as the non-list case, but in list form.
# -------------------------------------------------------------------
test_that("getGoDecisions: list boundary_rules each element must start with c()", {
  br <- list(quote(list(TRUE, TRUE, TRUE)))
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – list boundary_rules length must match #cohorts
# Input:
#   - boundary_rules = list(quote(c(TRUE, TRUE))) for 3 cohorts.
# Behaviour:
#   - Length mismatch per element triggers error.
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Ensures consistent vector length across all cohorts.
# -------------------------------------------------------------------
test_that("getGoDecisions: list boundary_rules must match number of cohorts", {
  br <- list(quote(c(TRUE, TRUE)))
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – boundary_rules list length vs method_names
# Input:
#   - boundary_rules list longer than number of methods.
# Behaviour:
#   - Length check should reject boundary_rules with length > length(method_names).
# Expectations:
#   - Error referencing 'boundary_rules'.
# Why:
#   - Each method can have at most one decision rule expression.
# -------------------------------------------------------------------
test_that("errors if boundary_rules list is longer than method_names", {
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  br   <- rep(list(quote(c(TRUE, TRUE, TRUE))), length(m) + 1L)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = coh,
      evidence_levels = ev,
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – evidence_levels list length vs method_names
# Input:
#   - evidence_levels list longer than number of methods.
# Behaviour:
#   - Length check should reject gamma_levels with length > length(method_names).
# Expectations:
#   - Error referencing 'evidence_levels'.
# Why:
#   - Each method can have at most one vector of evidence thresholds.
# -------------------------------------------------------------------
test_that("errors if evidence_levels list is longer than method_names", {
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  
  ev_vec <- c(0.5, 0.5, 0.5)
  ev     <- rep(list(ev_vec), length(m) + 1L)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = coh,
      evidence_levels = ev,
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    ),
    "evidence_levels"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – single boundary_rules expression recycled
# Input:
#   - boundary_rules = quote(c(TRUE, TRUE, TRUE)) and multiple methods.
# Behaviour:
#   - Single expression is replicated for each method internally.
# Expectations:
#   - decision_rules$boundary_rules has length = #methods,
#     each element identical to the original rule.
# Why:
#   - Confirms convenience recycling logic.
# -------------------------------------------------------------------
test_that("single boundary_rules expression is recycled to all methods", {
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  rule <- quote(c(TRUE, TRUE, TRUE))
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = coh,
    evidence_levels = ev,
    boundary_rules  = rule
  )
  
  br <- dec[[1]]$decision_rules$boundary_rules
  
  expect_identical(length(br), length(m))
  for (i in seq_along(br)) {
    expect_true(identical(br[[i]], rule))
  }
})


# -------------------------------------------------------------------
# Test: getGoDecisions – single evidence_levels vector recycled
# Input:
#   - evidence_levels = c(0.5, 0.5, 0.5) and multiple methods.
# Behaviour:
#   - Single numeric vector is replicated per method.
# Expectations:
#   - decision_rules$gamma_levels has length = #methods,
#     each element equal to the original vector.
# Why:
#   - Confirms convenience recycling for evidence thresholds.
# -------------------------------------------------------------------
test_that("single evidence_levels vector is recycled to all methods", {
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = coh,
    evidence_levels = ev,
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  gamma <- dec[[1]]$decision_rules$gamma_levels
  
  expect_identical(length(gamma), length(m))
  for (i in seq_along(gamma)) {
    expect_equal(gamma[[i]], ev)
  }
})


# -------------------------------------------------------------------
# Test: getGoDecisions – method names consistent across scenarios
# Input:
#   - analyses_list with identical method_names per scenario.
# Behaviour:
#   - Method name matrix check passes; decisions are computed.
# Expectations:
#   - Returned object is of class decision_list.
# Why:
#   - Sanity check for the main use-case of multiple scenarios with same methods.
# -------------------------------------------------------------------
test_that("getGoDecisions: succeeds when all scenarios use identical method_names", {
  res <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  expect_s3_class(res, "decision_list")
})


# -------------------------------------------------------------------
# Test: getGoDecisions – error when scenarios use different methods
# Input:
#   - analyses_list where scenario_2 has method_names reversed vs scenario_1.
# Behaviour:
#   - Method-name consistency check fails.
# Expectations:
#   - Error with message “analysed with different methods”.
# Why:
#   - Mixed-method setups per scenario are not supported.
# -------------------------------------------------------------------
test_that("getGoDecisions: errors when scenarios were analysed with different methods", {
  bad <- analyses_list
  bad$scenario_2 <- bad$scenario_1
  bad$scenario_2$analysis_parameters$method_names <-
    rev(bad$scenario_2$analysis_parameters$method_names)
  
  expect_error(
    getGoDecisions(
      analyses_list   = bad,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    ),
    "analysed with different methods"
  )
})


# -------------------------------------------------------------------
# Test: getGoDecisions – structure of decision_list and 'overall' column
# Input:
#   - Valid analyses_list, three cohorts, evidence levels.
# Behaviour:
#   - getGoDecisions returns decision_list with:
#     - decisions_list per method/scenario including 'overall' and cohort columns.
#     - decision_rules storing cohort_names and gamma_levels.
# Expectations:
#   - 'overall' column exists.
#   - At least one cohort-level decision column exists.
#   - decision_rules$cohort_names identical to input.
#   - stored gamma levels contain the supplied values.
# Why:
#   - Verifies structural contract of decision_list objects.
# -------------------------------------------------------------------
test_that("getGoDecisions: returns decision_list with overall and cohort decisions", {
  decisions <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  scen1  <- decisions[[1]]
  m_dec  <- as.matrix(scen1$decisions_list[[1]])
  
  expect_true("overall" %in% colnames(m_dec))
  cohort_cols <- setdiff(colnames(m_dec), "overall")
  expect_true(length(cohort_cols) >= 1)
  
  expect_identical(scen1$decision_rules$cohort_names, default_cohorts)
  
  stored_gamma_flat <- unlist(scen1$decision_rules$gamma_levels, use.names = FALSE)
  expect_true(all(c(0.5, 0.5, 0.8) %in% stored_gamma_flat))
})


# -------------------------------------------------------------------
# Test: getGoDecisions – semantics of overall_min_gos = 1
# Input:
#   - overall_min_gos = 1, boundary_rules always TRUE.
# Behaviour:
#   - overall decision is TRUE if at least one cohort is TRUE.
# Expectations:
#   - 'overall' column equals row-wise indicator of ≥1 TRUE among cohorts.
# Why:
#   - Checks threshold behaviour for overall Go decision.
# -------------------------------------------------------------------
test_that("overall_min_gos = 1 means overall Go if at least one cohort is Go", {
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = rep(0.5, length(default_cohorts)),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE)),
    overall_min_gos = 1L
  )
  
  m <- as.matrix(dec[[1]]$decisions_list[[1]])
  coh_cols <- setdiff(colnames(m), "overall")
  expect_true(length(coh_cols) >= 1)
  
  coh <- m[, coh_cols, drop = FALSE] > 0
  overall_calc <- apply(coh, 1, function(x) sum(x) >= 1L)
  
  overall <- as.logical(m[, "overall"])
  expect_identical(overall, overall_calc)
})


# -------------------------------------------------------------------
# Test: getGoDecisions – semantics of overall_min_gos = 2
# Input:
#   - overall_min_gos = 2, boundary_rules always TRUE.
# Behaviour:
#   - overall decision is TRUE if at least two cohorts are TRUE.
# Expectations:
#   - 'overall' column equals row-wise indicator of ≥2 TRUE among cohorts.
# Why:
#   - Confirms threshold is correctly applied for higher overall minimum.
# -------------------------------------------------------------------
test_that("overall_min_gos = 2 means overall Go if at least two cohorts are Go", {
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = rep(0.5, length(default_cohorts)),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE)),
    overall_min_gos = 2L
  )
  
  m <- as.matrix(dec[[1]]$decisions_list[[1]])
  coh_cols <- setdiff(colnames(m), "overall")
  expect_true(length(coh_cols) >= 1)
  
  coh <- m[, coh_cols, drop = FALSE] > 0
  overall_calc <- apply(coh, 1, function(x) sum(x) >= 2L)
  
  overall <- as.logical(m[, "overall"])
  expect_identical(overall, overall_calc)
})


# -------------------------------------------------------------------
# Test: getGoDecisions – no resurrection of previously stopped cohorts
# Input:
#   - analyses_list with scenario_data$previous_analyses$go_decisions present.
# Behaviour:
#   - New go_decisions are multiplied by previous_gos, so a previously stopped
#     cohort cannot become Go again.
# Expectations:
#   - No cell where new_go == TRUE and previous_go == FALSE.
# Why:
#   - Implements the intended “no resurrection” rule across analyses.
# -------------------------------------------------------------------
test_that("getGoDecisions: previous go_decisions prevent resurrection of stopped cohorts", {
  if (is.null(analyses_list[[1]]$scenario_data$previous_analyses)) {
    skip("previous_analyses not available in scenario_data")
  }
  
  decisions <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  new_mat     <- as.matrix(decisions[[1]]$decisions_list[[1]])
  new_cohcols <- setdiff(colnames(new_mat), "overall")
  new_gos     <- new_mat[, new_cohcols, drop = FALSE] > 0
  
  prev_mat <- analyses_list[[1]]$scenario_data$previous_analyses$go_decisions
  prev_gos <- as.matrix(prev_mat[, -1, drop = FALSE]) > 0
  
  expect_identical(dim(prev_gos), dim(new_gos))
  
  resurrected <- new_gos & !prev_gos
  expect_false(any(resurrected))
})


# Tests for getGoProbabilities -------------------------------------------------

set.seed(456)

scenarios_list_prob <- simulateScenarios(
  n_subjects_list     = list(c(10, 20)),
  response_rates_list = list(rep(0.9, 2)),
  n_trials            = 10
)

analyses_list_prob <- performAnalyses(
  scenario_list     = scenarios_list_prob,
  target_rates      = rep(0.5, 2),
  n_mcmc_iterations = 100
)

prob_cohorts <- c("p_1", "p_2")

go_decisions_list <- getGoDecisions(
  analyses_list   = analyses_list_prob,
  cohort_names    = prob_cohorts,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(x[1] > 0.8, x[2] > 0.6))
)

nogo_decisions_list <- getGoDecisions(
  analyses_list   = analyses_list_prob,
  cohort_names    = prob_cohorts,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(x[1] < 0.5, x[2] < 0.3))
)


# -------------------------------------------------------------------
# Test: getGoProbabilities – go_decisions_list must be decision_list
# Input:
#   - go_decisions_list = list().
# Behaviour:
#   - assertClass(go_decisions_list, "decision_list") should fail.
# Expectations:
#   - Error mentioning 'go_decisions_list'.
# Why:
#   - Ensures user passes the correct object type.
# -------------------------------------------------------------------
test_that("getGoProbabilities: errors if go_decisions_list is not a decision_list", {
  expect_error(
    getGoProbabilities(go_decisions_list = list()),
    "go_decisions_list"
  )
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – nogo_decisions_list must be decision_list or NULL
# Input:
#   - nogo_decisions_list = list().
# Behaviour:
#   - assertClass(..., null.ok = TRUE) fails when object is not decision_list.
# Expectations:
#   - Error mentioning 'decision_list'.
# Why:
#   - Prevents malformed NoGo objects.
# -------------------------------------------------------------------
test_that("getGoProbabilities: errors if nogo_decisions_list has wrong class", {
  expect_error(
    getGoProbabilities(
      go_decisions_list   = go_decisions_list,
      nogo_decisions_list = list()
    ),
    "decision_list"
  )
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – Go/NoGo matrices must have same dimensions
# Input:
#   - bad_nogo with first method’s matrix missing a column.
# Behaviour:
#   - Dimension check should fail if the shapes differ.
# Expectations:
#   - Error is thrown (no specific regex needed).
# Why:
#   - Go and NoGo must align element-wise.
# -------------------------------------------------------------------
test_that("getGoProbabilities: errors if Go and NoGo matrices have different dimensions", {
  bad_nogo <- nogo_decisions_list
  bad_nogo[[1]]$decisions_list[[1]] <-
    bad_nogo[[1]]$decisions_list[[1]][, -1, drop = FALSE]
  
  expect_error(
    getGoProbabilities(
      go_decisions_list   = go_decisions_list,
      nogo_decisions_list = bad_nogo
    )
  )
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – Go-only case structure
# Input:
#   - go_decisions_list, nogo_decisions_list = NULL.
# Behaviour:
#   - Returns list-per-method of matrices with single row "Go".
# Expectations:
#   - probs is a list of lists of matrices.
#   - First matrix has rownames "Go" and at least one column.
# Why:
#   - Checks basic output shape in the simpler Go-only mode.
# -------------------------------------------------------------------
test_that("getGoProbabilities: Go-only call returns list-of-lists of matrices", {
  probs <- getGoProbabilities(go_decisions_list)
  
  expect_type(probs, "list")
  expect_true(all(sapply(probs, is.list)))
  
  first_method   <- names(probs)[1]
  first_scenario <- names(probs[[first_method]])[1]
  mat            <- probs[[first_method]][[first_scenario]]
  
  expect_true(is.matrix(mat))
  expect_identical(rownames(mat), "Go")
  expect_true(ncol(mat) >= 1)
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – Go row equals column means of go_decisions
# Input:
#   - go_decisions_list with full decision matrices.
# Behaviour:
#   - Row "Go" is defined as colMeans(go_decisions).
# Expectations:
#   - For each method/scenario, mat_prob["Go", ] == colMeans(go_mat).
# Why:
#   - Validates the probability computation logic.
# -------------------------------------------------------------------
test_that("getGoProbabilities: Go row equals column means of go_decisions", {
  probs <- getGoProbabilities(go_decisions_list)
  
  method_names   <- names(go_decisions_list[[1]]$decisions_list)
  scenario_names <- names(go_decisions_list)
  
  for (m in method_names) {
    for (s in scenario_names) {
      mat_prob <- probs[[m]][[s]]
      go_mat   <- go_decisions_list[[s]]$decisions_list[[m]]
      
      expected <- colMeans(go_mat)
      expect_equal(
        mat_prob["Go", ],
        expected,
        info = paste("Mismatch in method", m, "scenario", s)
      )
    }
  }
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – Go and NoGo rows match colMeans of inputs
# Input:
#   - go_decisions_list and nogo_decisions_list for same scenarios.
# Behaviour:
#   - "Go" row = colMeans(go_mat), "NoGo" row = colMeans(nogo_mat).
# Expectations:
#   - Equality holds for all methods and scenarios.
# Why:
#   - Confirms correct aggregation of both Go and NoGo decisions.
# -------------------------------------------------------------------
test_that("getGoProbabilities: Go and NoGo rows match colMeans of input decisions", {
  probs <- getGoProbabilities(
    go_decisions_list   = go_decisions_list,
    nogo_decisions_list = nogo_decisions_list
  )
  
  method_names   <- names(go_decisions_list[[1]]$decisions_list)
  scenario_names <- names(go_decisions_list)
  
  for (m in method_names) {
    for (s in scenario_names) {
      mat_prob <- probs[[m]][[s]]
      
      go_mat   <- go_decisions_list[[s]]$decisions_list[[m]]
      nogo_mat <- nogo_decisions_list[[s]]$decisions_list[[m]]
      
      expected_go   <- colMeans(go_mat)
      expected_nogo <- colMeans(nogo_mat)
      
      expect_equal(
        mat_prob["Go", ],
        expected_go,
        info = paste("Go row mismatch in method", m, "scenario", s)
      )
      
      expect_equal(
        mat_prob["NoGo", ],
        expected_nogo,
        info = paste("NoGo row mismatch in method", m, "scenario", s)
      )
    }
  }
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – column sums to 1 with Go/NoGo/Consider
# Input:
#   - go_decisions_list and nogo_decisions_list.
# Behaviour:
#   - "Consider" row is computed as 1 - p(Go) - p(NoGo), so column sums = 1.
# Expectations:
#   - For each scenario, colSums(matrix) == 1 exactly.
# Why:
#   - Sanity check that probabilities are normalized.
# -------------------------------------------------------------------
test_that("getGoProbabilities: columns sum to 1 when NoGo decisions are provided", {
  probs <- getGoProbabilities(
    go_decisions_list   = go_decisions_list,
    nogo_decisions_list = nogo_decisions_list
  )
  
  first_method <- names(probs)[1]
  for (scen_name in names(probs[[first_method]])) {
    mat <- probs[[first_method]][[scen_name]]
    
    col_sums <- colSums(mat)
    expect_true(
      all(abs(col_sums - 1) == 0),
      info = paste("Column sums not 1 for scenario", scen_name)
    )
  }
})


# -------------------------------------------------------------------
# Test: getGoProbabilities – overlapping Go and NoGo decisions
# Input:
#   - Modified go_decisions_list and nogo_decisions_list where at least one
#     cell is TRUE in both.
# Behaviour:
#   - Overlap check should throw an error.
# Expectations:
#   - Error containing “both go and nogo decisions”.
# Why:
#   - A trial cannot be simultaneously Go and NoGo in the same cell.
# -------------------------------------------------------------------
test_that("getGoProbabilities: errors if any decision is both Go and NoGo", {
  overlap_go   <- go_decisions_list
  overlap_nogo <- nogo_decisions_list
  
  scen1_name <- names(overlap_go)[1]
  meth1_name <- names(overlap_go[[scen1_name]]$decisions_list)[1]
  
  overlap_go[[scen1_name]]$decisions_list[[meth1_name]][1, 1]   <- TRUE
  overlap_nogo[[scen1_name]]$decisions_list[[meth1_name]][1, 1] <- TRUE
  
  expect_error(
    getGoProbabilities(
      go_decisions_list   = overlap_go,
      nogo_decisions_list = overlap_nogo
    ),
    "both go and nogo decisions"
  )
})


# Tests for print.decision_list ------------------------------------------------

set.seed(789)

scenarios_print <- simulateScenarios(
  n_subjects_list     = list(c(10, 20), c(10, 20)),
  response_rates_list = list(c(0.3, 0.5), c(0.4, 0.6)),
  n_trials            = 5
)

analyses_print <- performAnalyses(
  scenario_list     = scenarios_print,
  target_rates      = c(0.3, 0.3),
  n_mcmc_iterations = 80
)

coh_print <- c("p_1", "p_2")

decisions_dl <- getGoDecisions(
  analyses_list   = analyses_print,
  cohort_names    = coh_print,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(TRUE, TRUE))
)

n_scenarios  <- length(decisions_dl)
n_methods    <- length(decisions_dl[[1]]$decisions_list)


# -------------------------------------------------------------------
# Test: print.decision_list – header shows scenario/method counts
# Input:
#   - decisions_dl with known number of scenarios and methods.
# Behaviour:
#   - Header line prints “decision_list of N scenario(s) with M method(s)”.
# Expectations:
#   - Output contains the constructed header pattern.
# Why:
#   - Confirms informative summary in the print method.
# -------------------------------------------------------------------
test_that("print.decision_list: header shows correct number of scenarios and methods", {
  header_pattern <- sprintf(
    "decision_list of %d scenario%s with %d method%s",
    n_scenarios,
    ifelse(n_scenarios == 1, "", "s"),
    n_methods,
    ifelse(n_methods == 1, "", "s")
  )
  
  expect_output(
    print(decisions_dl),
    header_pattern
  )
})


# -------------------------------------------------------------------
# Test: print.decision_list – scenario sections appear
# Input:
#   - decisions_dl with scenario names.
# Behaviour:
#   - For each scenario, the print method outputs “  - <scenario_name>”.
# Expectations:
#   - Each scenario name appears at least once in captured output.
# Why:
#   - Confirms per-scenario blocks are printed.
# -------------------------------------------------------------------
test_that("print.decision_list: prints a section for each scenario", {
  scen_names <- names(decisions_dl)
  expect_true(length(scen_names) >= 1)
  
  output <- capture.output(print(decisions_dl))
  
  for (nm in scen_names) {
    expect_true(
      any(grepl(paste0("\\b", nm, "\\b"), output)),
      info = paste("Scenario name", nm, "not found in printed output")
    )
  }
})


# -------------------------------------------------------------------
# Test: print.decision_list – method labels appear
# Input:
#   - decisions_dl with method_names present.
# Behaviour:
#   - Row labels are built from method_names via firstUpper().
# Expectations:
#   - Each method name (firstUpper) appears somewhere in output.
# Why:
#   - Verifies mapping from internal method_names to printed labels.
# -------------------------------------------------------------------
test_that("print.decision_list: each method name appears in printed matrix rows", {
  method_names <- names(decisions_dl[[1]]$decisions_list)
  expect_true(length(method_names) >= 1)
  
  output <- capture.output(print(decisions_dl))
  
  for (m in method_names) {
    m_upper <- paste0(toupper(substr(m, 1, 1)), substr(m, 2, nchar(m)))
    expect_true(
      any(grepl(m_upper, output, fixed = TRUE)),
      info = paste("Method label", m_upper, "not found in printed output")
    )
  }
})


# -------------------------------------------------------------------
# Test: print.decision_list – digits argument works
# Input:
#   - digits = 1 and digits = 4.
# Behaviour:
#   - Numeric entries are rounded with the given digits; no error.
# Expectations:
#   - Both calls produce output containing “decision_list of”.
# Why:
#   - Ensures digits parameter is respected but does not break printing.
# -------------------------------------------------------------------
test_that("print.decision_list: digits argument is accepted and does not error", {
  expect_output(
    print(decisions_dl, digits = 1),
    "decision_list of"
  )
  
  expect_output(
    print(decisions_dl, digits = 4),
    "decision_list of"
  )
})


# -------------------------------------------------------------------
# Test: print.decision_list – handles non-NULL decision_rules
# Input:
#   - decisions_dl with decision_rules present.
# Behaviour:
#   - Internal block rewrites boundary_rules for pretty printing.
# Expectations:
#   - Printing completes without error and includes header.
# Why:
#   - Confirms the complex decision_rules formatting branch works.
# -------------------------------------------------------------------
test_that("print.decision_list: handles non-NULL decision_rules without error", {
  expect_output(
    print(decisions_dl),
    "decision_list of"
  )
})


# -------------------------------------------------------------------
# Test: print.decision_list – works when decision_rules are NULL
# Input:
#   - decisions_dl with decision_rules removed.
# Behaviour:
#   - Skips decision_rules formatting and prints summary based only on decisions.
# Expectations:
#   - Header is printed as usual, no error occurs.
# Why:
#   - Ensures robust printing even if decision_rules slot is absent.
# -------------------------------------------------------------------
test_that("print.decision_list: works when decision_rules are NULL", {
  dl_no_rules <- decisions_dl
  for (i in seq_along(dl_no_rules)) {
    dl_no_rules[[i]]$decision_rules <- NULL
  }
  
  header_pattern <- sprintf(
    "decision_list of %d scenario%s with %d method%s",
    n_scenarios,
    ifelse(n_scenarios == 1, "", "s"),
    n_methods,
    ifelse(n_methods == 1, "", "s")
  )
  
  expect_output(
    print(dl_no_rules),
    header_pattern
  )
})


# -------------------------------------------------------------------
# Test: print.decision_list – one row per method per scenario
# Input:
#   - decisions_dl with >1 scenario and potentially multiple methods.
# Behaviour:
#   - For each scenario, mat_out = rbind over methods; rownames become
#     “    - <method>”.
# Expectations:
#   - Total number of lines starting with “    - ” equals n_scenarios * #methods.
# Why:
#   - Confirms the aggregation over methods matches the printed labels.
# -------------------------------------------------------------------
test_that("print.decision_list: for multiple methods, one row per method is printed per scenario", {
  go_probs <- getGoProbabilities(decisions_dl)
  
  scenario_index <- 1L
  mat_ref <- do.call(
    rbind,
    lapply(go_probs, function(y) y[[scenario_index]])
  )
  
  output <- capture.output(print(decisions_dl))
  
  method_label_lines <- grep("^    - ", output, value = TRUE)
  expected_total_rows <- n_scenarios * nrow(mat_ref)
  
  expect_equal(
    length(method_label_lines),
    expected_total_rows,
    info = paste(
      "Number of printed method label lines (", length(method_label_lines),
      ") does not match expected total rows (", expected_total_rows, ")"
    )
  )
})


# ===================================================================
# Tests for getAverageNSubjects -------------------------------------
# ===================================================================

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30), c(15, 25, 35)),
  response_rates_list = list(c(0.1, 0.2, 0.3), c(0.15, 0.25, 0.35)),
  n_trials            = 2
)


# -------------------------------------------------------------------
# Test: getAverageNSubjects – structure for multiple scenarios
# Input:
#   - scenario_list with two scenarios.
# Behaviour:
#   - Returns list of numeric vectors, one per scenario.
# Expectations:
#   - Result is a list with length equal to number of scenarios.
#   - Each element is numeric.
# Why:
#   - Basic contract of getAverageNSubjects.
# -------------------------------------------------------------------
test_that("returns correct structure for multiple scenarios", {
  result <- getAverageNSubjects(scenarios_list)
  expect_type(result, "list")
  expect_equal(length(result), length(scenarios_list))
  expect_true(all(sapply(result, is.numeric)))
})


# -------------------------------------------------------------------
# Test: getAverageNSubjects – computes column means of n_subjects
# Input:
#   - scenarios_list with known n_subjects for scenario_1.
# Behaviour:
#   - For each scenario, colMeans() of n_subjects is returned.
# Expectations:
#   - result[[1]] equals colMeans(scenarios_list[[1]]$n_subjects).
# Why:
#   - Verifies the actual calculation performed.
# -------------------------------------------------------------------
test_that("computes correct column means", {
  expected <- colMeans(scenarios_list[[1]]$n_subjects)
  result <- getAverageNSubjects(scenarios_list)
  expect_equal(result[[1]], expected)
})


# -------------------------------------------------------------------
# Test: getAverageNSubjects – works for single scenario
# Input:
#   - single_scenario scenario_list with one element.
# Behaviour:
#   - Returns list of length 1, named as input, with column means.
# Expectations:
#   - names(result) match names(single_scenario).
#   - result[[1]] equals colMeans(single_scenario[[1]]$n_subjects).
# Why:
#   - Handles degenerate one-scenario case correctly.
# -------------------------------------------------------------------
test_that("works for single scenario", {
  single_scenario <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(c(0.1, 0.2, 0.3)),
    n_trials            = 3
  )
  result <- getAverageNSubjects(single_scenario)
  expect_equal(names(result), names(single_scenario))
  expect_equal(result[[1]], colMeans(single_scenario[[1]]$n_subjects))
})


# -------------------------------------------------------------------
# Test: getAverageNSubjects – invalid input class
# Input:
#   - list() without class "scenario_list".
# Behaviour:
#   - assertClass(scenario_list, "scenario_list") fails.
# Expectations:
#   - Error mentioning 'scenario_list'.
# Why:
#   - Prevents accidental misuse with arbitrary lists.
# -------------------------------------------------------------------
test_that("throws error for invalid input class", {
  expect_error(getAverageNSubjects(list()), "scenario_list")
})


# -------------------------------------------------------------------
# Test: getAverageNSubjects – empty scenario_list
# Input:
#   - empty list with class "scenario_list".
# Behaviour:
#   - lapply over empty list produces empty list.
# Expectations:
#   - result is list().
# Why:
#   - Defines behaviour for degenerate empty input gracefully.
# -------------------------------------------------------------------
test_that("returns empty list for empty scenario_list", {
  empty_list <- structure(list(), class = "scenario_list")
  result <- getAverageNSubjects(empty_list)
  expect_equal(result, list())
})


# ===================================================================
# Tests for negateGoDecisions ---------------------------------------
# ===================================================================

set.seed(101)

scenarios_neg <- simulateScenarios(
  n_subjects_list     = list(c(10, 15)),
  response_rates_list = list(c(0.4, 0.7)),
  n_trials            = 6
)

analyses_neg <- performAnalyses(
  scenario_list     = scenarios_neg,
  target_rates      = c(0.4, 0.4),
  n_mcmc_iterations = 80
)

coh_neg <- c("p_1", "p_2")

go_decisions_list_neg <- getGoDecisions(
  analyses_list   = analyses_neg,
  cohort_names    = coh_neg,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(TRUE, TRUE))
)

n_scen_neg <- length(go_decisions_list_neg)
n_meth_neg <- length(go_decisions_list_neg[[1]]$decisions_list)


# -------------------------------------------------------------------
# Test: negateGoDecisions – go_decisions_list must be decision_list
# Input:
#   - go_decisions_list = list().
# Behaviour:
#   - assertClass(go_decisions_list, "decision_list") fails.
# Expectations:
#   - Error mentioning 'go_decisions_list'.
# Why:
#   - Ensures correct object type as input.
# -------------------------------------------------------------------
test_that("negateGoDecisions: go_decisions_list must be a decision_list", {
  expect_error(
    negateGoDecisions(go_decisions_list = list()),
    "go_decisions_list"
  )
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – overall_min_nogos validation
# Input:
#   - overall_min_nogos = -1L, and overall_min_nogos = "foo".
# Behaviour:
#   - Combined check: either 'all' or non-negative integer.
# Expectations:
#   - Both invalid values trigger error mentioning 'overall_min_nogos'.
# Why:
#   - Ensures threshold is meaningful or 'all'.
# -------------------------------------------------------------------
test_that("negateGoDecisions: overall_min_nogos must be 'all' or non-negative integer", {
  expect_error(
    negateGoDecisions(
      go_decisions_list = go_decisions_list_neg,
      overall_min_nogos = -1L
    ),
    "overall_min_nogos"
  )
  
  expect_error(
    negateGoDecisions(
      go_decisions_list = go_decisions_list_neg,
      overall_min_nogos = "foo"
    ),
    "overall_min_nogos"
  )
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – cohort-level entries are logically negated
# Input:
#   - go_decisions_list_neg with multiple scenarios/methods.
# Behaviour:
#   - All non-overall columns are negated (TRUE->FALSE, FALSE->TRUE).
# Expectations:
#   - For each scenario/method, nogo_coh == !go_coh.
# Why:
#   - Core semantics of “negating” a decision_list.
# -------------------------------------------------------------------
test_that("negateGoDecisions: cohort-level entries are logically negated", {
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = "all"
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      go_mat   <- go_decisions_list_neg[[s]]$decisions_list[[m]]
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(go_mat) > 1L) {
        go_coh   <- go_mat[, -1, drop = FALSE]
        nogo_coh <- nogo_mat[, -1, drop = FALSE]
        
        expect_identical(
          nogo_coh,
          !go_coh,
          info = paste("Cohort decisions not negated in scenario", s, "method", m)
        )
      }
    }
  }
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – overall_min_nogos = "all"
# Input:
#   - overall_min_nogos = "all".
# Behaviour:
#   - Threshold internally becomes n_decisions - 1 (all cohort columns).
#   - overall TRUE iff all cohort NoGo columns are TRUE.
# Expectations:
#   - overall column equals row-wise indicator of “all cohorts NoGo”.
# Why:
#   - Confirms the semantics of 'all' for overall NoGo decision.
# -------------------------------------------------------------------
test_that("negateGoDecisions: overall_min_nogos = 'all' means 'all cohorts NoGo'", {
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = "all"
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(nogo_mat) <= 1L) next
      
      n_decisions   <- ncol(nogo_mat)
      nogo_coh      <- nogo_mat[, -1, drop = FALSE]
      expected_over <- apply(
        nogo_coh, 1,
        function(x) sum(x) >= (n_decisions - 1L)
      )
      actual_over   <- as.logical(nogo_mat[, 1])
      
      expect_identical(
        actual_over,
        expected_over
      )
    }
  }
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – overall_min_nogos = 0
# Input:
#   - overall_min_nogos = 0L.
# Behaviour:
#   - sum(x) >= 0 is always TRUE, so overall should be TRUE in all rows.
# Expectations:
#   - overall column is entirely TRUE for all scenarios/methods.
# Why:
#   - Edge-case semantics of minimal threshold.
# -------------------------------------------------------------------
test_that("negateGoDecisions: overall_min_nogos = 0 makes overall always TRUE", {
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = 0L
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(nogo_mat) <= 1L) next
      
      overall_col <- as.logical(nogo_mat[, 1])
      expect_true(
        all(overall_col)
      )
    }
  }
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – numeric overall_min_nogos = k
# Input:
#   - overall_min_nogos = 1L.
# Behaviour:
#   - overall TRUE iff at least k cohort NoGo decisions are TRUE.
# Expectations:
#   - overall column equals row-wise indicator of sum(no-go) >= k.
# Why:
#   - Confirms general threshold logic for k ≥ 1.
# -------------------------------------------------------------------
test_that("negateGoDecisions: numeric overall_min_nogos = k means 'at least k cohorts NoGo'", {
  k <- 1L
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = k
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      if (ncol(nogo_mat) <= 1L) next
      
      nogo_coh      <- nogo_mat[, -1, drop = FALSE]
      expected_over <- apply(nogo_coh, 1, function(x) sum(x) >= k)
      actual_over   <- as.logical(nogo_mat[, 1])
      
      expect_identical(
        actual_over,
        expected_over,
        info = paste("Overall NoGo != '>= k cohorts' for k =", k,
                     "in scenario", s, "method", m)
      )
    }
  }
})


# -------------------------------------------------------------------
# Test: negateGoDecisions – single-column decision matrices
# Input:
#   - decision_list with only one column and no explicit overall logic.
# Behaviour:
#   - All entries are simply negated; no overall recomputation.
# Expectations:
#   - new_mat == !orig_mat.
# Why:
#   - Confirms graceful behaviour when there is no explicit overall column.
# -------------------------------------------------------------------
test_that("negateGoDecisions: single-column decisions are just negated, no overall logic", {
  simple_dl <- list(
    scenario_1 = list(
      decisions_list = list(
        method_1 = matrix(c(TRUE, FALSE, TRUE), ncol = 1)
      )
    )
  )
  class(simple_dl) <- "decision_list"
  
  nogo_simple <- negateGoDecisions(
    go_decisions_list = simple_dl,
    overall_min_nogos = "all"
  )
  
  orig_mat <- simple_dl$scenario_1$decisions_list$method_1
  new_mat  <- nogo_simple$scenario_1$decisions_list$method_1
  
  expect_identical(
    new_mat,
    !orig_mat,
    info = "Single-column decision matrix not fully negated"
  )
})
