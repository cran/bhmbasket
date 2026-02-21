# Tests for applicablePreviousTrials -------------------------------------------

app_prev <- setup_applicable_prev_trials()

scen_initial       <- app_prev$scen_initial
analysis_initial   <- app_prev$analysis_initial
go_initial         <- app_prev$go_initial
scen_next          <- app_prev$scen_next
analysis_with_diff <- app_prev$analysis_with_diff
go_with_diff       <- app_prev$go_with_diff
scen_next_diff     <- app_prev$scen_next_diff
method_names_prev  <- app_prev$method_names_prev
quantiles_prev     <- app_prev$quantiles_prev
n_coh_prev         <- app_prev$n_coh_prev

# ------------------------------------------------------------------
# Test: TRUE when all preconditions for reusing previous trials are met
# Input:
#   - scenario_list = scen_next
#   - method_names = method_names_prev
#   - quantiles = quantiles_prev
#   - n_cohorts = n_coh_prev
#   - calc_differences = NULL
# Behaviour:
#   - No diff columns requested; structure and methods consistent.
# Expectations:
#   - Returns a logical scalar TRUE.
# Why:
#   - Baseline positive case for applicablePreviousTrials().
# ------------------------------------------------------------------
test_that("applicablePreviousTrials returns TRUE when all conditions are met (no differences)", {
  res <- applicablePreviousTrials(
    scenario_list    = scen_next,
    method_names     = method_names_prev,  # "pooled"
    quantiles        = quantiles_prev,
    n_cohorts        = n_coh_prev,
    calc_differences = NULL
  )
  
  expect_true(is.logical(res))
  expect_length(res, 1L)
  expect_true(res)
})

# ------------------------------------------------------------------
# Test: TRUE when diff columns are requested and present in previous_analyses
# Input:
#   - scenario_list = scen_next_diff (contains p_diff_* columns)
#   - calc_differences = matrix(c(3,2)) -> p_diff_32
# Behaviour:
#   - Checks that all requested differences are available in previous quantiles.
# Expectations:
#   - Returns TRUE.
# Why:
#   - Confirms that the function correctly recognises usable previous diffs.
# ------------------------------------------------------------------
test_that("applicablePreviousTrials returns TRUE when required diff columns exist", {
  quantiles_diff <- analysis_with_diff[[1]]$analysis_parameters$quantiles
  calc_diff      <- matrix(c(3, 2), ncol = 2)   # corresponds to "p_diff_32"
  
  res <- applicablePreviousTrials(
    scenario_list    = scen_next_diff,
    method_names     = method_names_prev,
    quantiles        = quantiles_diff,
    n_cohorts        = n_coh_prev,
    calc_differences = calc_diff
  )
  
  expect_true(is.logical(res))
  expect_length(res, 1L)
  expect_true(res)
})

# ------------------------------------------------------------------
# Test: FALSE when at least one required diff column is missing
# Input:
#   - scenario_list = scen_next_diff but with a p_diff_* column removed
#   - calc_differences still asks for that difference
# Behaviour:
#   - Detection of missing p_diff_ij should invalidate reuse.
# Expectations:
#   - Returns FALSE.
# Why:
#   - Ensures robustness to incomplete previous quantile objects.
# ------------------------------------------------------------------
test_that("applicablePreviousTrials returns FALSE when required diff columns are missing", {
  scen_broken <- scen_next_diff
  
  # Remove the diff column from the first method's first trial matrix
  pq      <- scen_broken[[1]]$previous_analyses$post_quantiles
  method1 <- names(pq)[1]
  mat1    <- pq[[method1]][[1]]
  diff_cols <- grepl("^p_diff_", colnames(mat1))
  
  if (any(diff_cols)) {
    mat1 <- mat1[, !diff_cols, drop = FALSE]
    pq[[method1]][[1]] <- mat1
    scen_broken[[1]]$previous_analyses$post_quantiles <- pq
  } else {
    skip("No p_diff_* column to remove")
  }
  
  quantiles_diff <- analysis_with_diff[[1]]$analysis_parameters$quantiles
  calc_diff      <- matrix(c(3, 2), ncol = 2)   # still asking for p_diff_32
  
  res <- applicablePreviousTrials(
    scenario_list    = scen_broken,
    method_names     = method_names_prev,
    quantiles        = quantiles_diff,
    n_cohorts        = n_coh_prev,
    calc_differences = calc_diff
  )
  
  expect_true(is.logical(res))
  expect_length(res, 1L)
  expect_false(res)
})

# ------------------------------------------------------------------
# Test: FALSE when method_names differ across scenarios
# Input:
#   - Two scenarios with different names in post_quantiles
# Behaviour:
#   - Applicable only if all scenarios share identical method sets.
# Expectations:
#   - Returns FALSE.
# Why:
#   - Prevents reusing previous trials when method structure is inconsistent.
# ------------------------------------------------------------------
test_that("applicablePreviousTrials returns FALSE when method_names differ across scenarios", {
  # Copy scen_next to two scenarios, then break method-name consistency in scenario_2
  scen_multi <- scen_next
  scen_multi$scenario_2 <- scen_multi$scenario_1
  names(scen_multi) <- c("scenario_1", "scenario_2")
  class(scen_multi) <- "scenario_list"
  
  names(scen_multi$scenario_2$previous_analyses$post_quantiles) <-
    paste0("alt_", names(scen_multi$scenario_2$previous_analyses$post_quantiles))
  
  res <- applicablePreviousTrials(
    scenario_list    = scen_multi,
    method_names     = method_names_prev,
    quantiles        = quantiles_prev,
    n_cohorts        = n_coh_prev,
    calc_differences = NULL
  )
  
  expect_true(is.logical(res))
  expect_length(res, 1L)
  expect_false(res)
})


# Tests for calcDiffsMCMC ------------------------------------------------------

# ------------------------------------------------------------------
# Test: calcDiffsMCMC adds p_diff_ij columns with correct values
# Input:
#   - posterior_samples with p_1, p_2, p_3, mu
#   - calc_differences = (1,2) and (3,1)
# Behaviour:
#   - Adds p_diff_12 = p_1 - p_2, p_diff_31 = p_3 - p_1.
# Expectations:
#   - Original columns retained; new diff columns added with correct numeric values.
# Why:
#   - Validates difference computation used downstream in decision rules.
# ------------------------------------------------------------------
test_that("calcDiffsMCMC: adds correctly named difference columns with correct values", {
  posterior_samples <- cbind(
    p_1 = c(0.1, 0.2),
    p_2 = c(0.5, 0.6),
    p_3 = c(0.8, 0.9),
    mu  = c(1.0, 2.0)
  )
  
  calc_differences <- rbind(
    c(1, 2),
    c(3, 1)
  )
  
  out <- calcDiffsMCMC(
    posterior_samples = posterior_samples,
    calc_differences  = calc_differences
  )
  
  expect_true(all(c("p_1", "p_2", "p_3", "mu") %in% colnames(out)))
  expect_true(all(c("p_diff_12", "p_diff_31") %in% colnames(out)))
  
  expected_diff_12 <- posterior_samples[, "p_1"] - posterior_samples[, "p_2"]
  expected_diff_31 <- posterior_samples[, "p_3"] - posterior_samples[, "p_1"]
  
  expect_equal(out[, "p_diff_12"], expected_diff_12)
  expect_equal(out[, "p_diff_31"], expected_diff_31)
})


# Tests for performJags --------------------------------------------------------

model_file <- tempfile(fileext = ".bug")
writeLines(
  c(
    "model {",
    "  theta ~ dbeta(1, 1)",
    "  for (i in 1:N) {",
    "    y[i] ~ dbern(theta)",
    "  }",
    "}"
  ),
  con = model_file
)

data_list <- list(
  N = 20L,
  y = c(rep(1L, 14), rep(0L, 6))
)

# ------------------------------------------------------------------
# Test: performJags runs a simple Bernoulli model and returns reasonable samples
# Input:
#   - Beta(1,1) prior, 14 successes / 6 failures
# Behaviour:
#   - Posterior is Beta(15,7); posterior mean ~ 15/22.
# Expectations:
#   - Returns matrix with column "theta", expected number of rows,
#     and sample mean close to 15/22.
# Why:
#   - Integration test for JAGS wrapper and posterior extraction.
# ------------------------------------------------------------------
test_that("performJags runs a simple Bernoulli model and returns a sensible sample matrix", {
  skip_if_not_installed("rjags")
  
  n_chains <- 2L
  n_iter   <- 1000L
  n_burnin <- 100L
  
  samples <- performJags(
    data               = data_list,
    parameters_to_save = "theta",
    model_file         = model_file,
    n_chains           = n_chains,
    n_iter             = n_iter,
    n_burnin           = n_burnin
  )
  
  expect_true(is.matrix(samples))
  expect_identical(colnames(samples), "theta")
  
  expected_rows <- n_chains * (n_iter - n_burnin)
  expect_equal(nrow(samples), expected_rows)
  
  # Given 14 successes and 6 failures with a Beta(1,1) prior, the posterior is Beta(15, 7).
  # The analytic posterior mean is 15 / (15 + 7) = 15/22.
  expect_equal(mean(samples[, "theta"]), 15/22, tolerance = 0.01)
})

# ------------------------------------------------------------------
# Test: performJags handles n_burnin = 0
# Input:
#   - Same model/data, with n_burnin = 0.
# Behaviour:
#   - Uses all n_iter samples in the posterior.
# Expectations:
#   - Number of rows equals n_chains * n_iter;
#     sample mean still close to 15/22.
# Why:
#   - Verifies edge case where no explicit burn-in is requested.
# ------------------------------------------------------------------
test_that("performJags handles n_burning equals 0 correctly", {
  skip_if_not_installed("rjags")
  
  n_chains <- 2L
  n_iter   <- 1000L
  n_burnin <- 0
  
  samples <- performJags(
    data               = data_list,
    parameters_to_save = "theta",
    model_file         = model_file,
    n_chains           = n_chains,
    n_iter             = n_iter,
    n_burnin           = n_burnin
  )
  
  expect_true(is.matrix(samples))
  expect_identical(colnames(samples), "theta")
  
  expected_rows <- n_chains * (n_iter - n_burnin)
  expect_equal(nrow(samples), expected_rows)
  
  expect_equal(mean(samples[, "theta"]), 15/22, tolerance = 0.01)
})


# Tests for getPosteriors ------------------------------------------------------

# ------------------------------------------------------------------
# Test: getPosteriors basic behaviour and name cleaning (no exch weights)
# Input:
#   - Berry model prepared via prepareAnalysis, single trial scenario.
# Behaviour:
#   - Calls JAGS, returns posterior samples matrix with cleaned column names
#     (no brackets, no w_*, no "exch").
# Expectations:
#   - Matrix result, finite entries, column names stripped of indices,
#     no exNEX-specific columns.
# Why:
#   - Ensures generic posterior extraction cleans names consistently.
# ------------------------------------------------------------------
test_that("getPosteriors: basic posterior sampling and name cleaning (no exch weights)", {
  skip_if_not_installed("rjags")
  
  set.seed(123)
  
  scen_list <- simulateScenarios(
    n_subjects_list     = list(c(10, 20)),
    response_rates_list = list(c(0.5, 0.6)),
    n_trials            = 1
  )
  scen1 <- scen_list$scenario_1
  
  target_rates <- c(0.5, 0.5)
  priors       <- getPriorParameters(
    method_names = "berry",
    target_rates = target_rates
  )
  
  prep <- prepareAnalysis(
    method_name      = "berry",
    target_rates     = target_rates,
    prior_parameters = priors[["berry"]]
  )
  
  j_data <- prep$j_data
  j_data$r <- as.numeric(scen1$n_responders[1, ])
  j_data$n <- as.numeric(scen1$n_subjects[1, ])
  
  post <- getPosteriors(
    j_parameters      = prep$j_parameters,
    j_model_file      = prep$j_model_file,
    j_data            = j_data,
    n_mcmc_iterations = 20
  )
  
  expect_true(is.matrix(post))
  expect_gt(nrow(post), 0)
  expect_true(all(is.finite(post)))
  
  expect_false(any(grepl("\\[", colnames(post))))
  expect_false(any(grepl("\\]", colnames(post))))
  
  expect_false(any(grepl("^w_", colnames(post))))
  expect_false(any(grepl("exch", colnames(post))))
})

# ------------------------------------------------------------------
# Test: getPosteriors for exNEX renames exch weights and strips extras
# Input:
#   - exNEX analysis prepared via prepareAnalysis on a small scenario.
# Behaviour:
#   - Columns for exchangeability weights are renamed to w_j*;
#     exNEX internal parameters and redundant labels removed.
# Expectations:
#   - Posterior matrix, finite; no "exch" in names;
#     exactly J cohort w_* columns; no ",1" suffixes.
# Why:
#   - Confirms special handling of mixture weights in exNEX models.
# ------------------------------------------------------------------
test_that("getPosteriors: exNEX exchangeability weights are renamed to w_j and extras dropped", {
  skip_if_not_installed("rjags")
  
  set.seed(456)
  
  scen_list <- simulateScenarios(
    n_subjects_list     = list(c(10, 20)),
    response_rates_list = list(c(0.6, 0.7)),
    n_trials            = 1
  )
  scen1 <- scen_list$scenario_1
  
  target_rates <- c(0.5, 0.5)
  priors       <- getPriorParameters(
    method_names = "exnex",
    target_rates = target_rates
  )
  
  prep <- prepareAnalysis(
    method_name      = "exnex",
    target_rates     = target_rates,
    prior_parameters = priors[["exnex"]]
  )
  
  j_data <- prep$j_data
  j_data$r <- as.numeric(scen1$n_responders[1, ])
  j_data$n <- as.numeric(scen1$n_subjects[1, ])
  
  post <- getPosteriors(
    j_parameters      = prep$j_parameters,
    j_model_file      = prep$j_model_file,
    j_data            = j_data,
    n_mcmc_iterations = 20
  )
  
  expect_true(is.matrix(post))
  expect_gt(nrow(post), 0)
  expect_true(all(is.finite(post)))
  
  expect_false(any(grepl("exch", colnames(post))))
  
  w_cols <- grep("^w_", colnames(post))
  expect_equal(length(w_cols), length(j_data$n))
  
  expect_false(any(grepl(",1", colnames(post))))
})


# Tests for getPostQuantiles ---------------------------------------------------

# ------------------------------------------------------------------
# Test: pooled backend with calc_differences adds zero p_diff_* columns
# Input:
#   - scenario_data with one trial / 2 cohorts
#   - method_name = "pooled", j_data = Beta prior, calc_differences for 1↔2
# Behaviour:
#   - Pooled backend duplicates same marginal posterior for p_1 and p_2;
#     p_diff_* columns are filled with zeros by design.
# Expectations:
#   - One matrix in list; exact column/row names;
#     p_diff_* = 0; pooled Beta quantiles/mean/SD match analytic forms.
# Why:
#   - Verifies deterministic pooled implementation and diff-column handling.
# ------------------------------------------------------------------
test_that("getPostQuantiles (pooled): calc_differences adds zero columns with correct names and shapes", {
  n_subjects_mat   <- matrix(c(10, 20), nrow = 1)
  n_responders_mat <- matrix(c(3,  5),  nrow = 1)
  
  scenario_data <- list(
    n_subjects   = n_subjects_mat,
    n_responders = n_responders_mat
  )
  
  j_data    <- list(a = 1, b = 1)
  quantiles <- c(0.025, 0.5, 0.975)
  
  # Differences to compute (pooled backend fills them with zeros)
  calc_differences <- rbind(c(1, 2), c(2, 1))
  
  out <- getPostQuantiles(
    method_name       = "pooled",
    quantiles         = quantiles,
    scenario_data     = scenario_data,
    calc_differences  = calc_differences,
    j_parameters      = NULL,
    j_model_file      = NULL,
    j_data            = j_data,
    n_mcmc_iterations = 20,
    save_path         = NULL,
    save_trial        = NULL
  )
  
  expect_type(out, "list")
  expect_length(out, 1)
  
  mat <- out[[1]]
  expect_true(is.matrix(mat))
  
  expect_identical(
    colnames(mat),
    c("p_1", "p_2", "p_diff_12", "p_diff_21")
  )
  expect_identical(
    rownames(mat),
    c("2.5%", "50%", "97.5%", "Mean", "SD")
  )
  
  expect_true(all(mat[, c("p_diff_12", "p_diff_21")] == 0))
  
  r_tot   <- sum(n_responders_mat[1, ])
  n_tot   <- sum(n_subjects_mat[1, ])
  shape_1 <- j_data$a + r_tot
  shape_2 <- j_data$b + (n_tot - r_tot)
  
  expected_q  <- stats::qbeta(quantiles, shape1 = shape_1, shape2 = shape_2)
  expected_mu <- shape_1 / (shape_1 + shape_2)
  expected_sd <- sqrt((shape_1 * shape_2) / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1)))
  
  expect_equal(unname(mat[c("2.5%", "50%", "97.5%"), "p_1"]), expected_q, tolerance = 1e-12)
  expect_equal(unname(mat[c("2.5%", "50%", "97.5%"), "p_2"]), expected_q, tolerance = 1e-12)
  expect_equal(unname(mat["Mean", c("p_1", "p_2")]), rep(expected_mu, 2), tolerance = 1e-12)
  expect_equal(unname(mat["SD",   c("p_1", "p_2")]), rep(expected_sd, 2), tolerance = 1e-12)
})

# ------------------------------------------------------------------
# Test: stratified backend with multiple trials and calc_differences
# Input:
#   - 2 trials x 2 cohorts; method_name = "stratified"
#   - calc_differences = (1,2)
# Behaviour:
#   - Returns per-trial quantile matrices with p_1, p_2 and p_diff_12,
#     finite values for marginal posterior summaries.
# Expectations:
#   - Length of out equals number of trials; each matrix has expected
#     row/column labels and finite p_j entries.
# Why:
#   - Confirms stratified Beta backend and diff logic over multiple trials.
# ------------------------------------------------------------------
test_that("getPostQuantiles (stratified): multiple trials and calc_differences produce p_j and p_diff_*", {
  skip_if_not_installed("foreach")
  skip_if_not_installed("doRNG")
  
  foreach::registerDoSEQ()
  
  n_subjects_mat <- rbind(
    c(10, 10),
    c(20, 20)
  )
  n_responders_mat <- rbind(
    c(3, 5),
    c(6, 10)
  )
  
  scenario_data <- list(
    n_subjects   = n_subjects_mat,
    n_responders = n_responders_mat
  )
  
  j_data <- list(
    a_j = c(1, 1),
    b_j = c(1, 1)
  )
  
  quantiles        <- c(0.25, 0.5, 0.75)
  calc_differences <- matrix(c(1, 2), ncol = 2)
  
  out <- getPostQuantiles(
    method_name       = "stratified",
    quantiles         = quantiles,
    scenario_data     = scenario_data,
    calc_differences  = calc_differences,
    j_parameters      = NULL,
    j_model_file      = NULL,
    j_data            = j_data,
    n_mcmc_iterations = 20,
    save_path         = NULL,
    save_trial        = NULL
  )
  
  expect_type(out, "list")
  expect_length(out, 2)
  
  for (mat in out) {
    expect_true(is.matrix(mat))
    
    expect_setequal(
      colnames(mat),
      c("p_1", "p_2", "p_diff_12")
    )
    
    expect_setequal(
      rownames(mat),
      c(paste0(quantiles * 100, "%"), "Mean", "SD")
    )
    
    expect_true(all(is.finite(mat[, c("p_1", "p_2")])))
    
    diff_col <- mat[, "p_diff_12"]
    expect_true(is.numeric(diff_col))
  }
})


# Tests for prepareAnalysis ----------------------------------------------------

# ------------------------------------------------------------------
# Test: berry branch builds j_data and parameter names
# Input:
#   - target_rates = c(0.5, 0.7); priors from getPriorParameters("berry")
# Behaviour:
#   - j_data fields populated from priors and target_rates;
#     J equals number of cohorts; JAGS parameters and model file exist.
# Expectations:
#   - mean_mu, precision_mu, precision_tau, p_t, J consistent;
#     j_parameters = c("p","mu","tau"); model file exists.
# Why:
#   - Validates mapping from high-level Berry prior to JAGS inputs.
# ------------------------------------------------------------------
test_that("prepareAnalysis: berry branch builds correct j_data and parameters", {
  target_rates <- c(0.5, 0.7)
  priors_list  <- getPriorParameters(method_names = "berry", target_rates = target_rates)
  priors_berry <- priors_list[["berry"]]
  
  prep <- prepareAnalysis(
    method_name      = "berry",
    prior_parameters = priors_berry,
    target_rates     = target_rates
  )
  
  j_data <- prep$j_data
  
  expect_equal(j_data$mean_mu,       priors_berry$mu_mean)
  expect_equal(j_data$precision_mu,  priors_berry$mu_sd^-2)
  expect_equal(j_data$precision_tau, priors_berry$tau_scale^-2)
  expect_equal(j_data$p_t,           target_rates)
  expect_equal(j_data$J,             length(target_rates))
  
  expect_equal(prep$j_parameters, c("p", "mu", "tau"))
  expect_true(is.character(prep$j_model_file))
  expect_true(file.exists(prep$j_model_file))
})

# ------------------------------------------------------------------
# Test: exnex branch constructs mixture hyperparameters and pMix
# Input:
#   - target_rates = c(0.4,0.6,0.8); priors from getPriorParameters("exnex")
# Behaviour:
#   - j_data contains Nexch, Nmix, Nstrata, mu_mean, mu_prec,
#     tau_HN_scale, nex_mean, nex_prec and pMix that sums to 1.
# Expectations:
#   - Dimensions consistent with priors; pMix length Nexch+1 and sums to 1.
# Why:
#   - Ensures exNEX prior is correctly translated into JAGS hyperparameters.
# ------------------------------------------------------------------
test_that("prepareAnalysis: exnex branch builds mixture priors and pMix", {
  target_rates <- c(0.4, 0.6, 0.8)
  priors_list  <- getPriorParameters(method_names = "exnex", target_rates = target_rates)
  priors_exnex <- priors_list[["exnex"]]
  
  prep <- prepareAnalysis(
    method_name      = "exnex",
    prior_parameters = priors_exnex,
    target_rates     = target_rates
  )
  
  j_data <- prep$j_data
  
  Nexch   <- length(priors_exnex$mu_mean)
  Nstrata <- length(priors_exnex$mu_j)
  
  expect_equal(j_data$Nexch,   Nexch)
  expect_equal(j_data$Nmix,    Nexch + 1L)
  expect_equal(j_data$Nstrata, Nstrata)
  
  expect_equal(j_data$mu_mean, priors_exnex$mu_mean)
  expect_equal(j_data$mu_prec, priors_exnex$mu_sd^-2)
  
  expect_equal(j_data$tau_HN_scale,
               rep(priors_exnex$tau_scale, Nexch))
  
  expect_equal(j_data$nex_mean, priors_exnex$mu_j)
  expect_equal(j_data$nex_prec, priors_exnex$tau_j^-2)
  
  expect_length(j_data$pMix, Nexch + 1L)
  expect_equal(sum(j_data$pMix), 1, tolerance = 1e-12)
  
  expect_equal(prep$j_parameters, c("p", "mu", "tau", "exch"))
  expect_true(file.exists(prep$j_model_file))
})

# ------------------------------------------------------------------
# Test: exnex_adj branch includes p_target and uses adjusted model
# Input:
#   - target_rates = c(0.5,0.6); priors from getPriorParameters("exnex_adj")
# Behaviour:
#   - j_data$p_target equals target_rates; exnex_adj JAGS model used;
#     parameters include "exch".
# Expectations:
#   - p_target set; j_parameters = c("p","mu","tau","exch");
#     model file path contains "exnex_adj".
# Why:
#   - Checks adjusted exNEX prior (centered) is wired correctly.
# ------------------------------------------------------------------
test_that("prepareAnalysis: exnex_adj branch includes p_target and uses exnex_adj model", {
  target_rates <- c(0.5, 0.6)
  priors_list  <- getPriorParameters(method_names = "exnex_adj", target_rates = target_rates)
  priors       <- priors_list[["exnex_adj"]]
  
  prep <- prepareAnalysis(
    method_name      = "exnex_adj",
    prior_parameters = priors,
    target_rates     = target_rates
  )
  
  j_data <- prep$j_data
  
  expect_equal(j_data$p_target, target_rates)
  expect_equal(prep$j_parameters, c("p", "mu", "tau", "exch"))
  expect_true(grepl("exnex_adj", prep$j_model_file))
  expect_true(file.exists(prep$j_model_file))
})

# ------------------------------------------------------------------
# Test: stratified and pooled branches pass priors through and use dummy JAGS info
# Input:
#   - target_rates = (0.5,0.5); priors for "stratified" and "pooled"
# Behaviour:
#   - j_data equals prior structures;
#     j_model_file and j_parameters are placeholders.
# Expectations:
#   - Dummy values returned; priors unchanged.
# Why:
#   - These methods do not require JAGS; the interface must still be consistent.
# ------------------------------------------------------------------
test_that("prepareAnalysis: stratified and pooled use dummy JAGS info and pass priors through", {
  target_rates <- c(0.5, 0.5)
  priors_list  <- getPriorParameters(
    method_names = c("stratified", "pooled"),
    target_rates = target_rates
  )
  
  prep_strat <- prepareAnalysis(
    method_name      = "stratified",
    prior_parameters = priors_list[["stratified"]],
    target_rates     = target_rates
  )
  expect_equal(prep_strat$j_model_file, "dummy path to JAGS model")
  expect_equal(prep_strat$j_parameters, "dummy JAGS parameters")
  expect_identical(prep_strat$j_data, priors_list[["stratified"]])
  
  prep_pooled <- prepareAnalysis(
    method_name      = "pooled",
    prior_parameters = priors_list[["pooled"]],
    target_rates     = target_rates
  )
  expect_equal(prep_pooled$j_model_file, "dummy path to JAGS model")
  expect_equal(prep_pooled$j_parameters, "dummy JAGS parameters")
  expect_identical(prep_pooled$j_data, priors_list[["pooled"]])
})

# ------------------------------------------------------------------
# Test: prepareAnalysis errors for unsupported method_name
# Input:
#   - method_name = "invalid_method"
# Behaviour:
#   - Guard clause should trigger with informative message.
# Expectations:
#   - Error containing "method_name must be one of".
# Why:
#   - Protects against typos / unsupported analysis methods.
# ------------------------------------------------------------------
test_that("prepareAnalysis: invalid method_name throws an error", {
  expect_error(
    prepareAnalysis(method_name = "invalid_method"),
    "method_name must be one of"
  )
})


# Tests for getUniqueRows ------------------------------------------------------

# ------------------------------------------------------------------
# Test: getUniqueRows returns unique row combinations with full column set
# Input:
#   - Matrix with duplicate rows intermixed.
# Behaviour:
#   - Removes duplicates, preserving all columns.
# Expectations:
#   - ncol(out) == ncol(input); rows correspond to unique() of data.frame(mat).
# Why:
#   - Low-level helper for deduplicating trial patterns.
# ------------------------------------------------------------------
test_that("getUniqueRows: returns unique row combinations with correct columns", {
  mat <- rbind(
    c(1, 2),
    c(1, 2),
    c(2, 3),
    c(2, 3),
    c(4, 5)
  )
  
  out <- getUniqueRows(mat)
  
  expect_equal(ncol(out), ncol(mat))
  
  out_df <- as.data.frame(out)
  expect_equal(nrow(out_df), nrow(unique(out_df)))
  
  exp_df <- as.data.frame(unique(mat))
  
  out_ord <- out_df[do.call(order, out_df), , drop = FALSE]
  exp_ord <- exp_df[do.call(order, exp_df), , drop = FALSE]
  
  expect_equal(unname(out_ord), unname(exp_ord))
})


# Tests for getUniqueTrials ----------------------------------------------------

# ------------------------------------------------------------------
# Test: getUniqueTrials combines scenarios and returns unique responder/subject/go rows
# Input:
#   - scenario_list with controlled n_responders/n_subjects and go_decisions.
# Behaviour:
#   - Stacks all trials, attaches go_flag, and returns unique rows.
# Expectations:
#   - 5 columns (2 responders + 2 subjects + 1 go_flag);
#     no duplicate rows; equal to unique() of combined matrix.
# Why:
#   - Verifies mapping from scenario_list to a unique-trial representation.
# ------------------------------------------------------------------
test_that("getUniqueTrials: combines scenarios and returns unique responder/subject/go rows", {
  set.seed(123)
  
  scenario_list <- simulateScenarios(
    n_subjects_list     = list(c(10, 10),
                               c(10, 10)),
    response_rates_list = list(c(0.1, 0.2),
                               c(0.1, 0.2)),
    n_trials            = 2
  )
  
  # scenario_1: two identical trials (1,2 | 10,10)
  scenario_list$scenario_1$n_responders <- rbind(
    c(1, 2),
    c(1, 2)
  )
  scenario_list$scenario_1$n_subjects <- rbind(
    c(10, 10),
    c(10, 10)
  )
  
  # scenario_2: one duplicate row (1,2 | 10,10) and one distinct row (2,1 | 10,10)
  scenario_list$scenario_2$n_responders <- rbind(
    c(1, 2),
    c(2, 1)
  )
  scenario_list$scenario_2$n_subjects <- rbind(
    c(10, 10),
    c(10, 10)
  )
  
  scenario_list$scenario_1$previous_analyses <- list(
    go_decisions = cbind(c(1, 0))
  )
  scenario_list$scenario_2$previous_analyses <- list(
    go_decisions = cbind(c(1, 1))
  )
  
  class(scenario_list) <- "scenario_list"
  
  out <- getUniqueTrials(scenario_list)
  
  expect_equal(ncol(out), 5)
  
  out_df <- as.data.frame(out)
  expect_equal(nrow(out_df), nrow(unique(out_df)))
  
  all_resp <- do.call(rbind, lapply(scenario_list, function(x) x$n_responders))
  all_subj <- do.call(rbind, lapply(scenario_list, function(x) x$n_subjects))
  all_go   <- do.call(
    rbind,
    lapply(scenario_list, function(x) x$previous_analyses$go_decisions)
  )[, 1]
  
  combined <- cbind(all_resp, all_subj, go_flag = all_go)
  
  exp_df <- as.data.frame(unique(combined))
  
  out_ord <- out_df[do.call(order, out_df), , drop = FALSE]
  exp_ord <- exp_df[do.call(order, exp_df), , drop = FALSE]
  
  expect_equal(
    unname(as.matrix(out_ord)),
    unname(as.matrix(exp_ord))
  )
})


# Tests for mapUniqueTrials ----------------------------------------------------

set.seed(123)

trial_data <- createTrial(
  n_subjects   = c(10, 10),
  n_responders = c(1, 1))

analysis_list <- performAnalyses(
  scenario_list      = trial_data,
  target_rates       = rep(0.5, 2),
  method_names       = "berry",
  n_mcmc_iterations  = 20)

trials_unique <- getUniqueTrials(trial_data)
n_cohorts     <- (ncol(trials_unique) - 1L) / 2

applicable_previous_trials <- applicablePreviousTrials(
  scenario_list    = trial_data,
  method_names     = "berry",
  quantiles        = analysis_list$scenario_1$analysis_parameters$quantiles,
  n_cohorts        = n_cohorts,
  calc_differences = NULL)

expect_false(applicable_previous_trials)

if (applicable_previous_trials) {
  calc_trial_indices <- trials_unique[, ncol(trials_unique)] > 0
} else {
  calc_trial_indices <- rep(TRUE, nrow(trials_unique))
}

trials_unique_calc <- trials_unique[calc_trial_indices, -ncol(trials_unique)]
n_responders       <- trials_unique_calc[, seq_len(n_cohorts)]
n_subjects         <- trials_unique_calc[, seq_len(n_cohorts) + n_cohorts]

method_quantiles_list        <- vector(mode = "list", length = 1)
names(method_quantiles_list) <- "berry"

prior_parameters_list <- getPriorParameters(
  method_names = "berry",
  target_rates = rep(0.5, 2))

prepare_analysis <- prepareAnalysis(
  method_name       = "berry",
  target_rates      = rep(0.5, 2),
  prior_parameters  = prior_parameters_list[["berry"]])

method_quantiles_list[["berry"]] <- getPostQuantiles(
  method_name       = "berry",
  quantiles         = analysis_list$scenario_1$analysis_parameters$quantiles,
  scenario_data     = list(n_subjects   = n_subjects,
                           n_responders = n_responders),
  calc_differences  = NULL,
  j_parameters      = prepare_analysis$j_parameters,
  j_model_file      = prepare_analysis$j_model_file,
  j_data            = prepare_analysis$j_data,
  n_mcmc_iterations = 20,
  save_path         = NULL,
  save_trial        = NULL)

scenario_method_quantiles_list <- mapUniqueTrials(
  scenario_list              = trial_data,
  method_quantiles_list      = method_quantiles_list,
  trials_unique_calc         = trials_unique_calc,
  applicable_previous_trials = applicable_previous_trials)

# ------------------------------------------------------------------
# Test: with no previous trials, mapUniqueTrials copies unique quantiles per scenario
# Input:
#   - trial_data with a single scenario and no applicable_previous_trials.
# Behaviour:
#   - All trials are treated as new; their matrices equal method_quantiles_list entries.
# Expectations:
#   - scenario_1 exists; "berry" method present; each trial's matrix identical
#     to precomputed method_quantiles_list$berry[[i]].
# Why:
#   - Baseline behaviour for the mapping when no reuse is possible.
# ------------------------------------------------------------------
test_that("mapUniqueTrials: without previous trials, maps unique trial quantiles back per scenario", {
  foreach::registerDoSEQ()
  
  out <- scenario_method_quantiles_list
  
  expect_type(out, "list")
  expect_identical(names(out), "scenario_1")
  
  scen1 <- out[["scenario_1"]]
  expect_true(is.list(scen1))
  expect_true("berry" %in% names(scen1))
  
  for (i in seq_along(method_quantiles_list$berry)) {
    expect_identical(scen1$berry[[i]], method_quantiles_list$berry[[i]])
  }
})

# ------------------------------------------------------------------
# Test: with previous trials, only GO trials are updated from hash tables
# Input:
#   - One scenario with 2 trials: trial 1 GO, trial 2 NoGo.
#   - previous_analyses$post_quantiles present; new values provided for GO pattern.
# Behaviour:
#   - GO trial gets updated quantiles; NoGo trial retains its old ones.
# Expectations:
#   - scen1$berry[[1]] == new_q_go; scen1$berry[[2]] == original old_q_list[[2]].
# Why:
#   - Validates fine-grained mapping logic between unique patterns and trials.
# ------------------------------------------------------------------
test_that("mapUniqueTrials: with previous trials, only GO trials are updated from hash tables", {
  foreach::registerDoSEQ()
  set.seed(456)
  
  scene <- getScenario(
    n_subjects     = c(10, 10),
    response_rates = c(0.1, 0.2),
    n_trials       = 2
  )
  
  trial_data <- list(scenario_1 = scene)
  class(trial_data) <- "scenario_list"
  if (is.null(trial_data$scenario_1$scenario_number)) {
    trial_data$scenario_1$scenario_number <- 1L
  }
  
  go_mat <- trial_data$scenario_1$previous_analyses$go_decisions
  go_mat[, 1] <- c(TRUE, FALSE)
  trial_data$scenario_1$previous_analyses$go_decisions <- go_mat
  
  n_trials   <- trial_data$scenario_1$n_trials
  n_cohorts  <- ncol(trial_data$scenario_1$n_subjects)
  quantiles  <- c(0.025, 0.5, 0.975)
  
  prior_parameters_list <- getPriorParameters(
    method_names = "berry",
    target_rates = rep(0.5, n_cohorts)
  )
  
  prep <- prepareAnalysis(
    method_name      = "berry",
    prior_parameters = prior_parameters_list[["berry"]],
    target_rates     = rep(0.5, n_cohorts)
  )
  
  old_q_list <- vector("list", length = n_trials)
  for (i in seq_len(n_trials)) {
    n_subj_i <- trial_data$scenario_1$n_subjects[i, , drop = FALSE]
    n_resp_i <- trial_data$scenario_1$n_responders[i, , drop = FALSE]
    
    old_q_list[[i]] <- getPostQuantiles(
      method_name       = "berry",
      quantiles         = quantiles,
      scenario_data     = list(
        n_subjects   = n_subj_i,
        n_responders = n_resp_i
      ),
      calc_differences  = NULL,
      j_parameters      = prep$j_parameters,
      j_model_file      = prep$j_model_file,
      j_data            = prep$j_data,
      n_mcmc_iterations = 20,
      save_path         = NULL,
      save_trial        = NULL
    )[[1]]
  }
  
  trial_data$scenario_1$previous_analyses$post_quantiles <- list(
    berry = old_q_list
  )
  
  trials_unique <- getUniqueTrials(trial_data)
  n_cohorts_u   <- (ncol(trials_unique) - 1L) / 2L
  
  applicable_previous_trials <- TRUE
  
  calc_trial_indices <- trials_unique[, ncol(trials_unique)] > 0
  trials_unique_calc <- trials_unique[calc_trial_indices, -ncol(trials_unique), drop = FALSE]
  
  new_q_go <- old_q_list[[1]] + 1
  
  method_quantiles_list <- list(
    berry = list(new_q_go)
  )
  
  out <- mapUniqueTrials(
    scenario_list              = trial_data,
    method_quantiles_list      = method_quantiles_list,
    trials_unique_calc         = trials_unique_calc,
    applicable_previous_trials = applicable_previous_trials
  )
  
  expect_type(out, "list")
  expect_identical(names(out), "scenario_1")
  
  scen1 <- out[["scenario_1"]]
  expect_true(is.list(scen1))
  expect_true("berry" %in% names(scen1))
  
  expect_identical(scen1$berry[[1]], new_q_go)
  expect_identical(scen1$berry[[2]], old_q_list[[2]])
})

## Test: mapUniqueTrials with pooled backend preserves naive per-trial quantiles
## Input:
##   - One scenario with several trials and 2 cohorts (generated via getScenario()).
##   - trials_unique_calc derived from getUniqueTrials(), with
##     applicable_previous_trials = FALSE so that all unique patterns are treated
##     as newly analysed.
##   - method_quantiles_list for method "pooled" containing posterior quantile
##     matrices computed only for the unique trial patterns using getPostQuantiles().
##
## Naive reference:
##   - For each original trial i (row i of n_responders, n_subjects), we call
##     getPostQuantiles(method_name = "pooled", ...) directly and store the
##     resulting single-trial quantile matrix in naive_list[[i]].
##
## Expected output:
##   - mapUniqueTrials returns a list with one element named "scenario_1".
##   - scenario_1$pooled is a list of length equal to the original number of trials.
##   - For each trial i, scenario_1$pooled[[i]] is exactly equal to
##     naive_list[[i]] (the per-trial quantile matrix from the naive analysis).
##
## Why this test:
##   - For the pooled design, getPostQuantiles() is deterministic given
##     (n_responders, n_subjects) and j_data.
##   - mapUniqueTrials must use hashing + unique trial patterns to avoid
##     redundant computations, but still reproduce the same per-trial quantiles
##     as the naive approach that analyses each trial separately.

test_that("mapUniqueTrials: pooled backend preserves naive per-trial Beta quantiles", {
  skip_if_not_installed("foreach")
  foreach::registerDoSEQ()
  set.seed(123)
  
  scene <- getScenario(
    n_subjects     = c(10, 20),
    response_rates = c(0.4, 0.6),
    n_trials       = 10
  )
  
  if (is.null(scene$scenario_number)) {
    scene$scenario_number <- 1L
  }
  
  scenario_list <- list(scenario_1 = scene)
  class(scenario_list) <- "scenario_list"
  
  n_subj   <- scene$n_subjects
  n_resp   <- scene$n_responders
  n_trials <- nrow(n_subj)
  n_coh    <- ncol(n_subj)
  
  trials_unique <- getUniqueTrials(scenario_list)
  n_coh_u       <- (ncol(trials_unique) - 1L) / 2L
  expect_equal(n_coh_u, n_coh)
  
  applicable_previous_trials <- FALSE
  
  calc_trial_indices <- rep(TRUE, nrow(trials_unique))
  trials_unique_calc <- trials_unique[calc_trial_indices, -ncol(trials_unique), drop = FALSE]
  
  n_resp_unique <- trials_unique_calc[, seq_len(n_coh),         drop = FALSE]
  n_subj_unique <- trials_unique_calc[, seq_len(n_coh) + n_coh, drop = FALSE]
  
  j_data    <- list(a = 1, b = 1)
  quantiles <- c(0.025, 0.5, 0.975)
  
  unique_quantiles <- getPostQuantiles(
    method_name       = "pooled",
    quantiles         = quantiles,
    scenario_data     = list(
      n_subjects   = n_subj_unique,
      n_responders = n_resp_unique
    ),
    calc_differences  = NULL,
    j_parameters      = NULL,
    j_model_file      = NULL,
    j_data            = j_data,
    n_mcmc_iterations = 20,
    save_path         = NULL,
    save_trial        = NULL
  )
  
  naive_list <- vector("list", length = n_trials)
  for (i in seq_len(n_trials)) {
    n_subj_i <- n_subj[i, , drop = FALSE]
    n_resp_i <- n_resp[i, , drop = FALSE]
    
    naive_list[[i]] <- getPostQuantiles(
      method_name       = "pooled",
      quantiles         = quantiles,
      scenario_data     = list(
        n_subjects   = n_subj_i,
        n_responders = n_resp_i
      ),
      calc_differences  = NULL,
      j_parameters      = NULL,
      j_model_file      = NULL,
      j_data            = j_data,
      n_mcmc_iterations = 20,
      save_path         = NULL,
      save_trial        = NULL
    )[[1]]
  }
  
  method_quantiles_list <- list(
    pooled = unique_quantiles
  )
  
  out <- mapUniqueTrials(
    scenario_list              = scenario_list,
    method_quantiles_list      = method_quantiles_list,
    trials_unique_calc         = trials_unique_calc,
    applicable_previous_trials = applicable_previous_trials
  )
  
  expect_type(out, "list")
  expect_identical(names(out), "scenario_1")
  
  scen1_out <- out[["scenario_1"]]
  expect_true(is.list(scen1_out))
  expect_true("pooled" %in% names(scen1_out))
  expect_equal(length(scen1_out$pooled), n_trials)
  
  for (i in seq_len(n_trials)) {
    expect_equal(
      scen1_out$pooled[[i]],
      naive_list[[i]],
      tolerance = 1e-12
    )
  }
})

# Tests for posteriors2Quantiles ----------------------------------------------

# ------------------------------------------------------------------
# Test: posteriors2Quantiles computes requested quantiles, mean and SD
# Input:
#   - Normal samples theta ~ N(2, 3^2); 5000 draws.
# Behaviour:
#   - Returns a matrix with rows: quantiles, Mean, SD.
# Expectations:
#   - Row/column names as specified; values close to empirical summaries.
# Why:
#   - Core summarisation helper for posterior draws.
# ------------------------------------------------------------------
test_that("posteriors2Quantiles: computes quantiles, mean, and sd for each column", {
  set.seed(123)
  
  theta <- rnorm(5000, mean = 2, sd = 3)
  post  <- cbind(theta = theta)
  
  quantiles <- c(0.25, 0.5, 0.75)
  
  out <- posteriors2Quantiles(
    quantiles  = quantiles,
    posteriors = post
  )
  
  expect_true(is.matrix(out))
  expect_identical(colnames(out), "theta")
  expect_identical(rownames(out),
                   c("25%", "50%", "75%", "Mean", "SD"))
  
  exp_q <- stats::quantile(theta, probs = quantiles)
  expect_equal(out[c("25%", "50%", "75%"), "theta"],
               exp_q,
               tolerance = 1e-12)
  
  expect_equal(out["Mean", "theta"], mean(theta), tolerance = 1e-12)
  expect_equal(out["SD",   "theta"], stats::sd(theta), tolerance = 1e-12)
})


# Tests for performAnalyses ----------------------------------------------------

# Shared scenario_list for performAnalyses tests
scenario_list_pa <- make_scenario_list_pa()

prior_parameters_list <- getPriorParameters(
  method_names = c("pooled", "exnex"),
  target_rates = c(0.5, 0.5, 0.5)
)

# ------------------------------------------------------------------
# Test: performAnalyses returns a well-formed analysis_list
# Input:
#   - scenario_list_pa; methods pooled + exnex; supplied priors.
# Behaviour:
#   - One analysis entry per scenario with quantiles_list, scenario_data,
#     and analysis_parameters.
# Expectations:
#   - S3 class "analysis_list"; names "scenario_k"; scenario_data unchanged.
# Why:
#   - High-level sanity check of analysis pipeline.
# ------------------------------------------------------------------
test_that("performAnalyses returns a well-formed analysis_list", {
  res <- performAnalyses(
    scenario_list         = scenario_list_pa,
    target_rates          = c(0.5, 0.5, 0.5),
    method_names          = c("pooled", "exnex"),
    prior_parameters_list = prior_parameters_list,
    n_mcmc_iterations     = 20,
    verbose               = FALSE
  )
  
  expect_s3_class(res, "analysis_list")
  expect_equal(length(res), length(scenario_list_pa))
  expect_identical(names(res), paste0("scenario_", sapply(scenario_list_pa, `[[`, "scenario_number")))
  
  scen1 <- res[[1]]
  expect_true(is.list(scen1))
  expect_true(all(c("quantiles_list", "scenario_data", "analysis_parameters") %in% names(scen1)))
  
  expect_identical(scen1$scenario_data, scenario_list_pa[[1]])
})

# ------------------------------------------------------------------
# Test: performAnalyses sorts method_names and stores in analysis_parameters
# Input:
#   - method_names in unsorted order: stratified, berry, pooled.
# Behaviour:
#   - Internally sorts method_names; quantiles_list and analysis_parameters
#     follow sorted order.
# Expectations:
#   - scen1$analysis_parameters$method_names sorted; quantiles_list names match.
# Why:
#   - Ensures deterministic method ordering across analyses.
# ------------------------------------------------------------------
test_that("performAnalyses sorts method_names and stores them in analysis_parameters", {
  res <- performAnalyses(
    scenario_list       = scenario_list_pa,
    target_rates        = c(0.5, 0.5, 0.5),
    method_names        = c("stratified", "berry", "pooled"),
    n_mcmc_iterations   = 20,
    verbose             = FALSE
  )
  
  scen1 <- res[[1]]
  
  expected_sorted <- sort(c("stratified", "berry", "pooled"))
  expect_identical(scen1$analysis_parameters$method_names, expected_sorted)
  
  expect_true(is.list(scen1$quantiles_list))
  expect_identical(names(scen1$quantiles_list), expected_sorted)
})

# ------------------------------------------------------------------
# Test: performAnalyses builds quantiles from defaults + evidence_levels
# Input:
#   - evidence_levels = c(0.1,0.2,0.3)
# Behaviour:
#   - analysis_parameters$quantiles is union of reversed evidence_levels
#     and default quantiles, sorted.
# Expectations:
#   - Numeric, sorted; contains all default and evidence-based quantiles.
# Why:
#   - Confirms evidence_levels are incorporated into posterior summaries.
# ------------------------------------------------------------------
test_that("performAnalyses constructs quantiles from defaults and evidence_levels", {
  ev_levels <- c(0.1, 0.2, 0.3)
  
  res <- performAnalyses(
    scenario_list       = scenario_list_pa,
    evidence_levels     = ev_levels,
    target_rates        = c(0.5, 0.5, 0.5),
    method_names        = c("pooled"),
    n_mcmc_iterations   = 20,
    verbose             = FALSE
  )
  
  q <- res[[1]]$analysis_parameters$quantiles
  expect_true(is.numeric(q))
  expect_false(is.unsorted(q))
  
  defaults <- c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975)
  expected_q <- sort(unique(round(1 - c(defaults, ev_levels), 9)))
  
  expect_true(all(expected_q %in% q))
})

# ------------------------------------------------------------------
# Test: performAnalyses auto-fills prior_parameters_list when omitted
# Input:
#   - method_names = c("berry","pooled"), prior_parameters_list = NULL.
# Behaviour:
#   - Internally calls getPriorParameters and stores result in analysis_parameters.
# Expectations:
#   - Non-null priors; list; contains both methods as names.
# Why:
#   - Convenience behaviour for users not providing priors explicitly.
# ------------------------------------------------------------------
test_that("performAnalyses fills prior_parameters_list when not supplied", {
  methods <- c("berry", "pooled")
  
  res <- performAnalyses(
    scenario_list         = scenario_list_pa,
    target_rates          = c(0.5, 0.5, 0.5),
    method_names          = methods,
    prior_parameters_list = NULL,
    n_mcmc_iterations     = 20,
    verbose               = FALSE
  )
  
  priors <- res[[1]]$analysis_parameters$prior_parameters_list
  
  expect_false(is.null(priors))
  expect_true(is.list(priors))
  expect_true(all(methods %in% names(priors)))
})

# ------------------------------------------------------------------
# Test: performAnalyses prints a progress message when verbose = TRUE
# Input:
#   - verbose = TRUE
# Behaviour:
#   - Progress message "Performing Analyses" printed to console.
# Expectations:
#   - expect_message detects that substring.
# Why:
#   - Confirms user-facing verbosity option works.
# ------------------------------------------------------------------
test_that("performAnalyses prints a progress message when verbose = TRUE", {
  expect_message(
    performAnalyses(
      scenario_list       = scenario_list_pa,
      target_rates        = c(0.5, 0.5, 0.5),
      method_names        = c("pooled"),
      n_mcmc_iterations   = 20,
      verbose             = TRUE
    ),
    "Performing Analyses"
  )
})


# Tests for loadAnalyses -------------------------------------------------------

# ------------------------------------------------------------------
# Test: loadAnalyses reads saved analyses and sets class/names
# Input:
#   - Two RDS files matching analysis_data_1_1 and analysis_data_2_2.
# Behaviour:
#   - Reads into an analysis_list, named scenario_1, scenario_2.
# Expectations:
#   - S3 class "analysis_list"; names c("scenario_1","scenario_2");
#     contents identical to original dummies.
# Why:
#   - Confirms correct file naming convention and object reconstruction.
# ------------------------------------------------------------------
test_that("loadAnalyses: loads saved analyses and sets class/names correctly", {
  tmpdir <- tempdir()
  
  scen_nums <- c(1, 2)
  anal_nums <- c(1, 2)
  
  dummy1 <- list(foo = 1)
  dummy2 <- list(bar = 2)
  
  saveRDS(dummy1, file = file.path(tmpdir, "analysis_data_1_1.rds"))
  saveRDS(dummy2, file = file.path(tmpdir, "analysis_data_2_2.rds"))
  
  loaded <- loadAnalyses(
    scenario_numbers = scen_nums,
    analysis_numbers = anal_nums,
    load_path        = tmpdir
  )
  
  expect_s3_class(loaded, "analysis_list")
  expect_identical(names(loaded), c("scenario_1", "scenario_2"))
  
  expect_identical(loaded[[1]], dummy1)
  expect_identical(loaded[[2]], dummy2)
})

# ------------------------------------------------------------------
# Test: loadAnalyses validates scenario_numbers as positive integers
# Input:
#   - scenario_numbers = "a" or c(0,1).
# Behaviour:
#   - Non-integer or non-positive scenario_numbers should fail.
# Expectations:
#   - Errors mentioning "scenario_numbers".
# Why:
#   - Enforces basic integrity on index arguments.
# ------------------------------------------------------------------
test_that("loadAnalyses: scenario_numbers must be positive integers", {
  tmpdir <- tempdir()
  
  expect_error(
    loadAnalyses(scenario_numbers = "a", load_path = tmpdir),
    "scenario_numbers"
  )
  
  expect_error(
    loadAnalyses(scenario_numbers = c(0, 1), load_path = tmpdir),
    "scenario_numbers"
  )
})

# ------------------------------------------------------------------
# Test: loadAnalyses validates analysis_numbers length and positivity
# Input:
#   - analysis_numbers length mismatched to scenario_numbers or non-positive.
# Behaviour:
#   - Should error with appropriate message.
# Expectations:
#   - Errors mentioning "analysis_numbers".
# Why:
#   - Prevents ambiguous mapping between scenarios and analyses.
# ------------------------------------------------------------------
test_that("loadAnalyses: analysis_numbers must be positive integers of same length", {
  tmpdir <- tempdir()
  
  expect_error(
    loadAnalyses(
      scenario_numbers = c(1, 2),
      analysis_numbers = c(1),
      load_path        = tmpdir
    ),
    "analysis_numbers"
  )
  
  expect_error(
    loadAnalyses(
      scenario_numbers = c(1, 2),
      analysis_numbers = c(1, -1),
      load_path        = tmpdir
    ),
    "analysis_numbers"
  )
})

# ------------------------------------------------------------------
# Test: loadAnalyses requires load_path be a single character string
# Input:
#   - load_path = 123 or c("a","b").
# Behaviour:
#   - Type and length checks must fail.
# Expectations:
#   - Errors mentioning "load_path".
# Why:
#   - Guards against invalid path inputs.
# ------------------------------------------------------------------
test_that("loadAnalyses: load_path must be a single character string", {
  expect_error(
    loadAnalyses(scenario_numbers = 1, load_path = 123),
    "load_path"
  )
  expect_error(
    loadAnalyses(scenario_numbers = 1, load_path = c("a", "b")),
    "load_path"
  )
})


# Tests for printAnalyses ------------------------------------------------------

# ------------------------------------------------------------------
# Test: print.analysis_list prints header, scenario blocks, labels, and estimates
# Input:
#   - analyses from pooled method on 2 scenarios; estimates via getEstimates().
# Behaviour:
#   - print() shows header, per-scenario listing, "Pooled" label,
#     and numeric summary values like mean and SD.
# Expectations:
#   - Output contains correct header; scenario lines; method label;
#     and numeric values matching rounded estimates.
# Why:
#   - Verifies user-facing summary of analysis_list.
# ------------------------------------------------------------------
test_that("print.analysis_list: prints header, scenario blocks, method label, and numeric estimates", {
  set.seed(123)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(10, 20),
                               c(10, 20)),
    response_rates_list = list(c(0.3, 0.6),
                               c(0.4, 0.7)),
    n_trials            = 50
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    method_names       = "pooled",
    target_rates       = c(0.5, 0.5),
    n_mcmc_iterations  = 20,
    verbose            = FALSE
  )
  
  est_list <- getEstimates(analyses)
  
  expect_true(is.list(est_list[[1]]))
  est_mat1 <- est_list[[1]][[1]]
  
  mean_p1 <- round(est_mat1["p_1", "Mean"], digits = 2)
  sd_p1   <- round(est_mat1["p_1", "SD"],   digits = 2)
  
  out <- capture.output(print(analyses))
  flat_out <- paste(out, collapse = " ")
  
  expect_true(
    any(grepl("analysis_list of 2 scenarios with 1 method", out))
  )
  
  expect_true(
    any(grepl("^  - scenario_1", out))
  )
  expect_true(
    any(grepl("^  - scenario_2", out))
  )
  
  expect_true(
    any(grepl("Pooled", out))
  )
  
  mean_str <- sprintf("%.2f", mean_p1)
  sd_str   <- sprintf("%.2f", sd_p1)
  
  expect_true(
    grepl(mean_str, flat_out)
  )
  expect_true(
    grepl(sd_str, flat_out)
  )
  
  expect_true(
    any(grepl("MCMC iterationns per BHM method", out))
  )
  expect_true(
    any(grepl("Available evidence levels:", out))
  )
})

# ------------------------------------------------------------------
# Test: print.analysis_list shows distinct scenario-specific estimates
# Input:
#   - Two scenarios with different response rates, pooled method.
# Behaviour:
#   - Printed output must include different means for p_1 per scenario.
# Expectations:
#   - Output includes both rounded means; these means differ numerically.
# Why:
#   - Ensures print summarises scenario-specific posterior differences.
# ------------------------------------------------------------------
test_that("print.analysis_list: multiple scenarios print distinct scenario-specific estimates", {
  set.seed(456)
  
  n_subj <- c(10, 20)
  rr1    <- c(0.3, 0.6)
  rr2    <- c(0.7, 0.2)
  
  scen_multi <- simulateScenarios(
    n_subjects_list     = list(n_subj, n_subj),
    response_rates_list = list(rr1, rr2),
    n_trials            = 30
  )
  
  analyses_multi <- performAnalyses(
    scenario_list      = scen_multi,
    method_names       = "pooled",
    target_rates       = c(0.5, 0.5),
    n_mcmc_iterations  = 20,
    verbose            = FALSE
  )
  
  est_multi <- getEstimates(analyses_multi)
  
  expect_true(is.list(est_multi[[1]]))
  m1 <- est_multi[[1]][[1]]
  m2 <- est_multi[[1]][[2]]
  
  mean1 <- round(m1["p_1", "Mean"], 2)
  mean2 <- round(m2["p_1", "Mean"], 2)
  
  out <- capture.output(print(analyses_multi))
  flat_out <- paste(out, collapse = " ")
  
  expect_true(
    any(grepl("analysis_list of 2 scenarios with 1 method", out))
  )
  
  expect_true(any(grepl("^  - scenario_1", out)))
  expect_true(any(grepl("^  - scenario_2", out)))
  
  mean1_str <- sprintf("%.2f", mean1)
  mean2_str <- sprintf("%.2f", mean2)
  
  expect_true(
    grepl(mean1_str, flat_out)
  )
  expect_true(
    grepl(mean2_str, flat_out)
  )
  
  expect_false(
    isTRUE(all.equal(mean1, mean2))
  )
})

# ------------------------------------------------------------------
# Test: digits argument controls numeric precision in printed output
# Input:
#   - digits = 2 and digits = 4; same analyses object.
# Behaviour:
#   - Numeric summaries printed with specified precision.
# Expectations:
#   - Output with digits=2 contains 2-decimal format; digits=4 contains 4-decimal format.
# Why:
#   - Confirms users can tune displayed precision.
# ------------------------------------------------------------------
test_that("print.analysis_list: digits argument controls printed numeric precision", {
  set.seed(789)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(10, 20),
                               c(10, 20)),
    response_rates_list = list(c(0.33, 0.66),
                               c(0.33, 0.66)),
    n_trials            = 40
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    method_names       = "pooled",
    target_rates       = c(0.5, 0.5),
    n_mcmc_iterations  = 20,
    verbose            = FALSE
  )
  
  est_list <- getEstimates(analyses)
  est_obj  <- est_list[[1]]
  expect_true(is.list(est_obj))
  est_mat1 <- est_obj[[1]]
  
  mean_p1 <- est_mat1["p_1", "Mean"]
  
  out_2 <- capture.output(print(analyses, digits = 2))
  out_4 <- capture.output(print(analyses, digits = 4))
  
  flat_2 <- paste(out_2, collapse = " ")
  flat_4 <- paste(out_4, collapse = " ")
  
  fmt2 <- sprintf("%.2f", round(mean_p1, 2))
  fmt4 <- sprintf("%.4f", round(mean_p1, 4))
  
  expect_true(
    grepl(fmt2, flat_2)
  )
  
  expect_true(
    grepl(fmt4, flat_4)
  )
})


# Tests for saveAnalyses -------------------------------------------------------

# ------------------------------------------------------------------
# Test: saveAnalyses writes analysis_list to disk and loadAnalyses reads it back
# Input:
#   - analyses from pooled method, multiple scenarios; temp directory.
# Behaviour:
#   - saveAnalyses returns path and scenario/analysis numbers;
#     files exist; loadAnalyses reproduces original object.
# Expectations:
#   - info list correctly populated; RDS files present;
#     loaded analysis_list identical in length/names to original.
# Why:
#   - Round-trip test for persistence layer.
# ------------------------------------------------------------------
test_that("saveAnalyses: saves analysis_list to disk and loadAnalyses can read it back", {
  set.seed(101)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(10, 20),
                               c(15, 25)),
    response_rates_list = list(c(0.3, 0.6),
                               c(0.5, 0.8)),
    n_trials            = 20
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    method_names       = "pooled",
    target_rates       = c(0.5, 0.5),
    n_mcmc_iterations  = 20,
    verbose            = FALSE
  )
  
  tmp_dir <- tempdir()
  
  info <- saveAnalyses(analyses, save_path = tmp_dir)
  
  expect_true(is.list(info))
  expect_equal(length(info$scenario_numbers), length(analyses))
  expect_equal(length(info$analysis_numbers), length(analyses))
  expect_identical(info$path, tmp_dir)
  
  internal_scen <- vapply(
    analyses,
    function(x) x$scenario_data$scenario_number,
    FUN.VALUE = numeric(1)
  )
  expect_equal(info$scenario_numbers, internal_scen)
  
  for (i in seq_along(info$scenario_numbers)) {
    f <- file.path(
      tmp_dir,
      paste0("analysis_data_",
             info$scenario_numbers[i], "_",
             info$analysis_numbers[i], ".rds")
    )
    expect_true(file.exists(f))
  }
  
  loaded <- loadAnalyses(
    scenario_numbers = info$scenario_numbers,
    analysis_numbers = info$analysis_numbers,
    load_path        = info$path
  )
  
  expect_s3_class(loaded, "analysis_list")
  expect_identical(length(loaded), length(analyses))
  expect_identical(names(loaded), names(analyses))
})

# ------------------------------------------------------------------
# Test: saveAnalyses requires an object of class 'analysis_list'
# Input:
#   - list(a = 1) instead of analysis_list.
# Behaviour:
#   - Class check must fail with informative message.
# Expectations:
#   - Error mentioning "analysis_list".
# Why:
#   - Prevents misuse with arbitrary lists.
# ------------------------------------------------------------------
test_that("saveAnalyses: non-analysis_list input triggers a class error", {
  expect_error(
    saveAnalyses(list(a = 1)),
    "analysis_list",
    ignore.case = TRUE
  )
})

# ------------------------------------------------------------------
# Test: saveAnalyses requires save_path to be a length-1 character vector
# Input:
#   - save_path = c("a","b").
# Behaviour:
#   - Argument validation should fail.
# Expectations:
#   - Error mentioning "save_path".
# Why:
#   - Enforces a single destination directory.
# ------------------------------------------------------------------
test_that("saveAnalyses: save_path must be a character vector of length 1", {
  set.seed(202)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(10, 20)),
    response_rates_list = list(c(0.3, 0.6)),
    n_trials            = 10
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    method_names       = "pooled",
    target_rates       = c(0.5, 0.5),
    n_mcmc_iterations  = 20,
    verbose            = FALSE
  )
  
  expect_error(
    saveAnalyses(analyses, save_path = c("a", "b")),
    "save_path",
    ignore.case = TRUE
  )
})