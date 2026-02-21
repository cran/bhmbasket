# Setup helpers for tests

setup_applicable_prev_trials <- function() {
  set.seed(123)
  
  scen_initial <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(c(0.4, 0.6, 0.8)),
    n_trials            = 10
  )
  
  analysis_initial <- performAnalyses(
    scenario_list       = scen_initial,
    target_rates        = c(0.5, 0.5, 0.5),
    method_names        = "pooled",
    n_mcmc_iterations   = 50,
    verbose             = FALSE
  )
  
  go_initial <- getGoDecisions(
    analyses_list   = analysis_initial,
    cohort_names    = c("p_1", "p_2", "p_3"),
    evidence_levels = c(0.5, 0.5, 0.5),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  scen_next <- continueRecruitment(
    n_subjects_add_list = list(c(5, 5, 5)),
    decisions_list      = go_initial,
    method_name         = "pooled"
  )
  
  analysis_with_diff <- performAnalyses(
    scenario_list       = scen_initial,
    target_rates        = c(0.5, 0.5, 0.5),
    method_names        = "pooled",
    calc_differences    = c(3, 2),
    n_mcmc_iterations   = 50,
    verbose             = FALSE
  )
  
  go_with_diff <- getGoDecisions(
    analyses_list   = analysis_with_diff,
    cohort_names    = c("p_1", "p_2", "p_3"),
    evidence_levels = c(0.5, 0.5, 0.5),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  scen_next_diff <- continueRecruitment(
    n_subjects_add_list = list(c(5, 5, 5)),
    decisions_list      = go_with_diff,
    method_name         = "pooled"
  )
  
  method_names_prev <- analysis_initial[[1]]$analysis_parameters$method_names
  quantiles_prev    <- analysis_initial[[1]]$analysis_parameters$quantiles
  n_coh_prev        <- ncol(scen_next[[1]]$n_subjects)
  
  list(
    scen_initial       = scen_initial,
    analysis_initial   = analysis_initial,
    go_initial         = go_initial,
    scen_next          = scen_next,
    analysis_with_diff = analysis_with_diff,
    go_with_diff       = go_with_diff,
    scen_next_diff     = scen_next_diff,
    method_names_prev  = method_names_prev,
    quantiles_prev     = quantiles_prev,
    n_coh_prev         = n_coh_prev
  )
}

make_scenario_list_pa <- function() {
  set.seed(123)
  simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(c(0.3, 0.5, 0.7)),
    n_trials            = 5
  )
}

setup_two_scenario_pooled <- function(
    n_subj,
    rr1,
    rr2,
    n_trials = 30,
    n_mcmc_iterations = 20,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  scen <- simulateScenarios(
    n_subjects_list     = list(n_subj, n_subj),
    response_rates_list = list(rr1, rr2),
    n_trials            = n_trials
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    method_names       = "pooled",
    target_rates       = rep(0.5, length(n_subj)),
    n_mcmc_iterations  = n_mcmc_iterations,
    verbose            = FALSE
  )
  
  list(
    scen     = scen,
    analyses = analyses
  )
}