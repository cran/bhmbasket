## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  eval     = TRUE,  # en- / disables R code evaluation globally
  cache    = FALSE,  # en- / disables R code caching globally
  collapse = TRUE,
  comment  = "#>"
)

## ----setup--------------------------------------------------------------------
library(bhmbasket)

rng_seed <- 123
set.seed(rng_seed)

## -----------------------------------------------------------------------------
chugh_trial <- createTrial(
  n_subjects   = c(15, 13, 12, 28, 29, 29, 26, 5, 2, 20),
  n_responders = c(2, 0, 1, 6, 7, 3, 5, 1, 0, 3))

## ----chugh--------------------------------------------------------------------
chugh_prior_parameters <- setPriorParametersExNex(
  mu_mean   = c(-1.735, 0.847),
  mu_sd     = c(0.146, 0.266)^-0.5,
  tau_scale = 1,
  mu_j      = rep(-1.734, 10),
  tau_j     = rep(0.128^-0.5, 10),
  w_j       = c(0.5, 0, 0.5))

if (file.exists("chugh_analysis.rds")) {
  
  chugh_analysis <- readRDS("chugh_analysis.rds")
  
} else {
  
  chugh_analysis <- performAnalyses(
  scenario_list         = chugh_trial,
  method_names          = "exnex",
  prior_parameters_list = chugh_prior_parameters,
  seed                  = rng_seed,
  n_mcmc_iterations     = 5e4,
  n_cores               = 2L,
  verbose               = FALSE)
  
  saveRDS(chugh_analysis, "chugh_analysis.rds")
  
}

## -----------------------------------------------------------------------------
getEstimates(analyses_list = chugh_analysis)

## ----simulateScenarios--------------------------------------------------------
scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(20, 20, 10, 10),
                             c(20, 20, 10, 10),
                             c(20, 20, 10, 10)),
  response_rates_list = list(c(0.1, 0.1, 0.1, 0.1),
                             c(0.1, 0.1, 0.3, 0.3),
                             c(0.1, 0.1, 0.1, 0.5)),
  scenario_numbers    = c(1, 3, 4),
  n_trials            = 1e4)

## ----performAnalyses----------------------------------------------------------
prior_parameters <- setPriorParametersExNex(
  mu_mean   = c(logit(0.1), logit(0.3)),
  mu_sd     = c(3.18, 1.94),
  tau_scale = 1,
  mu_j      = rep(logit(0.2), 4),
  tau_j     = rep(2.5, 4),
  w_j       = c(0.25, 0.25, 0.5))

if (file.exists("analyses_list.rds")) {
  
  analyses_list <- readRDS("analyses_list.rds")
  
} else {
  
  analyses_list <- performAnalyses(
    scenario_list         = scenarios_list,
    method_names          = "exnex",
    prior_parameters_list = prior_parameters,
    seed                  = rng_seed,
    n_mcmc_iterations     = 5e4,
    n_cores               = 2L)
  
  saveRDS(analyses_list, "analyses_list.rds")
  
}

## ----include = FALSE----------------------------------------------------------
## check for unnecessary data
if (any(grepl("mu_",
              colnames(analyses_list$scenario_1$quantiles_list$exnex[[1]])))) {
  
  ## remove unnecessary data
  
  ## ...in quantiles_list
  analyses_list$scenario_1$quantiles_list$exnex <-   lapply(analyses_list$scenario_1$quantiles_list$exnex, function   (x) {
    x[-c(2, 6), -c(1:6, 11, 12)]
  })
  analyses_list$scenario_3$quantiles_list$exnex <-   lapply(analyses_list$scenario_3$quantiles_list$exnex, function   (x) {
    x[-c(2, 6), -c(1:6, 11, 12)]
  })
  analyses_list$scenario_4$quantiles_list$exnex <-   lapply(analyses_list$scenario_4$quantiles_list$exnex, function   (x) {
    x[-c(2, 6), -c(1:6, 11, 12)]
  })
  
  ## ... on quantiles
  analyses_list$scenario_1$analysis_parameters$quantiles <-   c(0.025, 0.100, 0.200, 0.500, 0.975)
  analyses_list$scenario_3$analysis_parameters$quantiles <-   c(0.025, 0.100, 0.200, 0.500, 0.975)
  analyses_list$scenario_4$analysis_parameters$quantiles <-   c(0.025, 0.100, 0.200, 0.500, 0.975)
  
  ## ... on subjects and responders
  analyses_list$scenario_1$scenario_data$n_subjects <- NULL
  analyses_list$scenario_3$scenario_data$n_subjects <- NULL
  analyses_list$scenario_4$scenario_data$n_subjects <- NULL
  analyses_list$scenario_1$scenario_data$n_responders <- NULL
  analyses_list$scenario_3$scenario_data$n_responders <- NULL
  analyses_list$scenario_4$scenario_data$n_responders <- NULL
  
  analyses_list$scenario_1$quantiles_list$exnex[[1]]
    
  saveRDS(analyses_list, "analyses_list.rds", compress = "xz")
}

## ----getBiasMSE---------------------------------------------------------------
estimates <- getEstimates(
    analyses_list   = analyses_list,
    point_estimator = "mean")

scaleRoundList(
  list         = estimates,
  scale_param  = 100,
  round_digits = 2)

## ----eval = FALSE-------------------------------------------------------------
# cohort_names    = "p_1"
# boundary_rules  = quote(c(x[1] > 0.1))
# evidence_levels = 0.9

## ----eval = FALSE-------------------------------------------------------------
# cohort_names    = c("p_1", "p_1")
# boundary_rules  = quote(c(x[1] > 0.1 & x[2] > 0.2))
# evidence_levels = c(0.9, "mean")

## ----getGoDecisions-----------------------------------------------------------
decisions_list <- getGoDecisions(
  analyses_list   = analyses_list,
  cohort_names    = c("p_1", "p_1",
                      "p_2", "p_2",
                      "p_3", "p_3",
                      "p_4", "p_4"),
  evidence_levels = c(0.9, "mean",
                      0.9, "mean",
                      0.8, "mean",
                      0.8, "mean"),
  boundary_rules  = quote(c(x[1] > 0.1 & x[2] > 0.2,
                            x[3] > 0.1 & x[4] > 0.2,
                            x[5] > 0.1 & x[6] > 0.2,
                            x[7] > 0.1 & x[8] > 0.2)))

go_probabilities <- getGoProbabilities(decisions_list)

scaleRoundList(
  list         = go_probabilities,
  scale_param  = 100,
  round_digits = 0)

