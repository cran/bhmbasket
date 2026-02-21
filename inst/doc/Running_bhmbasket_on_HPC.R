## ----knitr_options, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  eval     = FALSE,  # en- / disables R code evaluation globally
  cache    = FALSE,  # en- / disables R code caching globally
  collapse = TRUE,
  comment  = "#>"
)

## ----setup--------------------------------------------------------------------
# library(bhmbasket)
# library(doFuture)
# library(future.batchtools)
# 
# rng_seed <- 5440
# set.seed(rng_seed)

## ----SLURM_Setup--------------------------------------------------------------
# ## Adapt the SLURM template to requirements
# job_time  <- 1   # time for job in hours
# n_workers <- 24  # number of worker nodes
# n_cpus    <- 16  # number of cpus per worker node
# gb_memory <- 2   # memory [GB] per cpu
# 
# slurm <- tweak(batchtools_slurm,
#            template  = system.file('templates/slurm-simple.tmpl',
#                                    package = 'batchtools'),
#            workers   = n_workers,
#            resources = list(
#              walltime  = 60 * 60 * job_time,
#              ncpus     = n_cpus,
#              memory    = 1000 * gb_memory))
# 
# ## Register the parallel backend
# registerDoFuture()
# 
# ## Specify how the futures should be resolved
# plan(list(slurm, multisession))

## -----------------------------------------------------------------------------
# scenarios_list <- simulateScenarios(
#     n_subjects_list     = list(c(10, 20, 30)),
#     response_rates_list = list(c(0.1, 0.2, 3)),
#     n_trials            = 10)
# 
#   analyses_list <- performAnalyses(
#     scenario_list       = scenarios_list,
#     target_rates        = c(0.1, 0.1, 0.1),
#     calc_differences    = matrix(c(3, 2, 2, 1), ncol = 2),
#     n_mcmc_iterations   = 100)
# 
#   getEstimates(analyses_list)

