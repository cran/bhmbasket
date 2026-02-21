applicablePreviousTrials <- function(
    
  scenario_list,
  method_names,
  quantiles,
  n_cohorts,
  calc_differences
  
) {
  
  ## analyze only unique trials that have not been previously analyzed,
  ## i.e. trials that were updated with continueRecruitment(),
  ## i.e. trials that have an overall go decision from a previous decision rule
  ## i.e. trials that have quantiles stored for each cohort that is to be analysed
  applicable_previous_trials <-
    ## check that in each scenario the same analysis methods were analyzed previously
    all(sapply(seq_along(scenario_list), function (i) {
      isTRUE(all.equal(names(scenario_list[[i]]$previous_analyses$post_quantiles),
                       names(scenario_list[[1]]$previous_analyses$post_quantiles)))
    })) &
    ## check that the current analysis method names match the method names of the previous analyses
    isTRUE(all.equal(names(scenario_list[[1]]$previous_analyses$post_quantiles), method_names)) &
    ## check that the stored quantiles are the same across all scenarios
    all(sapply(seq_along(scenario_list), function (i) {
      isTRUE(all.equal(rownames(scenario_list[[i]]$previous_analyses$post_quantiles[[1]][[1]]),
                       rownames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]])))
    })) &
    ## check that the new quantiles are within the stored quantiles
    all(paste0(as.character(quantiles * 100), "%") %in%
          rownames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]])) &
    ## check that there are stored quantiles for each cohort that is to be analysed
    all(paste0("p_", seq_len(n_cohorts)) %in%
          colnames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]]))
  
  ## check that all differences have been previously calculated
  if (!is.null(calc_differences)) {
    
    applicable_previous_trials <- applicable_previous_trials &
      all(apply(calc_differences, 1, function (x) {
        paste0("p_diff_", paste0(as.character(x), collapse = ""))
      }) %in% colnames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]]))
    
  }
  
  return (applicable_previous_trials)
  
}

calcDiffsMCMC <- function (
    
  posterior_samples,
  calc_differences
  
) {
  
  org_names <- colnames(posterior_samples)
  
  diffs <- apply(calc_differences, 1, function (x) {
    
    matrix(posterior_samples[, grepl(x[1], org_names) & grepl("p", org_names)] -
             posterior_samples[, grepl(x[2], org_names) & grepl("p", org_names)],
           ncol = 1)
    
  })
  
  diff_names <- apply(calc_differences, 1, function (x) {
    
    paste0("p_diff_", paste0(as.character(x), collapse = ""))
    
  })
  
  colnames(diffs) <- diff_names
  
  posterior_samples <- cbind(posterior_samples, diffs)
  
  return (posterior_samples)
  
}

getModelFile <- function (method_name) {
  
  if (method_name == "berry") {
    
    model_file <- "berry.txt"
    
  } else if (method_name == "berry_mix") {
    
    model_file <- "berry_mix.txt"
    
  } else if (method_name == "exnex") {
    
    model_file <- "exnex.txt"
    
  } else if (method_name == "exnex_adj") {
    
    model_file <- "exnex_adj.txt"
    
  } else {
    
    stop ("method_name must be one of berry, berry_mix, exnex, exnex_adj")
    
  }
  
  model_file <- system.file(package = "bhmbasket", "jags_models", model_file, mustWork = TRUE)
  
  return (model_file)
  
}

getPosteriors <- function (
    
  j_parameters,
  j_model_file,
  j_data,
  
  n_mcmc_iterations
  
) {
  
  ## Adaption and burn-in not included
  posterior_samples <- performJags(
    data               = j_data,
    parameters_to_save = j_parameters,
    model_file         = j_model_file, 
    n_iter             = n_mcmc_iterations)
  
  ## replace squarebrackets provided by rjags with workable characters
  colnames(posterior_samples) <- gsub("\\[", "_", colnames(posterior_samples))
  colnames(posterior_samples) <- gsub("\\]", "", colnames(posterior_samples))
  
  weights_indices <- grepl("exch", colnames(posterior_samples))
  if (any(weights_indices)) {
    
    superfluous_weights <- !grepl(",1", colnames(posterior_samples))
    
    colnames(posterior_samples)[weights_indices] <-
      paste0("w_", seq_along(j_data$n))
    
    posterior_samples <- posterior_samples[, !(weights_indices & superfluous_weights)]
    
  }
  
  return (posterior_samples)
  
}

getPostQuantiles <- function (
    
  ## The method to be applied to the likelihood and the quantiles of the posterior
  method_name,
  quantiles,
  
  ## Scenario data
  scenario_data,
  
  ## Differences between cohorts
  calc_differences  = NULL,
  
  ## JAGS parameters
  j_parameters,
  j_model_file,
  j_data,
  
  ## MCMC Parameters
  n_mcmc_iterations = 1e4,
  
  ## Where to save one of the posterior response rates approximations provided by JAGS
  save_path         = NULL,
  save_trial        = NULL
  
) {
  
  if (is.null(dim(scenario_data$n_responders))) {
    
    scenario_data$n_responders <- t(convertVector2Matrix(scenario_data$n_responders))
    scenario_data$n_subjects   <- t(convertVector2Matrix(scenario_data$n_subjects))
    
  }
  
  n_analyses <- nrow(scenario_data$n_responders)
  
  ## Create random index for saving one of the posterior response rates
  if (is.null(save_trial) && !is.null(save_path)) {
    # set.seed(seed)
    save_trial <- sample(seq_len(n_analyses), size = 1)
  }
  
  ## Run parallel loops
  ## prepare foreach loop over 
  
  exported_stuff <- c(
    "posteriors2Quantiles", "getPosteriors", "getPostQuantilesOfTrial",
    "qbetaDiff", "chunkVector")
  
  ## prepare chunking
  chunks_outer <- chunkVector(seq_len(n_analyses), foreach::getDoParWorkers())
  
  "%dorng%" <- doRNG::"%dorng%"
  "%dopar%" <- foreach::"%dopar%"
  posterior_quantiles_list <- suppressMessages(
    foreach::foreach(
      k = chunks_outer,
      .combine  = c,
      .verbose  = FALSE,
      .packages = c("rjags"),
      .export   = exported_stuff) %dorng% {
        
        chunks_inner <- chunkVector(k, foreach::getDoParWorkers())
        
        foreach::foreach(i = chunks_inner, .combine = c) %dorng% {
          
          lapply(i, function (j) {
            
            ## Calculate the posterior quantiles for the kth unique trial outcome
            getPostQuantilesOfTrial(
              n_responders      = as.numeric(scenario_data$n_responders[j, ]),
              n_subjects        = as.numeric(scenario_data$n_subjects[j, ]),
              j_data            = j_data,
              j_parameters      = j_parameters,
              j_model_file      = j_model_file,
              method_name       = method_name,
              quantiles         = quantiles,
              calc_differences  = calc_differences,
              n_mcmc_iterations = n_mcmc_iterations,
              save_path         = save_path,
              save_trial        = save_trial)
            
          })
          
        }
        
      })
  
  return (posterior_quantiles_list)
  
}

getPostQuantilesOfTrial <- function (
    
  n_responders,
  n_subjects,
  
  j_data,
  j_parameters,
  j_model_file,
  method_name,
  quantiles,
  calc_differences,
  n_mcmc_iterations,
  
  save_path,
  save_trial
  
) {
  
  j_data$r <- n_responders
  j_data$n <- n_subjects
  
  if (method_name == "stratified") {
    
    posterior_quantiles <- getPostQuantilesStratified(
      j_data            = j_data,
      quantiles         = quantiles,
      calc_differences  = calc_differences,
      n_mcmc_iterations = n_mcmc_iterations)
    
  } else if (method_name == "pooled") {
    
    posterior_quantiles <- getPostQuantilesPooled(
      j_data           = j_data,
      quantiles        = quantiles,
      calc_differences = calc_differences)
    
  } else {
    
    ## Get posterior response rates per indication
    posterior_samples <- getPosteriors(
      j_parameters      = j_parameters,
      j_model_file      = j_model_file,
      j_data            = j_data,
      n_mcmc_iterations = n_mcmc_iterations)
    
    ## Calculate differences between response rates of cohorts
    if (!is.null(calc_differences)) {
      
      posterior_samples <- calcDiffsMCMC(
        posterior_samples = posterior_samples,
        calc_differences  = calc_differences)
      
    }
    
    ## Save posterior response rates per indication for one randomly selected simulation,
    ## due to time and storage space constraints only one simulation
    if (!is.null(save_path)) {
      if (k == save_trial) {
        saveRDS(posterior_samples,
                file = file.path(save_path, paste0("posterior_samples_",
                                                   k, "_", method_name, "_rds")))
      }
    }
    
    ## Calculate the required quantiles for the decision rules
    posterior_quantiles <- posteriors2Quantiles(
      quantiles  = quantiles,
      posteriors = posterior_samples)
    
  }
  
  return (posterior_quantiles)
  
}

getPostQuantilesPooled <- function(
    
  j_data,
  quantiles,
  calc_differences
  
) {
  
  shape_1 <- j_data$a + sum(j_data$r)
  shape_2 <- j_data$b + sum(j_data$n) - sum(j_data$r)
  
  posterior_quantiles <- stats::qbeta(quantiles, shape1 = shape_1, shape2 = shape_2)
  
  posterior_quantiles <- matrix(posterior_quantiles,
                                ncol = length(j_data$r), nrow = length(quantiles))
  
  colnames(posterior_quantiles) <- paste0("p_", seq_along(j_data$r))
  rownames(posterior_quantiles) <- paste0(quantiles * 100, "%")
  
  posterior_mean      <- shape_1 / (shape_1 + shape_2)
  posterior_sd        <- ((shape_1 * shape_2) / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1)))^0.5
  posterior_quantiles <- rbind(posterior_quantiles,
                               Mean = posterior_mean,
                               SD   = posterior_sd)
  
  if (!is.null(calc_differences)) {
    
    diff_quantiles <- apply(calc_differences, 1, function (x) {
      
      matrix(rep(0, nrow(posterior_quantiles)), ncol = 1)
      
    })
    
    diff_names <- apply(calc_differences, 1, function (x) {
      
      paste0("p_diff_", paste0(as.character(x), collapse = ""))
      
    })
    
    colnames(diff_quantiles) <- diff_names
    
    posterior_quantiles <- cbind(posterior_quantiles, diff_quantiles)
    
  }
  
  return (posterior_quantiles)
  
}

getPostQuantilesStratified <- function(
    
  j_data,
  quantiles,
  calc_differences,
  n_mcmc_iterations
  
) {
  
  shape_1 <- j_data$a_j + j_data$r
  shape_2 <- j_data$b_j + j_data$n - j_data$r
  
  posterior_quantiles <- t(sapply(quantiles, function (x)
    stats::qbeta(x, shape1 = shape_1, shape2 = shape_2)))
  
  if (nrow(posterior_quantiles) == 1) {
    
    posterior_quantiles <- t(posterior_quantiles)
    
  }
  
  colnames(posterior_quantiles) <- paste0("p_", seq_along(j_data$a_j))
  rownames(posterior_quantiles) <- paste0(quantiles * 100, "%")
  
  posterior_mean      <- shape_1 / (shape_1 + shape_2)
  posterior_sd        <- ((shape_1 * shape_2) / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1)))^0.5
  posterior_quantiles <- rbind(posterior_quantiles,
                               Mean = posterior_mean,
                               SD   = posterior_sd)
  
  if (!is.null(calc_differences)) {
    
    diff_quantiles <- apply(calc_differences, 1, function (x) {
      
      matrix(qbetaDiff(
        quantiles  = quantiles,
        x_1_shape1 = shape_1[x[1]],
        x_1_shape2 = shape_2[x[1]],
        x_2_shape1 = shape_1[x[2]],
        x_2_shape2 = shape_2[x[2]],
        n_mcmc     = n_mcmc_iterations), ncol = 1)
      
    })
    
    diff_names <- apply(calc_differences, 1, function (x) {
      
      paste0("p_diff_", paste0(as.character(x), collapse = ""))
      
    })
    
    colnames(diff_quantiles) <- diff_names
    
    posterior_quantiles <- cbind(posterior_quantiles, diff_quantiles)
    
  }
  
  return (posterior_quantiles)
  
}

getUniqueRows <- function (
    
  matrix
  
) {
  
  n_rows      <- nrow(matrix)
  n_cols      <- ncol(matrix)
  
  unique_rows <- stats::aggregate(id ~ .,
                                  data = cbind(id = seq_along(n_rows), matrix),
                                  FUN  = length)
  
  return (unique_rows[, seq_len(n_cols)])
  
}

getUniqueTrials <- function (
    
  scenario_list
  
) {
  
  all_scenarios_n_responders <- do.call(rbind, lapply(scenario_list, function (x) x$n_responders))
  all_scenarios_n_subjects   <- do.call(rbind, lapply(scenario_list, function (x) x$n_subjects))
  all_scenarios_overall_gos  <- do.call(rbind, lapply(scenario_list, function (x) 
    x$previous_analyses$go_decisions))[, 1]
  
  return (getUniqueRows(cbind(all_scenarios_n_responders,
                              all_scenarios_n_subjects,
                              go_flag = all_scenarios_overall_gos)))
  
}

is.analysis_list <- function (x) {
  
  if (missing(x)) stop ("Please provide an object for the argument 'x'")
  
  inherits(x, "analysis_list")
  
}

#' @title loadAnalyses
#' @md
#' @description This function loads an analysis performed with
#' \code{\link[bhmbasket]{performAnalyses}}
#' @param load_path A string providing a path where the scenarios are being stored,
#' Default: \code{\link[base]{tempfile}}
#' @param scenario_numbers A (vector of) positive integer(s) for the scenario number(s)
#' @param analysis_numbers A (vector of) positive integer(s) for the analysis number(s),
#' Default: `rep(1, length(scenario_numbers))`
#' @return Returns an object of class `analysis_list`
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{saveAnalyses}}
#'  \code{\link[base]{tempfile}}
#' @rdname loadAnalyses
#' @examples
#'   trial_data <- createTrial(
#'     n_subjects   = c(10, 20, 30),
#'     n_responders = c(1, 2, 3))
#'
#'   analysis_list <- performAnalyses(
#'     scenario_list      = trial_data,
#'     target_rates       = rep(0.5, 3),
#'     n_mcmc_iterations  = 100)
#'
#'   save_info     <- saveAnalyses(analysis_list)
#'   analysis_list <- loadAnalyses(scenario_numbers = save_info$scenario_numbers,
#'                                 analysis_numbers = save_info$analysis_numbers,
#'                                 load_path        = save_info$path)
#' @author Stephan Wojciekowski
#' @export
loadAnalyses <- function (
    
  scenario_numbers,
  analysis_numbers = rep(1, length(scenario_numbers)),
  load_path        = tempdir()
  
) {
  
  error_scenario_numbers <-
    "Please provide a vector of positive integers for the argument 'scenario_numbers'"
  error_analysis_numbers <-
    "Please provide a vector of positive integers for the argument 'analysis_numbers'"
  error_load_path <-
    "Please provide a string containing a path for the argument 'load_path'"
  error_len_match <-
    "'scenario_numbers' and 'analysis_numbers' must have equal length"
  
  checkmate::assertNumeric(
    scenario_numbers,
    any.missing = FALSE,
    .var.name   = error_scenario_numbers
  )
  checkmate::assertIntegerish(
    scenario_numbers,
    lower     = 1,
    .var.name = error_scenario_numbers
  )
  
  checkmate::assertIntegerish(
    analysis_numbers,
    lower       = 1,
    any.missing = FALSE,
    .var.name   = error_analysis_numbers
  )
  
  checkmate::assertTRUE(
    identical(length(scenario_numbers), length(analysis_numbers)),
    .var.name = error_len_match
  )
  
  checkmate::assertCharacter(
    load_path,
    len         = 1,
    any.missing = FALSE,
    .var.name   = error_load_path
  )
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  analyses_list <- vector(mode = "list", length = length(scenario_numbers))
  
  for (s in seq_along(scenario_numbers)) {
    
    analyses_list[[s]] <- readRDS(
      file.path(load_path,
                paste0("analysis_data_", scenario_numbers[s], "_", analysis_numbers[s], ".rds"))
    )
    
  }
  
  names(analyses_list) <- paste0("scenario_", scenario_numbers)
  class(analyses_list) <- "analysis_list"
  
  return (analyses_list)
  
}

mapUniqueTrials <- function (
    
  scenario_list,
  method_quantiles_list,
  trials_unique_calc,
  applicable_previous_trials
  
) {
  
  method_names     <- names(method_quantiles_list)
  scenario_numbers <- sapply(scenario_list, function (x) x$scenario_number)
  
  ## Create hash tables for results for easy retrieval
  
  hash_keys        <- getHashKeys(trials_unique_calc)
  hash_tables_list <- vector(mode = "list", length = length(method_quantiles_list))
  
  for (n in seq_along(hash_tables_list)) {
    
    hash_tables_list[[n]] <- createHashTable(hash_keys, method_quantiles_list[[n]])
    
  }
  
  ## prepare foreach
  exported_stuff <- c("convertVector2Matrix")
  
  ## run foreach
  "%do%" <- foreach::"%do%"
  scenario_method_quantiles_list <- foreach::foreach(k = seq_along(scenario_numbers),
                                                     .verbose  = FALSE,
                                                     .export   = exported_stuff
  ) %do% {
    
    ## Find the indices of the trials of a specific scenario for go trials
    scenario_data_matrix <- cbind(scenario_list[[k]]$n_responders,
                                  scenario_list[[k]]$n_subjects)
    
    ## check whether there where previous analyses
    if (applicable_previous_trials) {
      
      scenario_go_flags         <- scenario_list[[k]]$previous_analyses$go_decisions[, 1] > 0
      scenario_method_quantiles <- scenario_list[[k]]$previous_analyses$post_quantiles
      
    } else {
      
      scenario_go_flags         <- rep(TRUE, length = nrow(scenario_data_matrix))
      scenario_method_quantiles <- vector(mode = "list", length = length(method_names))
      names(scenario_method_quantiles) <- method_names
      
    }
    
    ## In case there are trial realizations that need updating
    ## This should only not be the case if all trial realizations of a scenario have a NoGo decision
    ## and there are applicable previous trials.
    if (any(scenario_go_flags)) {
      
      ## Get search keys
      scenario_data_matrix_go <- convertVector2Matrix(scenario_data_matrix[scenario_go_flags, ])
      search_keys             <- getHashKeys(scenario_data_matrix_go)
      
      ## Save scenario specific posterior quantiles for each method
      for (n in seq_along(method_names)) {
        
        scenario_method_quantiles[[method_names[n]]][scenario_go_flags] <-
          getHashValues(search_keys, hash_tables_list[[n]])
        
      }
      
    } 
    
    return (scenario_method_quantiles)
    
  }
  
  names(scenario_method_quantiles_list) <- paste0("scenario_", scenario_numbers)
  
  return (scenario_method_quantiles_list)
  
}

## Wrapper of getPostQuantiles() for specifying several scenarios
## Returns a list of quantiles of posterior distributions according to supplied methods.
#' @title performAnalyses
#' @md
#' @description This function performs the analysis of simulated or observed trial data with the
#' specified methods
#' and returns the quantiles of the posterior response rates
#' @param scenario_list An object of class `scenario_list`,
#' as e.g. created with \code{\link[bhmbasket]{simulateScenarios}}
#' @param evidence_levels A vector of numerics in `(0, 1)` for the
#' `1-evidence_levels`-quantiles of the posterior response rates to be saved.
#' Default: `c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975)`
#' @param method_names A vector of strings for the names of the methods to be used. Must
#' be one of the default values, Default: `c("berry", "exnex", "exnex_adj", "pooled", "stratified")`
#' @param target_rates A vector of numerics in `(0, 1)` for the
#' target rates of each cohort, Default: `NULL`
#' @param prior_parameters_list An object of class `prior_parameters_list`,
#' as e.g. created with \code{\link[bhmbasket]{getPriorParameters}}
#' @param calc_differences A matrix of positive integers with 2 columns.
#' For each row the differences will be calculated.
#' Also a vector of positive integers can be provided for a single difference.
#' The integers are the numbers for the cohorts to be subtracted from one another.
#' E.g. providing `c(2, 1)` calculates the difference between cohort `2` and cohort `1`.
#' If `NULL`, no subtractions are performed, Default: `NULL`
#' @param n_mcmc_iterations A positive integer for the number of MCMC iterations,
#' see Details, Default: `10000`.
#' If `n_mcmc_iterations` is present in `.GlobalEnv` and `missing(n_mcmc_iterations)`,
#' the globally available value will be used.
#' @param n_cores Argument is deprecated and does nothing as of version 0.9.3.
#' A positive integer for the number of cores for the parallelization,
#' Default: `1`
#' @param seed Argument is deprecated and does nothing as of version 0.9.3.
#' A numeric for the random seed, Default: `1`
#' @param verbose A logical indicating whether messages should be printed, Default: `TRUE`
#' @return An object of class `analysis_list`.
#' @details
#' This function applies the following analysis models to (simulated) scenarios of class
#' `scenario_list`:
#' \itemize{
#'   \item Bayesian hierarchical model (BHM) proposed by Berry et al. (2013): `"berry"`
#'   \item BHM proposed by Neuenschwander et al. (2016): `"exnex"`
#'   \item BHM that combines above approaches: `"exnex_adj"`
#'   \item Pooled beta-binomial approach: `"pooled"`
#'   \item Stratified beta-binomial approach: `"stratified"`
#' }
#' The posterior distributions of the BHMs are approximated with Markov chain Monte Carlo (MCMC)
#' methods implemented in JAGS.
#' Two independent chains are used with each `n_mcmc_iterations` number of MCMC iterations.
#' The first `floor(n_mcmc_iterations / 3)` number of iterations are discarded as burn-in period.
#' No thinning is applied.
#'
#' Note that the value for `n_mcmc_iterations` required for a good approximation of the posterior
#' distributions depends on the analysis model, the investigated scenarios, and the use case.
#' The default value might be a good compromise between run-time and approximation for
#' the estimation of decision probabilities, but
#' it should definitively be increased for the analysis of a single trial's outcome.
#' 
#' The analysis models will only be applied to the unique trial realizations across 
#' all simulated scenarios.
#' The models can be applied in parallel by registering a parallel backend for the 'foreach'
#' framework, e.g. with `doFuture::registerDoFuture()` and `future::plan(future::multisession)`.
#' The parallelization is nested, so that the resources of a HPC environment can be used
#' efficiently.
#' For more on this topic, kindly see the respective vignette.
#' The tasks that are to be performed in parallel are chunked according to the number of workers
#' determined with `foreach::getDoParWorkers()`.
#'
#' The JAGS code for the BHM `"exnex"` was taken from Neuenschwander et al. (2016).
#' The JAGS code for the BHM `"exnex_adj"` is based on the JAGS code for `"exnex"`.
#' @seealso
#'  \code{\link[bhmbasket]{simulateScenarios}}
#'  \code{\link[bhmbasket]{createTrial}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#' @rdname performAnalyses
#' @author Stephan Wojciekowski
#' @examples
#'  trial_data <- createTrial(
#'    n_subjects   = c(10, 20, 30),
#'    n_responders = c(1, 2, 3))
#'
#'  analysis_list <- performAnalyses(
#'    scenario_list      = trial_data,
#'    target_rates       = rep(0.5, 3),
#'    calc_differences   = matrix(c(3, 2, 1, 1), ncol = 2),
#'    n_mcmc_iterations  = 100)
#' @references Berry, Scott M., et al. "Bayesian hierarchical modeling of patient subpopulations:
#' efficient designs of phase II oncology clinical trials."
#' \emph{Clinical Trials} 10.5 (2013): 720-734.
#' @references Neuenschwander, Beat, et al. "Robust exchangeability designs
#' for early phase clinical trials with multiple strata."
#' \emph{Pharmaceutical statistics} 15.2 (2016): 123-134.
#' @references Plummer, Martyn. "JAGS: A program for analysis of Bayesian graphical models
#' using Gibbs sampling."
#' \emph{Proceedings of the 3rd international workshop on distributed statistical computing.}
#' Vol. 124. No. 125.10. 2003.
#' @export
performAnalyses <- function (
    
  scenario_list,
  evidence_levels       = c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975),
  
  method_names          = c("berry", "exnex", "exnex_adj", "pooled", "stratified"),
  target_rates          = NULL,
  prior_parameters_list = NULL,
  
  calc_differences      = NULL,
  
  n_mcmc_iterations     = 1e4,
  n_cores               = 1,
  seed                  = 1,
  verbose               = TRUE
  
) {
  
  error_scenario_list <-
    "Please provide an object of class scenario_list for the argument 'scenario_list'"
  error_evidence_levels <-
    "Please provide a vector of numerics in (0, 1) for the argument 'evidence_levels'"
  error_method_names <-
    paste("Please provide a (vector of) strings for the argument 'method_names'\n",
          "Must be one of 'berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'")
  error_target_rates <-
    paste("Please provide either 'NULL' or a vector of numerics in (0, 1)",
          "for the argument 'target_rates'")
  error_prior_parameters_list <-
    "Please provide either 'NULL' or an object of class 'prior_parameters_list'"
  error_calc_differences <-
    paste("Please provide either 'NULL' or a matrix of integers with ncol = 2.",
          "The values of the integers must be less than or equal to the number",
          "of cohorts")
  error_n_mcmc_iterations <-
    "Please provide a positive integer for the argument 'n_mcmc_iterations'"
  error_verbose <-
    "Please provide a logical for the argument 'verbose'"
  
  error_need_target_rates_for_methods <-
    "Please provide 'target_rates' when using the methods 'berry' and/or 'exnex_adj'"
  error_need_one_of_prior_or_target <-
    "Please provide at least one of 'prior_parameters_list' or 'target_rates'"
  error_target_length <-
    "The length of 'target_rates' does not match the number of cohorts"
  error_prior_methods_mismatch <-
    paste("Not all specified methods in 'method_names'",
          "have prior parameters specified in 'prior_parameters_list'")
  error_prior_cohorts_mismatch <-
    paste("The number of cohorts specified in 'prior_parameters_list' does not match",
          "the number of cohorts specified in 'scenario_list'")
  
  warning_n_cores <- "The argument 'n_cores' is deprecated as of version 0.9.3."
  warning_seed    <- "The argument 'seed' is deprecated as of version 0.9.3."
  
  if (!missing(n_cores)) warning(warning_n_cores)
  if (!missing(seed))    warning(warning_seed)
  
  # scenario_list must be supplied and have the correct class
  checkmate::assertClass(
    scenario_list,
    "scenario_list",
    .var.name = error_scenario_list
  )
  
  # evidence_levels: numeric, all in (0, 1)
  checkmate::assertNumeric(
    evidence_levels,
    any.missing = FALSE,
    .var.name   = error_evidence_levels
  )
  checkmate::assertTRUE(
    all(evidence_levels > 0 & evidence_levels < 1),
    .var.name = error_evidence_levels
  )
  
  # method_names: character and must be one of the allowed methods
  checkmate::assertCharacter(
    method_names,
    any.missing = FALSE,
    min.len     = 1,
    .var.name   = error_method_names
  )
  method_names <- tryCatch(
    match.arg(
      method_names,
      choices    = c("berry", "exnex", "exnex_adj", "pooled", "stratified"),
      several.ok = TRUE
    ),
    error = function(e) {
      stop(error_method_names, call. = FALSE)
    }
  )
  
  # target_rates: either NULL or numeric in (0, 1)
  if (!is.null(target_rates)) {
    checkmate::assertNumeric(
      target_rates,
      any.missing = FALSE,
      .var.name   = error_target_rates
    )
    checkmate::assertTRUE(
      all(target_rates > 0 & target_rates < 1),
      .var.name = error_target_rates
    )
  }
  
  # prior_parameters_list: either NULL or object with class 'prior_parameters_list'
  if (!is.null(prior_parameters_list)) {
    checkmate::assertClass(
      prior_parameters_list,
      "prior_parameters_list",
      .var.name = error_prior_parameters_list
    )
  }
  
  
  if (is.null(target_rates)) {
    
    # If target_rates are missing, berry and exnex_adj are not allowed
    checkmate::assertTRUE(
      !any(c("berry", "exnex_adj") %in% method_names),
      .var.name = error_need_target_rates_for_methods
    )
    
    # Need at least one of prior_parameters_list or target_rates
    checkmate::assertTRUE(
      !is.null(prior_parameters_list),
      .var.name = error_need_one_of_prior_or_target
    )
    
  } else {
    # If target_rates is present, its length must match the number of cohorts
    n_coh <- ncol(scenario_list[[1]]$n_subjects)
    checkmate::assertTRUE(
      identical(length(target_rates), n_coh),
      .var.name = error_target_length
    )
  }
  
  
  if (!is.null(prior_parameters_list)) {
    
    # All methods used must have entries in prior_parameters_list
    checkmate::assertTRUE(
      all(method_names %in% names(prior_parameters_list)),
      .var.name = error_prior_methods_mismatch
    )
    
    # For exnex / exnex_adj / stratified, per-cohort prior lengths must
    # match the number of cohorts in scenario_list
    n_coh <- ncol(scenario_list[[1]]$n_subjects)
    
    inconsistent_cohorts <- any(
      sapply(
        intersect(names(prior_parameters_list),
                  c("exnex", "exnex_adj", "stratified")),
        function(name) {
          max(sapply(prior_parameters_list[[name]], length)) != n_coh
        }
      )
    )
    
    checkmate::assertFALSE(
      inconsistent_cohorts,
      .var.name = error_prior_cohorts_mismatch
    )
  }
  
  n_cohorts_min <- min(sapply(scenario_list, function(x) ncol(x$n_responders)))
  
  if (!is.null(calc_differences)) {
    
    checkmate::assertNumeric(
      calc_differences,
      any.missing = FALSE,
      .var.name   = error_calc_differences
    )
    
    is_len2   <- identical(length(calc_differences), 2L)
    has_2cols <- !is.null(dim(calc_differences)) &&
      identical(ncol(calc_differences), 2L)
    
    checkmate::assertTRUE(
      is_len2 || has_2cols,
      .var.name = error_calc_differences
    )
    
    checkmate::assertIntegerish(
      calc_differences,
      lower       = 1,
      any.missing = FALSE,
      .var.name   = error_calc_differences
    )
    
    checkmate::assertTRUE(
      max(calc_differences) <= n_cohorts_min,
      .var.name = error_calc_differences
    )
  }
  rm(n_cohorts_min)
  
  # If n_mcmc_iterations exists in .GlobalEnv and argument is missing, reuse it
  if ("n_mcmc_iterations" %in% ls(envir = .GlobalEnv) && missing(n_mcmc_iterations)) {
    n_mcmc_iterations <- get("n_mcmc_iterations", envir = .GlobalEnv)
  }
  
  checkmate::assertInt(
    n_mcmc_iterations,
    lower     = 1,
    .var.name = error_n_mcmc_iterations
  )
  
  checkmate::assertLogical(
    verbose,
    len         = 1L,
    any.missing = FALSE,
    .var.name   = error_verbose
  )
  
  # Only need a parallel backend if there is at least one method that is not pooled/stratified
  if (!all(method_names %in% c("stratified", "pooled"))) {
    checkForParallelBackend()
  }
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  ## message to user
  if (verbose) message(format(Sys.time(), "%d-%h-%Y"), " Performing Analyses")
  
  ## some housekeeping
  method_names <- sort(method_names)
  quantiles    <- sort(unique(round(1 - c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975,
                                          evidence_levels), 9)))
  if (!is.null(calc_differences)) {
    calc_differences <- convertVector2Matrix(calc_differences)
  }
  
  ## get scenario numbers
  scenario_numbers <- sapply(scenario_list, function (x) x$scenario_number)
  
  ## get unique trials over all scenarios
  trials_unique <- getUniqueTrials(scenario_list)
  n_cohorts     <- (ncol(trials_unique) - 1L) / 2
  
  ## analyze only unique trials that have not been previously analyzed
  applicable_previous_trials <- applicablePreviousTrials(
    scenario_list    = scenario_list,
    method_names     = method_names,
    quantiles        = quantiles,
    n_cohorts        = n_cohorts,
    calc_differences = calc_differences)
  
  ## only get previous go indices if all conditions for previous trials are met
  if (applicable_previous_trials) {
    calc_trial_indices <- trials_unique[, ncol(trials_unique)] > 0
  } else {
    calc_trial_indices <- rep(TRUE, nrow(trials_unique))
  }
  
  ## get resulting unique number of responders and number of subjects
  trials_unique_calc <- trials_unique[calc_trial_indices, -ncol(trials_unique)]
  n_responders       <- trials_unique_calc[, seq_len(n_cohorts)]
  n_subjects         <- trials_unique_calc[, seq_len(n_cohorts) + n_cohorts]
  
  ## message to user
  if (verbose) {
    
    message("         Analyzing ", length(scenario_numbers) ," scenario", ifelse(length(scenario_numbers) == 1, "", "s")," ",
            "(", nrow(n_responders), " unique", ifelse(applicable_previous_trials, " updated ", " "),
            "trial realization", ifelse(nrow(trials_unique) == 1, "", "s"),")")
    
  }
  
  ## get default prior parameters if needed
  if(is.null(prior_parameters_list)) {
    
    prior_parameters_list <- getPriorParameters(
      method_names = method_names,
      target_rates = target_rates)
    
  }
  
  ## Create lists to save all results of the unique trials
  method_quantiles_list        <- vector(mode = "list", length = length(method_names))
  names(method_quantiles_list) <- method_names
  
  ## For each method
  for (n in seq_along(method_names)) {
    
    ## message to user
    if (verbose) {
      start_time  <- Sys.time()
      out_message <- paste0(format(start_time, "   %H:%M", digits = 1),
                            " - with ", firstUpper(method_names[n]), " ...")
      message(out_message, rep(".", 33 - nchar(out_message)))
    }
    
    ## prepare analysis
    prepare_analysis <- prepareAnalysis(
      method_name       = method_names[n],
      target_rates      = target_rates,
      prior_parameters  = prior_parameters_list[[method_names[n]]])
    
    ## run analysis
    method_quantiles_list[[method_names[n]]] <- getPostQuantiles(
      method_name       = method_names[n],
      quantiles         = quantiles,
      scenario_data     = list(n_subjects   = n_subjects,
                               n_responders = n_responders),
      calc_differences  = calc_differences,
      j_parameters      = prepare_analysis$j_parameters,
      j_model_file      = prepare_analysis$j_model_file,
      j_data            = prepare_analysis$j_data,
      n_mcmc_iterations = n_mcmc_iterations,
      save_path         = NULL,
      save_trial        = NULL)
    
    ## message to user
    if (verbose) {
      message("             finished after ", round(Sys.time() - start_time, 1), " ",
              units(Sys.time() - start_time), ".")
      rm(start_time)
    }
    
  }
  
  ## message to user
  if (verbose) {
    start_time  <- Sys.time()
    message("         Processing scenarios ...")
  }
  
  ## Process scenarios
  scenario_method_quantiles_list <- mapUniqueTrials(
    scenario_list              = scenario_list,
    method_quantiles_list      = method_quantiles_list,
    trials_unique_calc         = trials_unique_calc,
    applicable_previous_trials = applicable_previous_trials)
  
  ## message to user
  if (verbose) {
    message("             finished after ", round(Sys.time() - start_time, 1), " ",
            units(Sys.time() - start_time), ".")
    rm(start_time)
  }
  
  ## combine results from all scenarios & return
  analyses_list        <- vector(mode = "list", length = length(scenario_numbers))
  names(analyses_list) <- paste0("scenario_", scenario_numbers)
  
  for (s in seq_along(scenario_numbers)) {
    
    analyses_list[[s]] <- list(
      quantiles_list      = scenario_method_quantiles_list[[s]],
      scenario_data       = scenario_list[[s]],
      analysis_parameters = list(
        quantiles             = quantiles,
        method_names          = method_names,
        prior_parameters_list = prior_parameters_list,
        n_mcmc_iterations     = n_mcmc_iterations))
    
  }
  
  class(analyses_list) <- "analysis_list"
  
  return (analyses_list)
  
}

## based on R2jags::jags
## stripped down to improve performance
performJags <- function (
    
  data,
  parameters_to_save,
  model_file, 
  n_chains = 2,
  n_iter   = 1e4,
  n_burnin = floor(n_iter/3)
  
) {
  
  n_adapt <- ifelse(n_burnin > 0, n_burnin, 100)
  
  inits <- vector("list", n_chains)
  for (i in 1:n_chains) {
    inits[[i]]$.RNG.name <- "base::Wichmann-Hill"
    inits[[i]]$.RNG.seed <- stats::runif(1, 0, 2^31)
  }
  
  j_model <- rjags::jags.model(file     = model_file,
                               data     = data,
                               inits    = inits, 
                               n.chains = n_chains,
                               n.adapt  = 0,
                               quiet    = TRUE)
  
  rjags::adapt(object         = j_model,
               n.iter         = n_adapt,
               progress.bar   = "none",
               end.adaptation = TRUE)
  
  samples <- rjags::coda.samples(model          = j_model,
                                 variable.names = parameters_to_save, 
                                 n.iter         = n_iter - n_burnin,
                                 thin           = 1, 
                                 progress.bar   = "none")
  
  return(do.call(rbind, samples))
  
}

posteriors2Quantiles <- function (
    
  quantiles,
  posteriors
  
) {
  
  posterior_quantiles <- apply(posteriors, 2, function (x) stats::quantile(x, probs = quantiles))
  
  posterior_mean      <- apply(posteriors, 2, mean)
  posterior_sd        <- apply(posteriors, 2, stats::sd)
  posterior_quantiles <- rbind(posterior_quantiles,
                               Mean = posterior_mean,
                               SD   = posterior_sd)
  
  return (posterior_quantiles)
  
}

prepareAnalysis <- function (
    
  method_name,
  
  prior_parameters = NULL,
  target_rates     = NULL
  
) {
  
  if (method_name == "berry") {
    
    j_data <- list(mean_mu       = prior_parameters$mu_mean,
                   precision_mu  = prior_parameters$mu_sd^-2,
                   precision_tau = prior_parameters$tau_scale^-2,
                   p_t           = target_rates,
                   J             = length(target_rates))
    
    # j_model_file <- writeTempModel(method_name = "berry")
    j_model_file <- getModelFile(method_name = "berry")
    
    j_parameters <- c("p", "mu", "tau")
    
  } else if (method_name == "exnex" | method_name == "exnex_adj") {
    # Nexch: number of exchangeable mixture components
    # Nmix:  number of mixture components
    j_data <- list(Nexch        = length(prior_parameters$mu_mean),
                   Nmix         = length(prior_parameters$mu_mean) + 1L,
                   Nstrata      = length(prior_parameters$mu_j),
                   mu_mean      = prior_parameters$mu_mean,
                   mu_prec      = prior_parameters$mu_sd^-2,
                   tau_HN_scale = rep(prior_parameters$tau_scale,
                                      length(prior_parameters$mu_mean)), # rep(..., Nexch)
                   nex_mean     = prior_parameters$mu_j,
                   nex_prec     = prior_parameters$tau_j^-2)
    
    if (identical(length(prior_parameters$w_j), 1L)) {
      j_data$pMix <- c(prior_parameters$w_j, 1 - prior_parameters$w_j)
    } else {
      j_data$pMix <- prior_parameters$w_j
    }
    
    if (method_name == "exnex") {
      
      j_model_file    <- getModelFile(method_name = "exnex")
      
    } else {
      
      j_data$p_target <- target_rates
      j_model_file    <- getModelFile(method_name = "exnex_adj")
      
    }
    
    j_parameters <- c("p", "mu", "tau", "exch")
    
  } else if (method_name == "stratified" | method_name == "pooled") {
    
    ## For methods "stratified" and "pooled" no MCMC simulations are necessary,
    ## as the posterior response rates of the cohorts follow known beta distributions.
    
    j_model_file <- "dummy path to JAGS model"
    j_parameters <- "dummy JAGS parameters"
    j_data       <- prior_parameters
    
  } else {
    
    stop ("method_name must be one of berry, exnex, exnex_adj, stratified, pooled")
    
  }
  
  return (list(j_parameters = j_parameters,
               j_model_file = j_model_file,
               j_data       = j_data))
}

#' @export
print.analysis_list <- function (x, digits = 2, ...) {
  
  n_scenarios    <- length(x)
  scenario_names <- names(x)
  
  n_methods      <- length(x[[1]]$quantiles_list)
  method_names   <- names(x[[1]]$quantiles_list)
  
  estimates          <- getEstimates(x)
  n_mcmc_interations <- x[[1]]$analysis_parameters$n_mcmc_iterations
  
  evidence_levels <- sort(1 - x[[1]]$analysis_parameters$quantiles)
  
  cat("analysis_list of ", n_scenarios, " scenario", ifelse(n_scenarios == 1, "", "s"),
      " with ", n_methods, " method", ifelse(n_methods == 1, "", "s"),"\n\n", sep = "")
  
  for (n in seq_along(scenario_names)) {
    
    if (n_scenarios == 1L) {
      
      expr <- quote(t(y[, 1:2]))
      
    } else {
      
      expr <- quote(t(y[[n]][, 1:2]))
      
    }
    
    mat_out <- do.call(rbind, lapply(estimates, function (y) eval(expr)))
    
    rownames(mat_out) <-  paste0(
      c("    - ", "      "),
      c(rbind(
        paste0(
          firstUpper(method_names),
          sapply(method_names, function (y) {
            getBlankString(max(nchar(method_names)) - nchar(y) + 1)
          })),
        rep(getBlankString(max(nchar(method_names)) + 3),
            length = length(method_names))
      )),
      rownames(mat_out))
    
    cat("  -", scenario_names[n], "\n")
    print(round(mat_out, digits = digits))
    
    cat("\n")
    
  }
  
  cat("  -", n_mcmc_interations, "MCMC iterationns per BHM method\n")
  cat("  - Available evidence levels:", evidence_levels, "\n")
  
}

qbetaDiff <- function (
    
  quantiles,
  
  x_1_shape1,
  x_1_shape2,
  
  x_2_shape1,
  x_2_shape2,
  
  n_mcmc = 1e6
  
) {
  
  sample_1   <- stats::rbeta(n_mcmc, shape1 = x_1_shape1, shape2 = x_1_shape2)
  sample_2   <- stats::rbeta(n_mcmc, shape1 = x_2_shape1, shape2 = x_2_shape2)
  difference <- sample_1 - sample_2
  
  quantiles_diff <- stats::quantile(difference, probs = quantiles)
  
  mean_diff      <- mean(difference)
  sd_diff        <- stats::sd(difference)
  quantiles_diff <- c(quantiles_diff, mean_diff, sd_diff)
  
  return (quantiles_diff)
  
}


#' @title saveAnalyses
#' @md
#' @description This function saves an object of class `analysis_list`
#' @param analyses_list An object of class `analysis_list`,
#' as created with \code{\link[bhmbasket]{performAnalyses}}
#' @param save_path A string for the path where the scenarios are being stored,
#' Default: \code{\link[base]{tempfile}}
#' @param analysis_numbers A positive integer naming the analysis number.
#' If `NULL`, the function will look for the number of saved analyses of the scenario
#' in the directory and add 1, Default: `NULL`
#' @return A named list of length 3 of vectors with scenario and analysis numbers and
#' the `save_path`
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{loadAnalyses}}
#'  \code{\link[base]{tempfile}}
#' @rdname saveAnalyses
#' @examples
#'   trial_data <- createTrial(
#'     n_subjects   = c(10, 20, 30),
#'     n_responders = c(1, 2, 3))
#'
#'   analysis_list <- performAnalyses(
#'     scenario_list      = trial_data,
#'     target_rates       = rep(0.5, 3),
#'     n_mcmc_iterations  = 100)
#'
#'   save_info     <- saveAnalyses(analysis_list)
#'   analysis_list <- loadAnalyses(scenario_numbers = save_info$scenario_numbers,
#'                                 analysis_numbers = save_info$analysis_numbers,
#'                                 load_path        = save_info$path)
#' @author Stephan Wojciekowski
#' @export
saveAnalyses <- function (
    
  analyses_list,
  save_path        = tempdir(),
  analysis_numbers = NULL
  
) {
  
  error_analyses_list <- 
    "Please provide an object of class analysis_list for the argument 'analyses_list'"
  error_save_path     <- 
    "Please provide a string containing a path for the argument 'save_path'"
  error_analysis_numbers <- paste(
    "Please provide a vector of positive integers for the argument 'analysis_numbers'",
    "with length equal to the length of 'analyses_list'"
  )
  
  checkmate::assertClass(
    analyses_list,
    "analysis_list",
    .var.name = error_analyses_list
  )
  checkmate::assertCharacter(
    save_path,
    len         = 1,
    any.missing = FALSE,
    .var.name   = error_save_path
  )
  
  if (!is.null(analysis_numbers)) {
    checkmate::assertIntegerish(
      analysis_numbers, 
      lower       = 1, 
      any.missing = FALSE,
      .var.name   = error_analysis_numbers
    )
    
    checkmate::assertTRUE(
      identical(length(analyses_list), length(analysis_numbers)),
      .var.name = error_analysis_numbers
    )
  }
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  scenario_numbers <- sapply(analyses_list, function (x) x$scenario_data$scenario_number)
  
  if (is.null(analysis_numbers)) {
    analysis_numbers <- rep(0, length(analyses_list))
  }
  
  for (s in seq_along(analyses_list)) {
    
    ## Get analysis number
    if (identical(analysis_numbers[s], 0)) {
      
      analysis_numbers[s] <- sum(grepl(paste0("analysis_data_", scenario_numbers[s], "_"),
                                       list.files(save_path))) + 1L
      
    }
    
    ## Save the analysis
    file_name <- paste0("analysis_data_", scenario_numbers[s], "_", analysis_numbers[s],".rds")
    saveRDS(analyses_list[[s]], file = file.path(save_path, file_name),
            compress = "xz")
    
  }
  
  return (list(scenario_numbers = scenario_numbers,
               analysis_numbers = analysis_numbers,
               path             = save_path))
  
}