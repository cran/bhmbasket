# library(sinew)
# makeOxygen(LoadScenario)

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

    j_model_file <- writeTempModel(method_name = "berry")

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

      j_model_file    <- writeTempModel(method_name = "exnex")

    } else {

      j_data$p_target <- target_rates

      j_model_file    <- writeTempModel(method_name = "exnex_adj")

    }

    j_parameters <- c("p", "mu", "tau", "exch")

  } else if (method_name == "stratified" | method_name == "pooled") {

    ## For methods "stratified" and "pooled" no MCMC simulations are necessary,
    ## as the posterior response rates of the cohorts follow known beta distributions.

    j_model_file <- "dummy path to model"
    j_parameters <- "dummy parameters"
    j_data       <- prior_parameters

  } else {

    stop ("method_name must be one of berry, exnex, exnex_adj, stratified, pooled")

  }

  return (list(j_parameters = j_parameters,
               j_model_file = j_model_file,
               j_data       = j_data))
}

getPosteriors <- function (

  j_parameters,
  j_model_file,
  j_data,

  n_mcmc_iterations

) {

  jags_fit <- R2jags::jags(data               = j_data,
                           parameters.to.save = j_parameters,
                           model.file         = j_model_file,
                           n.chains           = 2,
                           n.iter             = n_mcmc_iterations,
                           n.burnin           = floor(n_mcmc_iterations / 3),
                           n.thin             = 1,
                           DIC                = FALSE,
                           progress.bar       = "none",
                           jags.module        = "mix")

  ## Adaption and burn-in not included in sims.array
  posterior_distributions <- rbind(jags_fit$BUGSoutput$sims.array[, 1, ],
                                   jags_fit$BUGSoutput$sims.array[, 2, ])

  colnames(posterior_distributions) <- gsub("\\[", "_", colnames(posterior_distributions))
  colnames(posterior_distributions) <- gsub("\\]", "", colnames(posterior_distributions))

  weights_indices <- grepl("exch", colnames(posterior_distributions))
  if (any(weights_indices)) {

    superfluous_weights <- !grepl(",1", colnames(posterior_distributions))

    colnames(posterior_distributions)[weights_indices] <-
      paste0("w_", seq_along(j_data$n))

    posterior_distributions <- posterior_distributions[, !(weights_indices & superfluous_weights)]

  }

  return (posterior_distributions)

}

posteriors2Quantiles <- function (

  quantiles,
  posteriors,

  post_mean = TRUE

) {

  posterior_quantiles <- apply(posteriors, 2, function (x) stats::quantile(x, probs = quantiles))

  if (post_mean) {

    posterior_mean      <- apply(posteriors, 2, mean)
    posterior_sd        <- apply(posteriors, 2, stats::sd)
    posterior_quantiles <- rbind(posterior_quantiles,
                                 Mean = posterior_mean,
                                 SD   = posterior_sd)

  }

  return (posterior_quantiles)

}

qbetaDiff <- function (

  quantiles,
  post_mean,

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

  if (post_mean) {

    mean_diff      <- mean(difference)
    sd_diff        <- stats::sd(difference)
    quantiles_diff <- c(quantiles_diff, mean_diff, sd_diff)

  }

  return (quantiles_diff)

}

getPostQuantiles <- function (

  ## The method to be applied to the likelihood and the quantiles of the posterior
  method_name,
  quantiles,
  post_mean,

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
  seed              = as.numeric(Sys.time()),
  n_cores           = parallel::detectCores() - 1L,

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
    set.seed(seed)
    save_trial <- sample(seq_len(n_analyses), size = 1)
  }

  ## Run parallel loops
  ## prepare foreach loop over k
  "%do%"    <- foreach::"%do%"
  "%dopar%" <- foreach::"%dopar%"

  if (isTRUE(all.equal(n_analyses, 1))) {
    n_cores <- 1L
  }

  exported_stuff <- c("seed", "scenario_data", "j_data", "post_mean", "quantiles",
                      "calc_differences", "j_parameters", "j_model_file", "n_mcmc_iterations",
                      "getPosteriors", "save_path", "save_trial", "method_name",
                      "posteriors2Quantiles",
                      "qbetaDiff")

  foreach_expression <- quote({

    ## Set seed that changes for each iteration
    set.seed(seed + k)

    ##  Retrieve the likelihood data for the kth unique simulation
    j_data$r <- as.numeric(scenario_data$n_responders[k, ])
    j_data$n <- as.numeric(scenario_data$n_subjects[k, ])

    if (method_name == "stratified") {

      shape_1 <- j_data$a_j + j_data$r
      shape_2 <- j_data$b_j + j_data$n - j_data$r

      posterior_quantiles <- t(sapply(quantiles, function (x)
        stats::qbeta(x, shape1 = shape_1, shape2 = shape_2)))

      if (nrow(posterior_quantiles) == 1) {

        posterior_quantiles <- t(posterior_quantiles)

      }

      colnames(posterior_quantiles) <- paste0("p_", seq_along(j_data$a_j))
      rownames(posterior_quantiles) <- paste0(quantiles * 100, "%")

      if (post_mean) {

        posterior_mean      <- shape_1 / (shape_1 + shape_2)
        posterior_sd        <- shape_1 * shape_2 / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1))
        posterior_quantiles <- rbind(posterior_quantiles,
                                     Mean = posterior_mean,
                                     SD   = posterior_sd)

      }

      if (!is.null(calc_differences)) {

        diff_quantiles <- apply(calc_differences, 1, function (x) {

          matrix(qbetaDiff(
            quantiles  = quantiles,
            post_mean  = post_mean,
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

    } else if (method_name == "pooled") {

      shape_1 <- j_data$a + sum(j_data$r)
      shape_2 <- j_data$b + sum(j_data$n) - sum(j_data$r)

      posterior_quantiles <- stats::qbeta(quantiles, shape1 = shape_1, shape2 = shape_2)

      posterior_quantiles <- matrix(posterior_quantiles,
                                    ncol = length(j_data$r), nrow = length(quantiles))

      colnames(posterior_quantiles) <- paste0("p_", seq_along(j_data$r))
      rownames(posterior_quantiles) <- paste0(quantiles * 100, "%")

      if (post_mean) {

        posterior_mean      <- shape_1 / (shape_1 + shape_2)
        posterior_sd        <- shape_1 * shape_2 / ((shape_1 + shape_2)^2 * (shape_1 + shape_2 + 1))
        posterior_quantiles <- rbind(posterior_quantiles,
                                     Mean = posterior_mean,
                                     SD   = posterior_sd)

      }

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

    } else {

      ## Calculate posterior response rates per indication
      posterior_distributions <- getPosteriors(j_parameters      = j_parameters,
                                               j_model_file      = j_model_file,
                                               j_data            = j_data,
                                               n_mcmc_iterations = n_mcmc_iterations)

      ## Calculate differences between response rates of 2nd and 1st cohort
      if (!is.null(calc_differences)) {

        org_names <- colnames(posterior_distributions)

        diffs <- apply(calc_differences, 1, function (x) {

          matrix(posterior_distributions[, grepl(x[1], org_names) & grepl("p", org_names)] -
            posterior_distributions[, grepl(x[2], org_names) & grepl("p", org_names)],
            ncol = 1)

        })

        diff_names <- apply(calc_differences, 1, function (x) {

          paste0("p_diff_", paste0(as.character(x), collapse = ""))

        })

        colnames(diffs) <- diff_names

        posterior_distributions <- cbind(posterior_distributions, diffs)

      }

      ## Save posterior response rates per indication for one randomly selected simulation,
      ## due to time and storage space constraints only one simulation
      if (!is.null(save_path)) {
        if (k == save_trial) {
          saveRDS(posterior_distributions,
                  file = file.path(save_path, paste0("posterior_distributions_",
                                                     k, "_", method_name, "_rds")))
        }
      }

      ## Calculate the required quantiles for the decision rules
      posterior_quantiles <- posteriors2Quantiles(quantiles  = quantiles,
                                                  post_mean  = post_mean,
                                                  posteriors = posterior_distributions)

    }

    return (posterior_quantiles)

  })

  if (n_cores > 1) {

    doParallel::registerDoParallel(n_cores)
    on.exit(doParallel::stopImplicitCluster())

    posterior_quantiles_list <- foreach::foreach(k = seq_len(n_analyses),
                                                 .verbose  = FALSE,
                                                 .packages = c("R2jags"),
                                                 .export   = exported_stuff
    ) %dopar% {eval(foreach_expression)}

  } else {

    invisible(utils::capture.output({

      suppressMessages({

        posterior_quantiles_list <- foreach::foreach(k = seq_len(n_analyses),
                                                     .verbose  = FALSE,
                                                     .packages = c("R2jags"),
                                                     .export   = exported_stuff
        ) %do% {eval(foreach_expression)}

      })

    }))

  }

  return (posterior_quantiles_list)

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
#' Default: `10000`
#' @param n_cores A positive integer for the number of cores for the parallelization,
#' Default: `parallel::detectCores() - 1L`
#' @param seed A numeric for the random seed, Default: `as.numeric(Sys.time())`
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
#' The posterior distributions of the BHMs are approximated with Markov chain Monte Carlo methods.
#' The MCMC methods are implemented in JAGS.
#' The JAGS code for the BHM `"exnex"` was taken from Neuenschwander et al. (2016).
#' The JAGS code for the BHM `"exnex_adj"` is based on the JAGS code for `"exnex"`.
#' @seealso
#'  \code{\link[parallel]{detectCores}}
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
#'    n_mcmc_iterations  = 100,
#'    n_cores            = 1L)
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
#' @importFrom parallel detectCores
#' @export
performAnalyses <- function (

  scenario_list,
  evidence_levels       = c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975),

  method_names          = c("berry", "exnex", "exnex_adj", "pooled", "stratified"),
  target_rates          = NULL,
  prior_parameters_list = NULL,

  calc_differences      = NULL,

  n_mcmc_iterations     = 1e4,
  n_cores               = parallel::detectCores() - 1L,
  seed                  = as.numeric(Sys.time()),
  verbose               = TRUE

) {

  error_scenario_list <-
    simpleError("Please provide an object of class scenario_list for the argument 'scenario_list'")
  error_evidence_levels <-
    simpleError("Please provide a vector of numerics in (0, 1) for the argument 'evidence_levels'")
  error_method_names <-
    simpleError(paste("Please provide a (vector of) strings for the argument 'method_names'\n",
                      "Must be one of 'berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'"))
  error_target_rates <-
    simpleError(paste("Please provide either 'NULL' or a vector of numerics in (0, 1)",
                      "for the argument 'target_rates'"))
  error_prior_parameters_list <-
    simpleError("Please provide either 'NULL' or an object of class 'prior_parameters_list'")
  error_calc_differences <-
    simpleError(paste("Please provide either 'NULL' or a matrix of integers with ncol = 2.",
                      "The values of the integers must be less than or equal to the number",
                      "of cohorts"))
  error_n_mcmc_iterations <-
    simpleError("Please provide a positive integer for the argument 'n_mcmc_iterations'")
  error_n_cores <-
    simpleError("Please provide a positive integer for the argument 'n_cores'")
  error_seed <-
    simpleError("Please provide a numeric for the argument 'seed'")
  error_verbose <-
    simpleError("Please provide a logical for the argument 'verbose'")

  if (missing(scenario_list)) stop (error_scenario_list)

  method_names <- tryCatch({

    match.arg(
      method_names,
      choices    = c('berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'),
      several.ok = TRUE)

  }, error = function (e) e)

  if (!is.scenario_list(scenario_list))                   stop (error_scenario_list)
  if (!is.numeric.in.zero.one(evidence_levels))           stop (error_evidence_levels)
  if (inherits(method_names, "error"))                    stop (error_method_names)
  if (!is.null(target_rates) &&
      !is.numeric.in.zero.one(target_rates))              stop (error_target_rates)
  if (!is.null(prior_parameters_list) &&
      !is.prior_parameters_list(prior_parameters_list))   stop (prior_parameters_list)

  if (is.null(target_rates)) {
    if (any(c("berry", "exnex_adj") %in% method_names)) stop (simpleError(
      "Please provide 'target_rates' when using the methods 'berry' and/or 'exnex_adj'"))
    if (is.null(prior_parameters_list)) stop (simpleError(
      "Please provide at least one of 'prior_parameters_list' or 'target_rates'"))
  } else {
    if (!identical(length(target_rates), ncol(scenario_list[[1]]$n_subjects))) stop (simpleError(
      "The length of 'target_rates' does not match the number of cohorts"))
  }

  if (!is.null(prior_parameters_list)) {
    if (!all(method_names %in% names(prior_parameters_list))) stop (simpleError(
      paste("Not all specified methods in 'method_names'",
            "have prior parameters specified in 'prior_parameters_list'")))
    if (any(sapply(names(prior_parameters_list), function (name) {
      if (name %in% c('exnex', 'exnex_adj', 'stratified')) {
        !identical(max(sapply(prior_parameters_list[[name]], length)),
                   ncol(scenario_list[[1]]$n_subjects))
      } else FALSE
    }))) stop (simpleError(paste(
      "The number of cohorts specified in 'prior_parameters_list' does not match",
      "the number of cohorts specified in 'scenario_list'")))
  }

  n_cohorts_min <- min(sapply(scenario_list, function (x) {
    ncol(x$n_responders)
  }))
  if (!is.null(calc_differences) && (
       !is.numeric(calc_differences) ||
       !(identical(length(calc_differences), 2L) ||
         identical(ncol(calc_differences), 2L)) ||
       !is.positive.wholenumber(calc_differences) ||
       max(calc_differences) > n_cohorts_min))            stop (error_calc_differences)
  rm (n_cohorts_min)

  if (!is.single.positive.wholenumber(n_mcmc_iterations)) stop (error_n_mcmc_iterations)
  if (!is.single.positive.wholenumber(n_cores))           stop (error_n_cores)
  if (!is.single.numeric(seed))                           stop (error_seed)
  if (!is.logical(verbose))                               stop (error_verbose)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  ## message to user
  if (verbose) message(format(Sys.time(), "%m/%d/%y"), " Performing Analyses")

  ## some housekeeping
  method_names <- sort(method_names)
  quantiles    <- sort(unique(round(1 - c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975,
                                          evidence_levels), 9)))
  if (!is.null(calc_differences)) {
    calc_differences <- convertVector2Matrix(calc_differences)
  }

  ## get scenario numbers
  scenario_numbers <- sapply(scenario_list, function (x) x$scenario_number)

  ## get unique trials for all scenarios
  all_scenarios_n_responders <- do.call(rbind,
                                        lapply(scenario_list, function (x) x$n_responders))
  all_scenarios_n_subjects   <- do.call(rbind,
                                        lapply(scenario_list, function (x) x$n_subjects))
  all_scenarios_overall_gos  <-
    do.call(rbind, lapply(scenario_list, function (x) x$previous_analyses$go_decisions))[, 1]

  n_cohorts     <- ncol(all_scenarios_n_responders)
  trials_unique <- makeUniqueRows(cbind(all_scenarios_n_responders, all_scenarios_n_subjects,
                                        go_flag = all_scenarios_overall_gos))

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

  if (applicable_previous_trials) {

    go_trial_indices <- trials_unique[, ncol(trials_unique)] > 0

  } else {

    go_trial_indices <- rep(TRUE, nrow(trials_unique))

  }

  trials_unique <- trials_unique[, -ncol(trials_unique)]

  ## save resulting unique number of responders and number of subjects
  n_responders  <- trials_unique[go_trial_indices, seq_len(n_cohorts)]
  n_subjects    <- trials_unique[go_trial_indices, seq_len(n_cohorts) + n_cohorts]

  ## message to user
  if (verbose) {

    scenarios_message <- paste0(as.character(scenario_numbers), sep = ", ", collapse = "")
    message("         Analyzing Scenario", ifelse(length(scenario_numbers) == 1, "", "s")," ",
            substr(scenarios_message, 1, nchar(scenarios_message) - 2),
            " (", nrow(n_responders), " unique", ifelse(applicable_previous_trials, " updated ", " "),
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
      post_mean         = TRUE,
      scenario_data     = list(n_subjects   = n_subjects,
                               n_responders = n_responders),
      calc_differences  = calc_differences,
      j_parameters      = prepare_analysis$j_parameters,
      j_model_file      = prepare_analysis$j_model_file,
      j_data            = prepare_analysis$j_data,
      n_mcmc_iterations = n_mcmc_iterations,
      seed              = seed,
      n_cores           = n_cores,
      save_path         = NULL,
      save_trial        = NULL)

    ## message to user
    if (verbose) {

      message("             finished after ", round(Sys.time() - start_time, 1), " ",
              units(Sys.time() - start_time), ".")
      rm(start_time)

    }

  }

  ## Process scenarios

  ## message to user
  if (verbose) {

    start_time  <- Sys.time()
    message("         Processing Scenarios ...")

  }

  ## prepare foreach
  exported_stuff     <- c("scenario_list", "trials_unique", "go_trial_indices",
                          "applicable_previous_trials", "method_names", "method_quantiles_list",
                          "prior_parameters_list", "quantiles", "n_mcmc_iterations",
                          "getRowIndexOfVectorInMatrix")

  "%do%"             <- foreach::"%do%"
  "%dopar%"          <- foreach::"%dopar%"

  if (identical(length(scenario_numbers), 1L)) {
    n_cores <- 1L
  }

  foreach_expression <- quote({

    ## Find the indices of the trials of a specific scenario
    scenario_data_matrix <- cbind(scenario_list[[k]]$n_responders,
                                  scenario_list[[k]]$n_subjects)

    trial_indices <- unlist(sapply(seq_len(nrow(scenario_data_matrix)), function (i) {
      getRowIndexOfVectorInMatrix(
        vector_to_be_found    = scenario_data_matrix[i, ],
        matrix_to_be_searched = trials_unique)
    }))

    map_vector                   <- numeric(length(go_trial_indices))
    map_vector[go_trial_indices] <- seq_along(which(go_trial_indices))
    ## note: if applicable_previous_trials == FALSE, then map_indices == trial_indices
    map_indices                  <- map_vector[trial_indices]

    ## check whether there where previous analysis
    if (applicable_previous_trials) {

      scenario_method_quantiles <- scenario_list[[k]]$previous_analyses$post_quantiles

    } else {

      scenario_method_quantiles <- vector(mode = "list", length = length(method_names))
      names(scenario_method_quantiles) <- method_names

    }

    ## Save scenario specific posterior quantiles for each method
    for (n in seq_along(method_names)) {

      if (applicable_previous_trials) {

        ## This loop assigns only the updated/new trial realizations
        ## lapply() implementation slower than for-loop, probably due to re-assignments
        for (i in seq_along(map_indices)) {

          ## replace only updated trials
          if (!isTRUE(all.equal(0, map_indices[[i]]))) {

            scenario_method_quantiles[[method_names[n]]][[i]] <-
              method_quantiles_list[[method_names[n]]][[map_indices[[i]]]]

          }

        }

      } else {

        ## lapply() implementation that is not adapted to the inclusion of previously
        ## analyzed trial realizations, but much faster than one that was
        scenario_method_quantiles[[method_names[n]]] <- lapply(trial_indices, function (i) {
          method_quantiles_list[[method_names[n]]][[i]]
        })

      }

    }

    ## combine results
    scenario_analysis_list <- list(
      quantiles_list      = scenario_method_quantiles,
      scenario_data       = scenario_list[[k]],
      analysis_parameters = list(
        quantiles             = quantiles,
        method_names          = method_names,
        prior_parameters_list = prior_parameters_list,
        n_mcmc_iterations     = n_mcmc_iterations))

    return (scenario_analysis_list)

  })

  ## run foreach
  if (n_cores > 1) {
    doParallel::registerDoParallel(n_cores)
    on.exit(doParallel::stopImplicitCluster())

    analyses_list <- foreach::foreach(k = seq_along(scenario_numbers),
                                      .verbose  = FALSE,
                                      .export   = exported_stuff
    ) %dopar% {eval(foreach_expression)}

  } else {

    analyses_list <- foreach::foreach(k = seq_along(scenario_numbers),
                                      .verbose  = FALSE,
                                      .export   = exported_stuff
    ) %do% {eval(foreach_expression)}

  }

  ## message to user
  if (verbose) {

    message("             finished after ", round(Sys.time() - start_time, 1), " ",
            units(Sys.time() - start_time), ".")
    rm(start_time)

  }

  names(analyses_list) <- paste0("scenario_", scenario_numbers)
  class(analyses_list) <- "analysis_list"

  return (analyses_list)

}

is.analysis_list <- function (x) {

  if (missing(x)) stop ("Please provide an object for the argument 'x'")

  inherits(x, "analysis_list")

}

makeUniqueRows <- function (

  matrix

) {

  n_rows      <- nrow(matrix)
  n_cols      <- ncol(matrix)

  unique_rows <- stats::aggregate(formula = id ~ .,
                                  data    = cbind(id = seq_along(n_rows), matrix),
                                  FUN     = length)

  return (unique_rows[, seq_len(n_cols)])

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
#'     n_mcmc_iterations  = 100,
#'     n_cores            = 1L)
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

  error_analyses_list <- simpleError(
    "Please provide an object of class analysis_list for the argument 'analyses_list'")
  error_save_path     <- simpleError(
    "Please provide a string containing a path for the argument 'save_path'")
  error_analysis_numbers <- simpleError(paste(
    "Please provide a vector of positive integers for the argument 'analysis_numbers'",
    "with length equal to the length of 'analyses_list'"))

  if (missing(analyses_list))                                     stop (error_analyses_list)

  if (!is.analysis_list(analyses_list))                           stop (error_analyses_list)
  if (!is.character(save_path) || length(save_path) > 1)          stop (error_save_path)
  if (!is.null(analysis_numbers) && (
    any(!is.positive.wholenumber(analysis_numbers)) ||
    !identical(length(analyses_list), length(analysis_numbers)))) stop(error_analysis_numbers)

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

#' @title loadAnalyses
#' @md
#' @description This function loads an analysis performed with
#' \code{\link[bhmbasket]{performAnalyses}}
#' @param load_path A string providing a path where the scenarios are being stored,
#' Default: \code{\link[base]{tempfile}}
#' @param scenario_numbers A (vector of) positive integer(s) for the scenario number(s)
#' @param analysis_numbers A (vector of) positive integer(s) for the analysis number(s)
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
#'     n_mcmc_iterations  = 100,
#'     n_cores            = 1L)
#'
#'   save_info     <- saveAnalyses(analysis_list)
#'   analysis_list <- loadAnalyses(scenario_numbers = save_info$scenario_numbers,
#'                                 analysis_numbers = save_info$analysis_numbers,
#'                                 load_path        = save_info$path)
#' @author Stephan Wojciekowski
#' @export
loadAnalyses <- function (

  scenario_numbers,
  analysis_numbers,
  load_path = tempdir()

) {

  error_scenario_numbers <- simpleError(
    "Please provide a vector of positive integers for the argument 'scenario_numbers'")
  error_analysis_numbers <- simpleError(
    "Please provide a vector of positive integers for the argument 'analysis_numbers'")
  error_load_path        <- simpleError(
    "Please provide a string containing a path for the argument 'load_path'")

  if (missing(scenario_numbers)) stop (error_scenario_numbers)
  if (missing(analysis_numbers)) stop (error_analysis_numbers)

  if (!is.character(load_path) || length(load_path) > 1) stop (error_load_path)

  if (!identical(length(scenario_numbers), length(analysis_numbers))) {
    stop (simpleError("'scenario_numbers' and 'analysis_numbers' must have equal length"))
  }

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



