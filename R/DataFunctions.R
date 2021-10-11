# library(sinew)
# makeOxygen(getUniqueTrials)

#' @title createTrial
#' @description This function creates an object of class `scenario_list`
#' for a single trial outcome, which can subsequently be analyzed with other functions of
#' `bhmbasket`, e.g. \code{\link[bhmbasket]{performAnalyses}}
#' @param n_subjects A vector of integers for the number of subjects in the trial outcome
#' @param n_responders A vector of integers for the number of responders in the trial outcome
#' @return An object of class `scenario_list` with the scenario data for a single trial outcome.
#' @details This function is a wrapper for \code{\link[bhmbasket]{simulateScenarios}} with
#' ```
#' simulateScenarios(
#'  n_subjects_list     = list(n_subjects),
#'  response_rates_list = list(n_responders),
#'  n_trials            = 1)
#' ```
#' @seealso
#'  \code{\link[bhmbasket]{simulateScenarios}}
#'  \code{\link[bhmbasket]{performAnalyses}}
#' @author Stephan Wojciekowski
#' @rdname createTrial
#' @examples
#'  trial_outcome <- createTrial(n_subjects   = c(10, 20, 30, 40),
#'                               n_responders = c( 1,  2,  3,  4))
#' @export
#' @md
createTrial <- function (

  n_subjects,
  n_responders

) {

  error_n_subjects <- simpleError(
    "Please provide a vector of integers for the argument 'error_n_subjects'")
  error_n_responders <- simpleError(
    "Please provide a vector of integers for the argument 'n_responders'")

  if (missing(n_subjects))               stop (error_n_subjects)
  if (missing(n_responders))             stop (error_n_responders)

  if (any(!is.wholenumber(n_subjects)))   stop (error_n_subjects)
  if (any(!is.wholenumber(n_responders))) stop (error_n_responders)

  utils::capture.output({

    trial <- simulateScenarios(
      n_subjects_list     = list(n_subjects),
      response_rates_list = list(n_responders),
      n_trials            = 1)

  })

  return (trial)

}

getScenario <- function (

  n_subjects,
  response_rates,

  cohort_names = seq_along(n_subjects),

  n_trials     = 1e4

) {

  # n_cores      = parallel::detectCores() - 1L # old argument that is no longer needed
  # "%dopar%" <- foreach::"%dopar%"

  response_rates           <- convertVector2Matrix(response_rates)
  colnames(response_rates) <- paste0("rr_", cohort_names)

  if (any(response_rates < 1 & response_rates > 0)) {

    ## Simulations for new cohorts
    ## A cohort is new if it has a response rate greater than 0 and less than 1.
    ## Response rates equal to 0 or
    ## greater than or equal to 1 will be used as fixed responses.

    new_cohorts <- TRUE
    index_new   <- which(response_rates < 1 & response_rates > 0)

    if (length(n_subjects[index_new]) != length(response_rates[, index_new])) {
      stop ("n_subjects and response_rates must have same length")
    }

    ## Simulate the number of responses without interim analysis
    n_responders <- getRespondersNonParallel(
      response_rates = response_rates[, index_new],
      n_subjects     = n_subjects[index_new],
      n_trials       = n_trials)
      # n_cores        = n_cores,
      # seed           = seed)

  } else {

    ## No new cohorts, as response rates are all greater than or equal to 1

    new_cohorts <- FALSE

  }

  if (any(response_rates >= 1 | response_rates == 0)) {

    ## Historical cohorts

    hist_cohorts <- TRUE
    index_hist   <- which(response_rates >= 1 | response_rates == 0)

    if (length(n_subjects[index_hist]) != length(response_rates[, index_hist])) {
      stop ("n_subjects and response_rates must have same length")
    }

    if (new_cohorts) {

      ## In case of simulated (new) cohorts, the historic responses will be recycled to match
      ## the number of unique trials

      ## Make matrix if necessary
      n_responders_hist <-  matrix(rep(response_rates[index_hist],
                                       each = nrow(n_responders)),
                                   nrow = nrow(n_responders))

    } else {

      ## In case of no simulated cohorts, only one set of fixed responses is needed

      n_responders_hist <- matrix(response_rates[index_hist], nrow = 1)

    }

  } else {

    ## No historical cohorts, as response rates are all smaller than 1

    hist_cohorts <- FALSE

  }

  ## Combine new and historical cohorts as appropriate
  if (new_cohorts & hist_cohorts) {

    n_responders <- cbind(n_responders, n_responders_hist)

  } else if (hist_cohorts) {

    n_responders <- n_responders_hist

  }

  previous_gos <- matrix(TRUE, byrow = TRUE,
                         ncol = length(n_subjects) + 1L,
                         nrow = nrow(n_responders))
  colnames(previous_gos) <- c("overall", paste0("decision_", cohort_names))

  n_subjects <- matrix(n_subjects, byrow = TRUE,
                       ncol = length(n_subjects),
                       nrow = nrow(n_responders))

  colnames(n_subjects)   <- paste0("n_", cohort_names)
  colnames(n_responders) <- paste0("r_", cohort_names)


  ## Create list to return data
  scenario_data <- list(n_subjects        = n_subjects,
                        n_responders      = n_responders,
                        response_rates    = response_rates,
                        previous_analyses = list(go_decisions   = previous_gos,
                                                 post_quantiles = NULL),
                        n_trials          = n_trials)
                        # seed              = seed)

  # class(scenario_data) <- "scenario_data"

  return (scenario_data)

}

is.scenario_data <- function (x) {

  if (missing(x)) stop ("Please provide an object for the argument 'x'")

  inherits(x, "decision_list")

  # ## x must be list
  # if (!is.list(x)) {
  #   warning ("x must be list")
  #   return (FALSE)
  # }
  #
  # ## names of x must be as below
  # if (!identical(sort(names(x)), c("n_responders", "n_subjects", "n_trials",
  #                                  "response_rates", "scenario_number", "seed"))) {
  #   warning ("names of x must be c('n_responders', 'n_subjects', 'n_trials',
  #                                  'response_rates', 'scenario_number', 'seed')")
  #   return (FALSE)
  # }
  #
  # ## x$n_subjects must be matrix with non-negative integers
  # if (!is.matrix(x$n_subjects) || !is.numeric(x$n_subjects) ||
  #     !all(is.wholenumber(x$n_subjects)) || !all(x$n_subjects >= 0)) {
  #   warning ("x$n_subjects must be matrix with non-negative integers")
  #   return (FALSE)
  # }
  #
  # ## x$n_responders must be matrix with non-negative integers
  # if (!is.matrix(x$n_responders) || !is.numeric(x$n_responders) ||
  #     !all(is.wholenumber(x$n_responders)) || !all(x$n_responders >= 0)) {
  #   warning ("x$n_responders must be matrix with non-negative integers")
  #   return (FALSE)
  # }
  #
  # ## x$n_subjects and x$n_responders must have same dimension
  # if (!identical(dim(x$n_responders), dim(x$n_subjects))) {
  #   warning ("x$n_subjects and x$n_responders must have same dimension")
  #   return (FALSE)
  # }
  #
  # ## x$response_rates must be matrix with with integers in [0, 1]
  # if (!is.matrix(x$response_rates) || !is.numeric(x$response_rates) ||
  #     !all(x$response_rates >= 0 && x$response_rates <= 1)) {
  #   warning ("x$response_rates must be matrix with with integers in [0, 1]")
  #   return (FALSE)
  # }
  #
  # ## x$response_rates must have same number of columns as x$n_subjects
  # if (!identical(ncol(x$response_rates), ncol(x$n_subjects))) {
  #   warning ("x$response_rates must have same number of columns as x$n_subjects")
  #   return (FALSE)
  # }
  #
  # ## x$response_rates must have number of rows between 1 and nrow(x$n_subjects)
  # if (!nrow(x$response_rates) >= 1 || !nrow(x$response_rates) <= nrow(x$n_subjects)) {
  #   warning ("x$response_rates must have number of rows between 1 and nrow(x$n_subjects)")
  #   return (FALSE)
  # }
  #
  # ## x$n_trials must be positive integer
  # if (!is.numeric(x$n_trials) || !is.wholenumber(x$n_trials) || !x$n_trials > 0) {
  #   warning ("x$n_trials must be positive integer")
  #   return (FALSE)
  # }
  #
  # ## x$n_trials must be equal to nrow(x$n_subjects)
  # if (!isTRUE(all.equal(x$n_trials, nrow(x$n_subjects)))) {
  #   warning ("x$n_trials must be equal to nrow(x$n_subjects)")
  #   return (FALSE)
  # }
  #
  # ## x$seed must be non-negative integer
  # if (!is.numeric(x$seed) || !is.wholenumber(x$seed) || !x$seed >= 0) {
  #   warning ("x$seed must be non-negative integer")
  #   return (FALSE)
  # }
  #
  # ## x$scenario_number must be non-negative integer
  # if (!is.numeric(x$scenario_number) || !is.wholenumber(x$scenario_number) ||
  #     !x$scenario_number >= 0) {
  #   warning ("x$scenario_number must be non-negative integer")
  #   return (FALSE)
  # }
  #
  # return (TRUE)

}

getResponders <- function (

  n_subjects,
  response_rates

) {

  ## Adjust for working with apply()
  if (is.null(dim(response_rates))) {
    response_rates <- t(as.matrix(response_rates))
  }

  n_responders <- t(apply(response_rates, 1,
                          function (x) {stats::rbinom(n    = length(n_subjects),
                                                      size = n_subjects,
                                                      prob = x)}))

  return (n_responders)

}

#' @title simulateScenarios
#' @description This function creates scenarios for the analysis with
#' \code{\link[bhmbasket]{performAnalyses}}.
#' @param n_subjects_list A list that contains for each scenario a vector for
#' the number of subjects per cohort.
#' A single vector can be provided if all scenarios should have the same number of subjects.
#' @param response_rates_list A list that contains for each scenario a vector for
#' the response rates per cohort.
#' @param scenario_numbers A vector of positive integers naming the scenarios,
#' Default: `seq_along(response_rates_list)`
#' @param n_trials An integer indicating the number of trial simulations per response rates,
#' Default: `10000`. If `n_trials` is present in `.GlobalEnv` and `missing(n_trials)`,
#' the globally available value will be used.
#' @return An object of class `scenario_list` with the scenario data for each specified scenario.
#' @details The function simulates trials with binary outcome for each scenario.
#' Integer values for the response rates will be treated as observed outcomes.
#' @author Stephan Wojciekowski
#' @seealso
#'  \code{\link[bhmbasket]{saveScenarios}}
#'  \code{\link[bhmbasket]{createTrial}}
#'  \code{\link[bhmbasket]{performAnalyses}}
#' @rdname simulateScenarios
#' @examples
#'   n_subjects     <- c(10, 20, 30)
#'
#'   rr_negative    <- rep(0.1, 3)
#'   rr_nugget      <- c(0.9, 0.1, 0.1)
#'   rr_positive    <- rep(0.9, 3)
#'
#'   scenarios_list <- simulateScenarios(
#'     n_subjects_list     = list(n_subjects,
#'                                n_subjects,
#'                                n_subjects),
#'     response_rates_list = list(rr_negative,
#'                                rr_nugget,
#'                                rr_positive))
#' @export
#' @md
simulateScenarios <- function (

  n_subjects_list,
  response_rates_list,
  scenario_numbers = seq_along(response_rates_list),

  n_trials         = 1e4

) {

  error_n_subjects_list     <- simpleError(
    "Please provide a list of vectors of positive integers for the argument 'n_subjects_list'")
  error_response_rates_list <- simpleError(
    paste("Please provide a list of vectors of non-negative numerics for the argument",
          "'response_rates_list'\n", "Values outside of (0, 1) must be integers"))
  error_n_trials            <- simpleError(
    "Please provide a positive integer for the argument 'n_trials'")
  error_scenario_numbers    <- simpleError(
    "Please provide a vector of positive integers for the argument 'scenario_numbers'")

  if (missing(n_subjects_list))                                stop (error_n_subjects_list)
  if (missing(response_rates_list))                            stop (error_response_rates_list)

  if (!is.list(response_rates_list) ||
      any(!sapply(response_rates_list, is.numeric)))           stop (error_response_rates_list)

  ## put n_subjects as a list if provided as vector
  if (!is.list(n_subjects_list)) {
    n_subjects_list <- rep(list(n_subjects_list), length(response_rates_list))
  }

  if (!is.list(n_subjects_list) ||
      any(!sapply(n_subjects_list, is.positive.wholenumber)))  stop (error_n_subjects_list)

  if (!is.positive.wholenumber(scenario_numbers))              stop (error_scenario_numbers)

  if (!identical(length(scenario_numbers), length(n_subjects_list))) stop (simpleError(
    "'scenario_numbers' and 'n_subjects_list' must have same lenth"
  ))

  if (!identical(length(n_subjects_list), length(response_rates_list))) stop (simpleError(
    "'n_subjects_list' and 'response_rates_list' must have same length"))

  if (any(!sapply(response_rates_list, function (x) {
    identical(length(response_rates_list[[1]]), length(x))}))) stop (simpleError(
      "All scenarios within a set of scenarios must have the same number of cohorts"))

  if (identical(length(response_rates_list[[1]]), 1L)) stop (simpleError(
      "Each scenario must have at least 2 cohorts"))

  if (any(!sapply(n_subjects_list, function (x) {
    identical(length(response_rates_list[[1]]), length(x))}))) stop (simpleError(
      "All scenarios within a set of scenarios must have the same number of cohorts"))

  if (any(sapply(seq_along(n_subjects_list), function (x) {
    any(n_subjects_list[[x]] < response_rates_list[[x]])})))   stop (simpleError(
      "Values in 'response_rates_list' must not be greater than values in 'n_subjects_list'"))

  if (any(!sapply(response_rates_list, function (x) {
    is.numeric(x) & x > 0 & x < 1 |
      is.wholenumber(x) & (x == 0 | x >= 1)})))                stop(error_response_rates_list)

  ## check whether n_trials is present in global environment
  if ("n_trials" %in% ls(envir = .GlobalEnv) & missing(n_trials)) {
    n_trials <- get("n_trials", envir = .GlobalEnv)
  }

  if (!is.single.positive.wholenumber(n_trials))               stop (error_n_trials)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  # cat(format(Sys.time(), "%m/%d/%y"), " Simulating Scenarios\n", sep = "")

  scenario_list    <- vector(mode = "list", length = length(scenario_numbers))
  for (s in seq_along(scenario_numbers)) {

    # start_time  <- Sys.time()
    # out_message <- paste0(format(start_time, "   %H:%M"),
    #                       " Simulating Scenario ", scenario_numbers[s], " ...")
    # cat(out_message, rep(".", 33 + nchar(max(scenario_numbers)) - nchar(out_message)), sep = "")

    scenario_list[[s]] <- getScenario(
      n_subjects          = n_subjects_list[[s]],
      response_rates      = response_rates_list[[s]],
      n_trials            = n_trials)

    scenario_list[[s]]$scenario_number <- scenario_numbers[s]

    # cat(" Finished after ", round(Sys.time() - start_time, 1), " ",
    #     units(Sys.time() - start_time), ".\n", sep = "")
    # rm(start_time)

  }

  names(scenario_list) <- paste0("scenario_", scenario_numbers)
  class(scenario_list) <- "scenario_list"

  return (scenario_list)

}

is.scenario_list <- function (x) {

  if (missing(x)) stop ("Please provide an object for the argument 'x'")

  inherits(x, "scenario_list")

  # ## x must be list
  # if (!is.list(x)) {
  #   warning ("x must be list")
  #   return (FALSE)
  # }
  #
  # is_scenario_data <- sapply(x, is.scenario_data)
  #
  # if (any(!is_scenario_data)) {
  #   warning (paste0("The items of x with the following indices are not scenario_data", !which(is_scenario_data)))
  #   return (FALSE)
  # }
  #
  # return (TRUE)

}

#' @title saveScenarios
#' @md
#' @description Saves the scenario data in a newly created or existing directory
#' @param scenario_list An object of class `scenario_list`, e.g. created with `simulateScenarios()`
#' @param save_path A string providing the path for the directory in which the directory of the
#' scenario should be created, Default: \code{\link[base]{tempfile}}
#' @return A named list of length 2 with the scenario numbers and the `save_path`
#' @author Stephan Wojciekowski
#' @seealso
#'  \code{\link[bhmbasket]{simulateScenarios}}
#'  \code{\link[bhmbasket]{loadScenarios}}
#'  \code{\link[base]{tempfile}}
#' @examples
#'   scenarios_list <- simulateScenarios(
#'     n_subjects_list     = list(c(10, 20, 30)),
#'     response_rates_list = list(rep(0.9, 3)),
#'     n_trials            = 10)
#'
#'   save_info      <- saveScenarios(scenarios_list)
#'   scenarios_list <- loadScenarios(scenario_numbers = save_info$scenario_numbers,
#'                                   load_path        = save_info$path)
#' @rdname saveScenarios
#' @export
saveScenarios <- function (

  scenario_list,
  save_path = tempdir()

) {

  if (missing(scenario_list))
    stop ("Please provide an object of class scenario_list for the argument 'scenario_list'")

  if (!is.scenario_list(scenario_list))
    stop ("Please provide an object of class scenario_list for the argument 'scenario_list'")
  if (!is.character(save_path) || length(save_path) > 1)
    stop ("Please provide a string for the argument 'save_path'")

  if (!dir.exists(save_path)) {
    dir.create(save_path)
  }

  scenario_numbers <- sapply(scenario_list, function (x) x$scenario_number)

  for (s in seq_along(scenario_list)) {

    saveRDS(scenario_list[[s]],
            file = paste0(save_path, "/scenario_data_", scenario_numbers[s], ".rds"))

  }

  return (list(scenario_numbers = scenario_numbers, path = save_path))

}

#' @title loadScenarios
#' @md
#' @description This function loads scenarios saved with \code{\link[bhmbasket]{saveScenarios}}
#' @param scenario_numbers A vector of integers naming the scenario to be loaded
#' @param load_path A string for the directory where the scenarios are being stored,
#' Default: \code{\link[base]{tempfile}}
#' @return Returns an object of class `scenario_list`
#' @rdname loadScenarios
#' @author Stephan Wojciekowski
#' @seealso
#'  \code{\link[bhmbasket]{simulateScenarios}}
#'  \code{\link[bhmbasket]{saveScenarios}}
#'  \code{\link[base]{tempfile}}
#' @examples
#'   scenarios_list <- simulateScenarios(
#'     n_subjects_list     = list(c(10, 20, 30)),
#'     response_rates_list = list(rep(0.9, 3)),
#'     n_trials            = 10)
#'
#'   save_info      <- saveScenarios(scenarios_list)
#'   scenarios_list <- loadScenarios(scenario_numbers = save_info$scenario_numbers,
#'                                   load_path        = save_info$path)
#' @export
loadScenarios <- function (

  scenario_numbers,
  load_path = tempdir()

) {

  if (missing(scenario_numbers))
    stop ("Please provide a vector of positive integers for the argument 'scenario_numbers'")

  if (!is.numeric(scenario_numbers) || !is.wholenumber(scenario_numbers) || scenario_numbers <= 0)
    stop ("Please provide a positive integer for the argument 'scenario_numbers'")
  if (!is.character(load_path) || length(load_path) > 1)
    stop ("Please provide a string for the argument 'load_path'")

  scenario_list <- vector(mode = "list", length = length(scenario_numbers))

  for (s in seq_along(scenario_numbers)) {

    paste0("scenario_data_", scenario_numbers[s], ".rds")

    scenario_list[[s]] <- readRDS(
      file.path(load_path, paste0("scenario_data_", scenario_numbers[s], ".rds")))

    names(scenario_list) <- paste0("scenario_", scenario_numbers)

  }

  class(scenario_list) <- "scenario_list"

  return (scenario_list)

}

getInterimPassed <- function (

  interim_n_responders,
  interim_n_min

) {

  interim_passed <- t(apply(interim_n_responders, 1, function (x) {
    x >= interim_n_min
  }))

  return (interim_passed)

}

getRespondersNonParallel <- function (

  response_rates,
  n_subjects,

  n_trials

) {

  n_responders <- t(replicate(n = n_trials, {

    stats::rbinom(n    = length(n_subjects), # n_cohorts
                  size = n_subjects,
                  prob = response_rates)



  }))

  return (n_responders)

}


getRespondersParallel <- function (

  response_rates,
  n_subjects,

  n_trials,
  n_cores = parallel::detectCores() - 1L
  # seed    = as.numeric(Sys.time())

) {

  "%dopar%" <- foreach::"%dopar%"

  ## somehow %dopar% does not default to %do% when no connection open
  doParallel::registerDoParallel(n_cores)
  on.exit(doParallel::stopImplicitCluster())
  # if (n_cores > 1) {
  #   doParallel::registerDoParallel(n_cores)
  #   on.exit(doParallel::stopImplicitCluster())
  # }
  n_responders <- foreach::foreach(k = seq_len(n_trials),
                                   .combine = rbind,
                                   .export = c("getResponders")) %dopar% {
                                     getResponders(response_rates = response_rates,
                                                   n_subjects     = n_subjects)}
                                                   # seed           = seed + k)}

  return (n_responders)

}

#' @title continueRecruitment
#' @md
#' @description This function continues the recruitment of subjects for a set of scenarios
#' based on the Go / NoGo decisions in the simulated trial outcomes of said scenarios.
#' @param n_subjects_add_list A list that contains for each scenario an integer vector for
#' the number of subjects per cohort to be additionally recruited.
#' @param decisions_list A list with decisions per scenario created with
#' \code{\link[bhmbasket]{getGoDecisions}}
#' @param method_name A string for the method name of the analysis the decisions are based on
#' @return An object of class `scenario_list` with the scenario data for each specified scenario.
#' @details
#' This function is intended to be used for analyses with the following work flow:\cr
#' `simulateScenarios()` -> `performAnalyses()` -> `getGoDecisions()`-> \cr
#' `continueRecruitment()` -> `performAnalyses()` -> `getGoDecisions()`-> \cr
#' `continueRecruitment()` -> ...
#'
#' Note that `n_subjects_add_list` takes the additional number of subjects to be recruited,
#' not the overall number of subjects.
#' This way the work flow can be repeated as often as
#' required, which can be useful e.g. for interim analyses.
#' @examples
#' interim_scenarios <- simulateScenarios(
#'   n_subjects_list     = list(c(10, 20, 30)),
#'   response_rates_list = list(rep(0.9, 3)),
#'   n_trials            = 10)
#'
#' interim_analyses <- performAnalyses(
#'   scenario_list       = interim_scenarios,
#'   target_rates        = rep(0.5, 3),
#'   n_mcmc_iterations   = 100,
#'   n_cores             = 1L)
#'
#' interim_gos <- getGoDecisions(
#'   analyses_list       = interim_analyses,
#'   cohort_names        = c("p_1", "p_2", "p_3"),
#'   evidence_levels     = c(0.5, 0.8, 0.5),
#'   boundary_rules      = quote(c(x[1] > 0.8, x[2] > 0.6, x[3] > 0.7)))
#'
#' scenarios_list <- continueRecruitment(
#'   n_subjects_add_list = list(c(30, 20, 10)),
#'   decisions_list      = interim_gos,
#'   method_name         = "exnex_adj")
#' @seealso
#'  \code{\link[bhmbasket]{simulateScenarios}}
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getGoDecisions}}
#' @rdname continueRecruitment
#' @author Stephan Wojciekowski
#' @export
continueRecruitment <- function (

  n_subjects_add_list,
  decisions_list,

  method_name

) {

  error_n_subjects_add_list <- simpleError(
    "Please provide a list of vectors of positive integers for the argument 'n_subjects_add_list'")
  error_decisions_list <- simpleError(
    "Please provide an object of class decision_list for the argument 'decisions_list'")
  error_method_name <- simpleError(paste(
    "Please provide a string naming an analysis method for the argument 'method_name'",
    "Must be one of 'berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'"))

  if (missing(n_subjects_add_list))            stop (error_n_subjects_add_list)
  if (missing(decisions_list))                 stop (error_decisions_list)
  if (missing(method_name))                    stop (error_method_name)

  method_name <- tryCatch({

    match.arg(
      method_name,
      choices    = c('berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'),
      several.ok = FALSE)

  }, error = function (e) e)

  if (!is.decision_list(decisions_list))       stop (error_decisions_list)
  if (inherits(method_name, "error"))          stop (error_method_name)

  if (!is.list(n_subjects_add_list)) {
    n_subjects_add_list <- rep(list(n_subjects_add_list), length(decisions_list))
  }

  if (!is.list(n_subjects_add_list) ||
      any(!sapply(n_subjects_add_list, is.non.negative.wholenumber)))
                                               stop (error_n_subjects_add_list)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ## Note: also some checks hereafter

  ## get scenario numbers
  scenario_numbers <- as.numeric(sub("scenario_", "", names(decisions_list)))

  if (length(n_subjects_add_list) != length(decisions_list)) {
    stop (simpleError("The lengths of 'n_subjects_add_list' and 'decisions_list' must be equal"))
  }

  # cat(format(Sys.time(), "%m/%d/%y"), " Continuing Recruitment\n", sep = "")

  scenario_list <- vector(mode = "list", length = length(decisions_list))
  names(scenario_list) <- paste0("scenario_", scenario_numbers)
  for (s in seq_along(scenario_list)) {

    if (!(method_name %in%
          decisions_list[[s]]$analysis_data$analysis_parameters$method_names)) {

      stop (simpleError("Selected method_name not analyzed"))

    }

    # start_time  <- Sys.time()
    # out_message <- paste0(format(start_time, "   %H:%M"),
    #                       " Simulating Scenario ", scenario_numbers[s], " ...")
    # cat(out_message, rep(".", 33 + nchar(max(scenario_numbers)) - nchar(out_message)), sep = "")

    ## get new data, i.e. get new responders and new number of subjects per trial

    n_subjects_add <- n_subjects_add_list[[s]]
    response_rates <- decisions_list[[s]]$scenario_data$response_rates
    cohort_names   <- sub("rr_", "", colnames(response_rates))

    if (any(response_rates > 0 & response_rates < 1)) {

      index_new <- which(response_rates > 0 & response_rates < 1)

    } else {

      stop (simpleError(paste0(
        "Only historical cohorts in scenario ",
        decisions_list[[s]]$scenario_data$scenario_number)))

    }

    response_rates_new <- response_rates[, index_new]
    cohort_names_new   <- cohort_names[index_new]

    if (!identical(length(n_subjects_add), length(response_rates_new))) {
      stop (simpleError(paste0(
        "The length of n_subjects_add must be equal ",
        "to the length of the response rates that are in (0, 1)")))
    }

    n_trials <- decisions_list[[s]]$scenario_data$n_trials
    # seed     <- decisions_list[[s]]$scenario_data$seed + n_trials + 1

    add_scenario <- getScenario(
      n_subjects     = n_subjects_add,
      response_rates = response_rates_new,
      cohort_names   = cohort_names_new,
      n_trials       = n_trials)
      # seed           = seed)

    n_responders_add <- add_scenario$n_responders
    n_subjects_add   <- add_scenario$n_subjects

    ## Combine with existing data

    go_decisions <- decisions_list[[s]]$decisions_list[[method_name]]
    previous_gos <- go_decisions

    if ("overall" %in% colnames(go_decisions)) {
      go_decisions <- go_decisions[, -which(colnames(go_decisions) == "overall")]
    }
    if (!all(index_new %in% as.numeric(sub("decision_", "", colnames(go_decisions))))) {
      stop (simpleError(
        "There must be a decision for each recruiting cohort in the 'decisions_list'"))
    }

    go_decisions <- go_decisions[, index_new]

    n_responders <- decisions_list[[s]]$scenario_data$n_responders
    n_subjects   <- decisions_list[[s]]$scenario_data$n_subjects

    n_responders[, index_new] <- n_responders[, index_new] + go_decisions * n_responders_add
    n_subjects[, index_new]   <- n_subjects[, index_new] + go_decisions * n_subjects_add

    ## Saving Scenario

    scenario_list[[s]] <- list(
      n_subjects        = n_subjects,
      n_responders      = n_responders,
      response_rates    = response_rates,
      previous_analyses = list(go_decisions   = previous_gos,
                               post_quantiles = decisions_list[[s]]$analysis_data$quantiles_list),
      n_trials          = n_trials)
      # seed              = seed)

    scenario_list[[s]]$scenario_number <-
      decisions_list[[s]]$scenario_data$scenario_number

    # cat(" Finished after ", round(Sys.time() - start_time, 1), " ",
    #     units(Sys.time() - start_time), ".\n", sep = "")
    # rm(start_time)

  }

  class(scenario_list) <- "scenario_list"

  return (scenario_list)

}

getNSubjects <- function (

  recruitment_per_month,
  start_date,
  analysis_dates,

  date_format = "%m/%d/%Y"

) {

  recruitment_per_day <- recruitment_per_month * (12 / 365)

  start_date     <- as.Date(start_date, format = date_format)
  analysis_dates <- as.Date(analysis_dates, format = date_format)

  n_subjects_matrix <- floor((analysis_dates - start_date) %o% recruitment_per_day)
  rownames(n_subjects_matrix) <- as.character(analysis_dates, format = date_format)
  colnames(n_subjects_matrix) <- paste0("cohort_", seq_len(ncol(n_subjects_matrix)))

  return (n_subjects_matrix)

}

getRecruitment <- function (

  n_subjects_required,
  recruitment_per_month,
  start_date,

  date_format = "%m/%d/%Y"

) {

  if (missing(n_subjects_required))
    stop ("Please provide a matrix of non-negative integers for the argument 'n_subjects_required'")
  if (missing(recruitment_per_month))
    stop ("Please provide a vector of non-negative numerics for the argument 'recruitment_per_month'")
  if (missing(start_date))
    stop ("Please provide a string in the format 'date_format' of for the argument 'start_date'")

  if (!is.numeric(n_subjects_required) || #!is.matrix(n_subjects_required) ||
      any(!is.wholenumber(n_subjects_required)) || any(n_subjects_required < 0))
    stop ("Please provide a matrix of non-negative integers for the argument 'n_subjects_required'")
  if (!is.numeric(recruitment_per_month) || any(recruitment_per_month < 0))
    stop ("Please provide a vector of non-negative numerics for the argument 'recruitment_per_month'")

  ## Get indices of historical cohorts
  hist_index <- recruitment_per_month == 0

  recruitment_per_day <- recruitment_per_month[!hist_index] * (12 / 365)

  start_date <- as.Date(start_date, format = date_format)

  n_subjects_required <- convertVector2Matrix(n_subjects_required)

  if (!identical(length(recruitment_per_month), ncol(n_subjects_required)))
    stop (paste0("The number of columns of 'n_subjects_required' ",
                 "and the length of 'recruitment_per_month' must be equal"))

  days_required <- apply(convertVector2Matrix(n_subjects_required[, !hist_index]),
                         1, function (n_subj) {
    max(ceiling(n_subj * recruitment_per_day^-1))
  })

  n_subjects_matrix_0 <- getNSubjects(
    recruitment_per_month = recruitment_per_month[!hist_index],
    start_date            = start_date,
    analysis_dates        = start_date + days_required,
    date_format           = date_format)

  n_subjects_matrix <- matrix(NA,
                              ncol = ncol(n_subjects_required),
                              nrow = nrow(n_subjects_required))

  n_subjects_matrix[, !hist_index] <- n_subjects_matrix_0
  n_subjects_matrix[, hist_index]  <- n_subjects_required[, hist_index]

  rownames(n_subjects_matrix) <- rownames(n_subjects_matrix_0)
  colnames(n_subjects_matrix) <- paste0("cohort_", seq_len(ncol(n_subjects_required)))

  return (n_subjects_matrix)

}
