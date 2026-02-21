#' @title continueRecruitment
#' @md
#' @description This function continues the recruitment of subjects for a set of scenarios
#' based on the Go / NoGo decisions in the simulated trial outcomes of said scenarios.
#' @param n_subjects_add_list A list that contains for each scenario an integer vector for
#' the number of subjects per cohort to be additionally recruited.
#' @param decisions_list A list with decisions per scenario created with
#' \code{\link[bhmbasket]{getGoDecisions}}
#' @param method_name A string for the method name of the analysis the decisions are based on.
#' Can be `NULL` if only one method has been used for analysis, Default: `NULL`
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
#'   n_mcmc_iterations   = 100)
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
  
  method_name = NULL
  
) {
  
  error_n_subjects_add_list <- 
    "Providing a list of vectors of positive integers for the argument 'n_subjects_add_list'"
  error_decisions_list <- 
    "Providing an object of class decision_list for the argument 'decisions_list'"
  error_method_name <- simpleError(paste(
    "Please provide a string naming an analysis method for the argument 'method_name'",
    "Must be one of 'berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'"))
  
  checkmate::assert_true(!missing(n_subjects_add_list), .var.name = error_n_subjects_add_list)
  
  checkmate::assert_true(!missing(decisions_list),      .var.name = error_decisions_list)
  
  checkmate::assert_class(decisions_list, "decision_list", .var.name = error_decisions_list)
  
  if (is.null(method_name)) {
    
    n_methods <- length(decisions_list$scenario_1$decisions_list)
    
    if (n_methods > 1) {
      
      stop (error_method_name)
      
    } else {
      
      method_name <- names(decisions_list$scenario_1$decisions_list)
      
    }
    
  } else {
    
    method_name <- tryCatch({
      
      match.arg(
        method_name,
        choices    = c('berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'),
        several.ok = FALSE)
      
    }, error = function (e) e)
    
    if (inherits(method_name, "error"))          stop (error_method_name)
    
  }
  
  if (!is.list(n_subjects_add_list)) {
    
    n_subjects_add_list <- rep(list(n_subjects_add_list), length(decisions_list))
    
  }
  
  
  checkmate::assert_list(
    
    n_subjects_add_list, types = c("integer", "numeric"),
    any.missing = FALSE, 
    .var.name = error_n_subjects_add_list
  )
  
  checkmate::assert_list(
    
    n_subjects_add_list, len = length(decisions_list), any.missing = FALSE, 
    .var.name = "The lengths of 'n_subjects_add_list' and 'decisions_list' must be equal"
  )
  
  checkmate::assert_true(
    
    all(vapply(n_subjects_add_list,
               checkmate::test_integerish,
               logical(1),
               lower = 0, any.missing = FALSE)
    ),
    .var.name = error_n_subjects_add_list
  )
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  scenario_numbers <- as.numeric(sub("scenario_", "", names(decisions_list)))
  
  scenario_list <- vector(mode = "list", length = length(decisions_list))
  names(scenario_list) <- paste0("scenario_", scenario_numbers)
  for (s in seq_along(scenario_list)) {
    
    if (!(method_name %in%
          decisions_list[[s]]$analysis_data$analysis_parameters$method_names)) {
      
      stop (simpleError("Selected method_name not analyzed"))
      
    }
    
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
    
    
    checkmate::assert_true(
      length(n_subjects_add) == length(response_rates_new),
      .var.name = "The length of n_subjects_add must be equal ,
      to the length of the response rates that are in (0, 1)"
    )
    
    n_trials <- decisions_list[[s]]$scenario_data$n_trials
    
    add_scenario <- getScenario(
      n_subjects     = n_subjects_add,
      response_rates = response_rates_new,
      cohort_names   = cohort_names_new,
      n_trials       = n_trials)
    
    go_decisions <- decisions_list[[s]]$decisions_list[[method_name]]
    previous_gos <- go_decisions
    
    if ("overall" %in% colnames(go_decisions)) {
      overall_gos  <- go_decisions[, which(colnames(go_decisions) == "overall")]
      go_decisions <- go_decisions[, -which(colnames(go_decisions) == "overall")]
    } else {
      overall_gos <- rep(TRUE, nrow(go_decisions))
    }
    if (!all(index_new %in% as.numeric(sub("decision_", "", colnames(go_decisions))))) {
      stop (simpleError(
        "There must be a decision for each recruiting cohort in the 'decisions_list'"))
    }
    
    go_decisions <- go_decisions[overall_gos, index_new]
    
    n_responders_add <- add_scenario$n_responders[overall_gos, ] * go_decisions
    n_subjects_add   <- add_scenario$n_subjects[overall_gos, ] * go_decisions
    
    n_responders <- decisions_list[[s]]$scenario_data$n_responders
    n_subjects   <- decisions_list[[s]]$scenario_data$n_subjects
    
    n_responders[overall_gos, index_new] <- n_responders[overall_gos, index_new] + n_responders_add
    n_subjects[overall_gos, index_new]   <- n_subjects[overall_gos, index_new] + n_subjects_add
    

    scenario_list[[s]] <- list(
      n_subjects        = n_subjects,
      n_responders      = n_responders,
      response_rates    = response_rates,
      previous_analyses = list(
        go_decisions   = previous_gos,
        post_quantiles = decisions_list[[s]]$analysis_data$quantiles_list),
      n_trials          = n_trials)
    
    scenario_list[[s]]$scenario_number <-
      decisions_list[[s]]$scenario_data$scenario_number
    
  }
  
  class(scenario_list) <- "scenario_list"
  
  return (scenario_list)
  
}

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
  
  error_n_subjects <- 
    "Providing a vector of integers for the argument 'error_n_subjects'"
  error_n_responders <-
    "Providing a vector of integers for the argument 'n_responders'"
  
  checkmate::assert_integerish(n_subjects, any.missing = FALSE, min.len = 1, 
                               .var.name = error_n_subjects)
  checkmate::assert_integerish(n_responders, any.missing = FALSE, min.len = 1,
                               .var.name = error_n_responders)
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  utils::capture.output({
    trial <- simulateScenarios(
      n_subjects_list     = list(as.integer(n_subjects)),
      response_rates_list = list(as.integer(n_responders)),
      n_trials            = 1
    )
  })
  
  return (trial)
  
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

  rownames(n_subjects_matrix) <- format(analysis_dates, format = date_format)
  colnames(n_subjects_matrix) <- paste0("cohort_", seq_len(ncol(n_subjects_matrix)))

  return (n_subjects_matrix)

}

getRecruitment <- function (

  n_subjects_required,
  recruitment_per_month,
  start_date,

  date_format = "%m/%d/%Y"

) {


  checkmate::assert_numeric(
    
    n_subjects_required, lower = 0, any.missing = FALSE,
  .var.name = "Providing a matrix of non-negative integers for the argument 'n_subjects_required'"
  )

  checkmate::assert_true(
    
    checkmate::test_integerish(n_subjects_required, lower = 0),
    .var.name = "Providing a matrix of non-negative integers for the argument 'n_subjects_required'"
    )

  checkmate::assert_numeric(
    
    recruitment_per_month, lower = 0, any.missing = FALSE,
    .var.name = "Providing a vector of non-negative numerics for the argument 'recruitment_per_month'")

  checkmate::assert_string(
    
    start_date, pattern = ".+",
    .var.name = "Please provide a string in the format 'date_format' of for the argument 'start_date'")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  hist_index <- recruitment_per_month == 0

  recruitment_per_day <- recruitment_per_month[!hist_index] * (12 / 365)

  start_date <- as.Date(start_date, format = date_format)

  n_subjects_required <- convertVector2Matrix(n_subjects_required)

  checkmate::assert_true(
    length(recruitment_per_month) == ncol(n_subjects_required))

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

getResponders <- function (
    
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

getScenario <- function(n_subjects, response_rates, 
                        cohort_names = seq_along(n_subjects), n_trials = 1e4
) {
  
  checkmate::assert_integerish(n_subjects, lower = 0, any.missing = FALSE)
  
  checkmate::assert_numeric(response_rates, any.missing = FALSE)
  
  checkmate::assert_integerish(n_trials, lower = 1, any.missing = FALSE, len = 1)
  
  checkmate::assert_true(length(cohort_names) == length(n_subjects),
                         .var.name = "cohort_names must match length of n_subjects")
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  response_rates            <- convertVector2Matrix(response_rates)
  cohort_names_chr          <- as.character(cohort_names)
  colnames(response_rates)  <- paste0("rr_", cohort_names_chr)
  
  new_cohorts  <- FALSE
  hist_cohorts <- FALSE
  
  # New cohorts: 0 < rr < 1
  if (any(response_rates < 1 & response_rates > 0)) {
    new_cohorts <- TRUE
    index_new   <- which(response_rates < 1 & response_rates > 0)
    
    # Match original check (length on sliced matrix/vector)
    checkmate::assert_true(
      
      length(n_subjects[index_new]) == length(response_rates[, index_new]),
      
      .var.name = "n_subjects and response_rates must have same length"
      
    )
    
    n_responders <- getResponders(
      response_rates = response_rates[, index_new, drop = FALSE],
      n_subjects     = n_subjects[index_new],
      n_trials       = n_trials
    )
  }
  
  # Historical cohorts: rr >= 1 OR rr == 0
  if (any(response_rates >= 1 | response_rates == 0)) {
    hist_cohorts <- TRUE
    index_hist   <- which(response_rates >= 1 | response_rates == 0)
    
    checkmate::assert_true(
      
      length(n_subjects[index_hist]) == length(response_rates[, index_hist]),
      
      .var.name = "n_subjects and response_rates must have same length"
      
    )
    
    if (new_cohorts) {
      n_responders_hist <- matrix(
        rep(response_rates[index_hist], each = nrow(n_responders)),
        nrow = nrow(n_responders)
      )
    } else {
      n_responders_hist <- matrix(response_rates[index_hist], nrow = 1)
    }
  }
  
  if (isTRUE(new_cohorts) && isTRUE(hist_cohorts)) {
    
    n_responders <- cbind(n_responders, n_responders_hist)
    
  } else if (isTRUE(hist_cohorts) && !isTRUE(new_cohorts)) {
    
    n_responders <- n_responders_hist
    
  }
  
  previous_gos <- matrix(
    TRUE, byrow = TRUE,
    ncol = length(n_subjects) + 1L,
    nrow = nrow(n_responders)
  )
  colnames(previous_gos) <- c("overall", paste0("decision_", cohort_names_chr))
  
  n_subjects_mat <- matrix(
    n_subjects, byrow = TRUE,
    ncol = length(n_subjects),
    nrow = nrow(n_responders)
  )
  
  colnames(n_subjects_mat) <- paste0("n_", cohort_names_chr)
  colnames(n_responders)   <- paste0("r_", cohort_names_chr)
  
  scenario_data <- list(
    n_subjects        = n_subjects_mat,
    n_responders      = n_responders,
    response_rates    = response_rates,
    previous_analyses = list(
      go_decisions   = previous_gos,
      post_quantiles = NULL
    ),
    n_trials          = n_trials
  )
  
  return(scenario_data)
}

is.scenario_list <- function (x) {
  if (missing(x)) stop ("Please provide an object for the argument 'x'")
  inherits(x, "scenario_list")
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
  
  checkmate::assert_integerish(
    scenario_numbers,
    lower = 1, any.missing = FALSE,
    .var.name = "Providing a vector of positive integers for the argument 'scenario_numbers'"
  )
  
  checkmate::assert_string(load_path, .var.name = "Providing a string for the argument 'load_path'")
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  files <- file.path(load_path, paste0("scenario_data_", scenario_numbers, ".rds"))
  scenario_list <- lapply(files, readRDS)
  
  names(scenario_list) <- paste0("scenario_", scenario_numbers)
  class(scenario_list) <- "scenario_list"
  
  return (scenario_list)
  
}

#' @export
print.scenario_list <- function(x, ...) {
  
  n_scenarios    <- length(x)
  scenario_names <- names(x)
  
  n_cohorts      <- length(x[[1]]$response_rates)
  cohort_names   <- paste0("c_", seq_len(n_cohorts))
  
  response_rates <- lapply(x, function (x) x$response_rates)
  n_subjects     <- getAverageNSubjects(x)
  
  n_trial_realizations  <- x[[1]]$n_trials
  n_unique_realizations <- nrow(getUniqueTrials(x))
  
  cat("scenario_list of ", n_scenarios, " scenario", ifelse(n_scenarios == 1, "", "s"),
      " with ", n_cohorts, " cohort", ifelse(n_cohorts == 1, "", "s"),"\n\n", sep = "")
  for (n in seq_along(scenario_names)) {
    
    df_out <- t(data.frame(c(response_rates[[n]]),
                           n_subjects[[n]]))
    rownames(df_out) <- c("    - true response rates:",
                          "    - average number of subjects:")
    colnames(df_out) <- cohort_names
    
    cat("  -", scenario_names[n], "\n")
    print(df_out)
    
    cat("\n")
    
  }
  
  cat("  -", n_trial_realizations, "trial realizations per scenario\n")
  cat("  -", n_unique_realizations, "unique trial realizations overall\n")
  
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
  
  checkmate::assert(
    checkmate::check_class(scenario_list, "scenario_list"),
    checkmate::check_character(save_path, len = 1),
    combine = "and",
    .var.name = "Please provide an object of class scenario_list for the argument 'scenario_list'"
  )
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  if (!dir.exists(save_path)) {
    dir.create(save_path)
  }
  
  scenario_numbers <- sapply(scenario_list, function (x) x$scenario_number)
  
  for (s in seq_along(scenario_list)) {
    saveRDS(
      scenario_list[[s]],
      file = paste0(save_path, "/scenario_data_", scenario_numbers[s], ".rds")
    )
  }
  
  return (list(scenario_numbers = scenario_numbers, path = save_path))
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
  
  error_n_subjects_list     <- 
    "Providing a list of vectors of positive integers for the argument 'n_subjects_list'"
  error_response_rates_list <- 
    paste("Providing a list of vectors of non-negative numerics for the argument",
          "'response_rates_list'\n", "Values outside of (0, 1) must be integers")
  error_n_trials            <- 
    "Providing a positive integer for the argument 'n_trials'"
  error_scenario_numbers    <- 
    "Providing a vector of positive integers for the argument 'scenario_numbers'"

  checkmate::assert_list(
    response_rates_list,
    types       = "numeric",
    any.missing = FALSE,
    .var.name   = error_response_rates_list
  )
  
  if (!is.list(n_subjects_list)) {
    n_subjects_list <- rep(list(n_subjects_list), length(response_rates_list))
  }
  
  checkmate::assert_list(
    n_subjects_list,
    types       = "numeric",
    any.missing = FALSE,
    .var.name   = error_n_subjects_list
  )
  
  checkmate::assert_true(
    all(vapply(
      n_subjects_list,
      checkmate::test_integerish,
      logical(1),
      lower        = 1,
      any.missing  = FALSE
    )),
    .var.name = error_n_subjects_list
  )
  
  checkmate::assert_integerish(
    scenario_numbers,
    lower       = 1,
    any.missing = FALSE,
    .var.name   = error_scenario_numbers
  )
  
  checkmate::assert_true(
    length(scenario_numbers) == length(n_subjects_list),
    .var.name = "'scenario_numbers' and 'n_subjects_list' must have same length"
  )
  
  checkmate::assert_true(
    length(n_subjects_list) == length(response_rates_list),
    .var.name = "'n_subjects_list' and 'response_rates_list' must have same length"
  )
  
  n_cohorts <- length(response_rates_list[[1L]])
  
  checkmate::assert_true(
    n_cohorts >= 2L,
    .var.name = "Each scenario having at least 2 cohorts"
  )
  
  cohort_lengths_rr <- vapply(response_rates_list, length, integer(1))
  checkmate::assert_true(
    all(cohort_lengths_rr == n_cohorts),
    .var.name = "All scenarios having same number of cohorts in 'response_rates_list'"
  )
  
  cohort_lengths_ns <- vapply(n_subjects_list, length, integer(1))
  checkmate::assert_true(
    all(cohort_lengths_ns == n_cohorts),
    .var.name = "All scenarios having same number of cohorts in 'n_subjects_list'"
  )
  
  for (rates in response_rates_list) {
    
    is_whole <- abs(rates - round(rates)) < .Machine$double.eps^0.5

    ok <- (rates > 0 & rates < 1) |
      (is_whole & (rates == 0 | rates >= 1))
    
    checkmate::assert_true(
      all(ok),
      .var.name = error_response_rates_list
    )
  }
  
  for (i in seq_along(n_subjects_list)) {
    checkmate::assert_true(
      all(n_subjects_list[[i]] >= response_rates_list[[i]]),
      .var.name = "Values in 'response_rates_list' must not exceed 'n_subjects_list'"
    )
  }
  

  if ("n_trials" %in% ls(envir = .GlobalEnv) && missing(n_trials)) {
    n_trials <- get("n_trials", envir = .GlobalEnv)
  }
  
  checkmate::assert_count(
    n_trials,
    positive = TRUE,
    .var.name = error_n_trials
  )
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  scenario_list <- vector(mode = "list", length = length(scenario_numbers))
  
  for (s in seq_along(scenario_numbers)) {
    scenario_list[[s]] <- getScenario(
      n_subjects     = n_subjects_list[[s]],
      response_rates = response_rates_list[[s]],
      n_trials       = n_trials
    )
    scenario_list[[s]]$scenario_number <- scenario_numbers[s]
  }
  
  names(scenario_list) <- paste0("scenario_", scenario_numbers)
  class(scenario_list) <- "scenario_list"
  
  return(scenario_list)
}
