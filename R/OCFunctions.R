
getGoDecisionsByCohort <- function (

  gamma_quantiles,

  boundary_rates = NULL,
  decision_rule  = quote(x > boundary_rates)

  # n_cohorts

) {

  go_decisions_list <- lapply(gamma_quantiles,
                              function (x) {eval(decision_rule)})

  go_decisions <- matrix(unlist(go_decisions_list),
                         nrow = length(go_decisions_list), byrow = TRUE)

  # if (dim(go_decisions)[2] == 1) {
  #   go_decisions <- matrix(rep(go_decisions, n_cohorts),
  #                          nrow = length(go_decisions_list), byrow = FALSE)
  # }

  colnames(go_decisions) <- paste0("decision_", seq_len(ncol(go_decisions)))

  return (go_decisions)

}

getNumericGammaIndex <- function (

  g_numeric,
  quantiles

) {

  if (is.numeric.in.zero.one(g_numeric)) {

    gamma_index <- which(round(1 - quantiles, 5) == round(g_numeric, 5))

    if (!(is.numeric(gamma_index) && length(gamma_index) > 0)) {

      stop (simpleError(paste0(
        "gamma must be one of ",
        paste(round(1 - quantiles, 5), collapse = ", "))))

    }

  } else {

    stop (simpleError(
      "gamma_levels must consist of posterior quantiles, 'mean' or 'sd'"))

  }

  return (gamma_index)

}

getGammaIndices <- function (

  gamma_levels,
  quantiles

) {

  gamma_indices <- sapply(gamma_levels, function (g) {

    if (is.character(g)) {

      if (g == "mean") {

        gamma_index <- length(quantiles) + 1L

      } else if (g == "sd") {

        gamma_index <- length(quantiles) + 2L

      } else {

        g_numeric <- tryCatch({as.numeric(g)}, warning = function(w) w)

        if (inherits(g_numeric, "warning")) stop(simpleError(paste(
          "The only strings allowed for the argument 'evidence_levels' are",
          "'mean' and 'sd'")))

        gamma_index <- getNumericGammaIndex(g_numeric, quantiles)

      }

    } else {

      gamma_index <- getNumericGammaIndex(g, quantiles)

    }

    # if (is.numeric(g)) {
    #
    #   gamma_index <- which(round(1 - quantiles, 5) == round(g, 5))
    #
    #   if (!(is.numeric(gamma_index) && length(gamma_index) > 0)) {
    #
    #     stop (simpleError(paste0(
    #       'gamma must be one of ',
    #       paste(round(1 - quantiles, 5), collapse = ", "))))
    #
    #   }
    #
    # } else if (g == "mean") {
    #
    #   gamma_index <- length(quantiles) + 1L
    #
    # } else if (g == "sd") {
    #
    #   gamma_index <- length(quantiles) + 2L
    #
    # } else {
    #
    #   stop ("gamma_levels must consist of posterior quantiles,
    #         'mean' or 'sd'")
    #
    # }

    return (gamma_index)

  })

  return (gamma_indices)

}

getPosteriorGammaQuantiles <- function (

  method_name,
  posterior_quantiles_list,
  gamma_levels,
  quantiles,
  cohort_names

) {

  gamma_indices <- getGammaIndices(gamma_levels = gamma_levels,
                                   quantiles    = quantiles)

  posterior_quantiles <-
    posterior_quantiles_list[[gsub("_mu", "", method_name)]]

  cohort_indices <- sapply(cohort_names, function (n) {
    grep(n, colnames(posterior_quantiles[[1]]), fixed = TRUE)[1]
  })

  if (length(gamma_indices) != length(cohort_indices)) {

    stop (paste("The number of gamma indices must be equal",
                 "to the number of cohorts to be analyzed."))

  }

  if (grepl("berry", method_name) | grepl("exnex", method_name)) {

    posterior_gamma_quantiles <- lapply(posterior_quantiles, function (x) {
      x[gamma_indices, cohort_indices]
    })

  } else if (method_name == "stratified" | method_name == "pooled") {

    posterior_gamma_quantiles <- lapply(posterior_quantiles, function (x) {
      x[gamma_indices, cohort_indices]
    })

  }

  if (length(gamma_indices) > 1) {
    posterior_gamma_quantiles <- lapply(posterior_gamma_quantiles, diag)
  }

  for (i in seq_along(posterior_gamma_quantiles)) {

    names(posterior_gamma_quantiles[[i]]) <-
      colnames(posterior_quantiles[[1]])[cohort_indices]

  }

  return (posterior_gamma_quantiles)

}

#' @title getGoProbabilities
#' @md
#' @description Calculates the Go probabilities for given decisions
#' @param go_decisions_list An object of class `decision_list`,
#' as returned by \code{\link[bhmbasket]{getGoDecisions}}
#' @param nogo_decisions_list An object of class `decision_list`,
#' as returned by \code{\link[bhmbasket]{getGoDecisions}}, Default: `NULL`
#' @return A list of matrices of Go (and Consider and NoGo) probabilities
#' @details If only `go_decisions_list` is provided
#' (i.e. `nogo_decisions_list` is `NULL`),
#' only Go probabilities will be calculated.
#' If both `go_decisions_list` and `nogo_decisions_list` are provided,
#' Go, Consider, and NoGo probabilities will be calculated.
#' @seealso
#'  \code{\link[bhmbasket]{getGoDecisions}}
#' @rdname getGoProbabilities
#' @examples
#' scenarios_list <- simulateScenarios(
#'   n_subjects_list     = list(c(10, 20)),
#'   response_rates_list = list(rep(0.9, 2)),
#'   n_trials            = 10)
#'
#' analyses_list <- performAnalyses(
#'   scenario_list       = scenarios_list,
#'   target_rates        = rep(0.5, 2),
#'   n_mcmc_iterations   = 100,
#'   n_cores             = 1L)
#'
#' go_decisions_list <- getGoDecisions(
#'   analyses_list       = analyses_list,
#'   cohort_names        = c("p_1", "p_2"),
#'   evidence_levels     = c(0.5, 0.8),
#'   boundary_rules      = quote(c(x[1] > 0.8, x[2] > 0.6)))
#'
#' nogo_decisions_list <- getGoDecisions(
#'   analyses_list       = analyses_list,
#'   cohort_names        = c("p_1", "p_2"),
#'   evidence_levels     = c(0.5, 0.8),
#'   boundary_rules      = quote(c(x[1] < 0.5, x[2] < 0.3)))
#'
#' getGoProbabilities(go_decisions_list)
#' getGoProbabilities(go_decisions_list, nogo_decisions_list)
#' @author Stephan Wojciekowski
#' @export
getGoProbabilities <- function (

  go_decisions_list,
  nogo_decisions_list = NULL

) {

  error_go_decisions_list <-
    simpleError( paste("Please provide an object of class decision_list",
                       "for the argument 'go_decisions_list'"))
  error_nogo_decisions_list <- simpleError(paste(
    "Please provide either NULL or an object of class decision_list",
    "for the argument 'nogo_decisions_list'"))

  if (missing(go_decisions_list))             stop (error_go_decisions_list)

  if (!is.decision_list(go_decisions_list))   stop (error_go_decisions_list)
  if (!is.null(nogo_decisions_list) &&
      !is.decision_list(nogo_decisions_list)) stop (error_nogo_decisions_list)

  if (!is.null(nogo_decisions_list) &&
      !identical(dim(go_decisions_list[[1]]$decisions_list[[1]]),
                 dim(nogo_decisions_list[[1]]$decisions_list[[1]]))) stop(
                   simpleError(paste("The decision_lists 'go_decisions_list'",
                                     "and 'go_decisions_list'",
                                     "do not follow the same format")))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

  go_probs_per_scenario <- vector(mode = "list",
                                  length = length(go_decisions_list))
  names(go_probs_per_scenario) <- names(go_decisions_list)

  for (s in seq_along(go_probs_per_scenario)) {

    method_names <- names(go_decisions_list[[s]]$decisions_list)

    go_probs_per_method_list <- vector(mode = "list",
                                       length = length(method_names))
    names(go_probs_per_method_list) <- method_names

    for (n in seq_along(method_names)) {

      go_decisions     <- go_decisions_list[[s]]$decisions_list[[n]]
      decisions_matrix <- t(as.matrix(colMeans(go_decisions)))
      row.names(decisions_matrix) <- "Go"

      if (!is.null(nogo_decisions_list)) {

        nogo_decisions <- nogo_decisions_list[[s]]$decisions_list[[n]]
        nogo_probs     <- t(as.matrix(colMeans(nogo_decisions)))
        consider_probs <- round(1 - nogo_probs - decisions_matrix, 9)

        # calculate consider decision probabilities directly:
        # consider_decisions <- sapply(2:5, function (i) {
        #   !apply(cbind(nogo_decisions[, i], go_decisions[, i]), 1, any)
        # })
        # consider_decisions <- cbind(apply(consider_decisions, 1, any) &
        #                               !go_decisions[, 1],
        #                             consider_decisions)
        # consider_probs     <- t(as.matrix(colMeans(consider_decisions)))

        if (!isTRUE(all.equal(sum(go_decisions * nogo_decisions), 0)))
          stop (paste("There are cohorts for which both go and nogo decisions",
                      "are TRUE. Please revise your decision rules."))

        decisions_matrix <- rbind(decisions_matrix, consider_probs, nogo_probs)
        row.names(decisions_matrix) <- c("Go", "Consider", "NoGo")

      }

      go_probs_per_method_list[[n]] <- decisions_matrix

    }

    go_probs_per_scenario[[s]] <- go_probs_per_method_list

  }

  go_probs_per_method <- listPerMethod(go_probs_per_scenario)

  return (go_probs_per_method)

}

#' @title negateGoDecisions
#' @md
#' @description Negates the go decisions derived with
#' \code{\link[bhmbasket]{getGoDecisions}}.
#' @param go_decisions_list An object of class `decision_list`,
#' as returned by \code{\link[bhmbasket]{getGoDecisions}}
#' @return A list of NoGo decisions of class `decision_list`
#' @details This function is intended for implementing decision rules with a
#' consider zone as
#' e.g. proposed in "Bayesian design of proof-of-concept trials" by
#' Fisch et al. (2015).
#' This approach involves two criteria, Significance and Relevance.
#' \itemize{
#'   \item Significance: high evidence that the treatment effect is greater
#'   than some smaller value (e.g. treatment effect under H0)
#'   \item Relevance: moderate evidence that the treatment effect is greater
#'   than some larger value (e.g. treatment effect under a certain alternative)
#' }
#' The decision for a cohort is then taken as follows:
#' \itemize{
#'   \item Go decision: Significance and Relevance
#'   \item Consider decision: either Significance, or Relevance, but not both
#'   \item NoGo decision: no Significance and no Relevance
#' }
#' In the example below, the following criteria for are implemented for each of
#' the three cohorts:
#' \itemize{
#'   \item Significance: \eqn{P(p_j > 0.4) > 0.95}
#'   \item Relevance: \eqn{P(p_j > 0.8) > 0.5}
#' }
#' @seealso
#'  \code{\link[bhmbasket]{getGoDecisions}}
#' @rdname negateGoDecisions
#' @examples
#' scenarios_list <- simulateScenarios(
#'   n_subjects_list     = list(c(10, 20, 30)),
#'   response_rates_list = list(rep(0.9, 3)),
#'   n_trials            = 10)
#'
#' analysis_list <- performAnalyses(
#'   scenario_list      = scenarios_list,
#'   target_rates       = rep(0.5, 3),
#'   n_mcmc_iterations  = 100,
#'   n_cores            = 1L)
#'
#' go_decisions_list <- getGoDecisions(
#'   analyses_list   = analysis_list,
#'   cohort_names    = c("p_1", "p_2", "p_3",
#'                       "p_1", "p_2", "p_3"),
#'   evidence_levels = c(0.5,  0.5,  0.5,
#'                       0.95, 0.95, 0.95),
#'   boundary_rules  = quote(c(x[1] > 0.8 & x[4] > 0.4,
#'                             x[2] > 0.8 & x[5] > 0.4,
#'                             x[3] > 0.8 & x[6] > 0.4)))
#'
#' nogo_decisions <- negateGoDecisions(getGoDecisions(
#'   analyses_list   = analysis_list,
#'   cohort_names    = c("p_1", "p_2", "p_3",
#'                       "p_1", "p_2", "p_3"),
#'   evidence_levels = c(0.5,  0.5,  0.5,
#'                       0.95, 0.95, 0.95),
#'   boundary_rules  = quote(c(x[1] > 0.8 | x[4] > 0.4,
#'                             x[2] > 0.8 | x[5] > 0.4,
#'                             x[3] > 0.8 | x[6] > 0.4))))
#'
#' getGoProbabilities(go_decisions_list, nogo_decisions)
#' @author Stephan Wojciekowski
#' @references Fisch, Roland, et al.
#' "Bayesian design of proof-of-concept trials."
#' \emph{Therapeutic innovation & regulatory science} 49.1 (2015): 155-162.
#' @export
negateGoDecisions <- function (

  go_decisions_list

) {

  error_go_decisions_list <-
    simpleError(paste("Please provide an object of class decision_list",
                      "for the argument 'go_decisions_list'"))

  if (missing(go_decisions_list))           stop (error_go_decisions_list)

  if (!is.decision_list(go_decisions_list)) stop (error_go_decisions_list)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

  for (s in seq_along(go_decisions_list)) {

    for (m in seq_along(go_decisions_list[[s]]$decisions_list)) {

      go_decisions_list[[s]]$decisions_list[[m]] <-
        !go_decisions_list[[s]]$decisions_list[[m]]

      if (ncol(go_decisions_list[[s]]$decisions_list[[m]]) > 1) {

        go_decisions_list[[s]]$decisions_list[[m]][, 1] <-
          apply(go_decisions_list[[s]]$decisions_list[[m]][, -1], 1, all)

      }

    }

  }

  return (go_decisions_list)

}

getScenarioNumbers <- function (analyses_list) {

  as.numeric(sub("scenario_", "", names(analyses_list)))

}

getGoBoundaries <- function (

    scenario_analysis_list,
    cohort_names,
    go_rates,
    gamma_levels_list,
    method_names

) {"dummy function"}
# getGoBoundaries <- function (
#
#   scenario_analysis_list,
#   cohort_names,
#   go_rates,
#   gamma_levels_list,
#   method_names
#
# ) {
#
#   # getBoundaries <- function (
#   #
#   #   scenario_analysis_list,
#   #   cohort_names,
#   #   go_rates,
#   #   gamma_levels_list,
#   #   method_name
#   #
#   # ) {
#   #
#   #   boundaries <- sapply(seq_along(cohort_names), function (n) {
#   #
#   #     stats::uniroot(
#   #
#   #       f = function (x) {
#   #
#   #     go_probs_list <- getGoProbabilities(
#   #       getGoDecisions(
#   #         analyses_list = scenario_analysis_list,
#   #         cohort_names           = cohort_names[n],
#   #         boundary_rules_list    = bquote(x[1] > .(x)),
#   #         gamma_levels_list      = bquote(list(.(gamma_levels_list[[n]])))))
#   #
#   #     return (go_probs_list[[method_name]][[1]][1, 1] - go_rates[n])
#   #
#   #       },
#   #
#   #       interval = c(0, 1)
#   #
#   #     )$root
#   #
#   #   })
#   #
#   #   return (boundaries)
#   #
#   # }
#
#   # boundaries <- getBoundaries(
#   #   scenario_analysis_list = scenario_analysis_list,
#   #   cohort_names           = cohort_names,
#   #   go_rates               = go_rates,
#   #   gamma_levels_list      = gamma_levels_list,
#   #   method_name            = method_name)
#   #
#   # boundary_rules <- str2expression(paste0("c(", paste0(
#   #   paste0("x[", seq_along(boundaries), "] > ", boundaries),
#   #          collapse = ", "), ")"))
#   #
#   # decisions <- getGoDecisions(
#   #   analyses_list = scenario_analysis_list,
#   #   cohort_names           = cohort_names,
#   #   boundary_rules_list    = boundary_rules,
#   #   gamma_levels_list      = gamma_levels_list)
#   #
#   # getGoProbabilities(decisions)[[method_name]][[1]][1, 1]
#
#   boundary_list <- vector(mode = "list", length = length(method_names))
#   names(boundary_list) <- method_names
#
#   for (k in seq_along(method_names)) {
#
#     boundary_list[[k]] <- sapply(seq_along(cohort_names), function (n) {
#
#       stats::uniroot(
#
#         f = function (x) {
#
#         go_probs_list <- getGoProbabilities(
#           getGoDecisions(
#             analyses_list = scenario_analysis_list,
#             cohort_names           = cohort_names[n],
#             boundary_rules_list    = bquote(x[1] > .(x)),
#             gamma_levels_list      = bquote(list(.(gamma_levels_list[[n]])))))
#
#         return (go_probs_list[[method_names[k]]][[1]][1, 1] - go_rates[n])
#
#         },
#
#         interval = c(0, 1)
#
#       )$root
#
#     })
#
#   }
#
#   if (length(boundary_list) == 1) {
#
#     boundary_list <- boundary_list[[1]]
#
#   }
#
#   return (boundary_list)
#
# }

#' @title getGoDecisions
#' @description This function applies decision rules to the analyzed trials.
#' The resulting \code{decision_list} can be further processed with
#' \code{\link[bhmbasket]{getGoProbabilities}} or
#' \code{\link[bhmbasket]{continueRecruitment}}.
#' @param analyses_list An object of class \code{analysis_list},
#' as created with \code{\link[bhmbasket]{performAnalyses}}
#' @param cohort_names A vector of strings with the names of the cohorts, e.g.
#' \code{c('p_1', 'p_2')}
#' @param evidence_levels A vector of numerics in \code{(0, 1)} for the
#' posterior probability thresholds for the cohorts.
#' Will be recycled to match the number of methods in the \code{analyses_list}
#' @param boundary_rules A quote of a vector for the boundary rules,
#' \code{quote(c(...))}, see details.
#' The number of decisions to be taken must match the number of cohorts.
#' Will be recycled to match the number of methods in the \code{analyses_list}
#' @param overall_min_gos A positive integer for the minimum number of
#' cohort-wise go decisions required for an overall go decision
#'  Default: \code{1}
#' @return An object of class \code{decision_list}
#' @details This function applies decision rules of the following type to the
#' outcomes of (simulated) basket trials with binary endpoints:
#' \deqn{P(p_j|data > p_{B,j}) > \gamma,}
#' where \eqn{p_j|data} is the posterior response rate of cohort \eqn{j},
#' \eqn{p_{B,j}} is the response rate boundary of cohort \eqn{j},
#' and \eqn{\gamma} is the evidence level.
#' This rule can equivalently be written as \deqn{q_{1-\gamma,j} > p_{B,j},}
#' where \eqn{q_{1-\gamma,j}} is the \eqn{1-\gamma}-quantile of the posterior
#' response rate of cohort \eqn{j}.
#'
#' The arguments \code{cohort_names} and \code{evidence_levels} determine
#' \eqn{q_{1-\gamma,j}}, where the entries of \code{cohort_names} and
#' \code{evidence_levels} are matched corresponding to their order.
#'
#' The argument \code{boundary_rules} provides the rules that describe what
#' should happen with  the posterior quantiles \eqn{q_{1-\gamma,j}}.
#' The first posterior quantile determined by the first items of
#' \code{cohort_names} and \code{evidence_levels} is referred to as \code{x[1]},
#' the second as \code{x[2]}, etc.
#' Using the \code{quote(c(...))}-notation,
#' many different rules can be implemented.
#' A decision rule for only one cohort would be
#' \code{boundary_rules = quote(c(x[1] > 0.1))},
#' \code{cohort_names = 'p_1'}, and \code{evidence_levels = 0.5},
#' which implements the rule \eqn{P(p_1|data > 0.1) > 0.5}.
#' The number of decisions to be taken must match the number of cohorts, i.e.
#' for each cohort there must be a decision rule in the vector separated by a comma.
#' See the example section for a decision rule for more than one cohort and
#' the example of \code{\link[bhmbasket]{negateGoDecisions}}
#' for the implementation of a more complex decision rule.
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getGoProbabilities}}
#'  \code{\link[bhmbasket]{negateGoDecisions}}
#'  \code{\link[bhmbasket]{continueRecruitment}}
#' @rdname getGoDecisions
#' @examples
#' scenarios_list <- simulateScenarios(
#'   n_subjects_list     = list(c(10, 20, 30)),
#'   response_rates_list = list(c(0.1, 0.1, 0.9)),
#'   n_trials            = 10)
#'
#' analyses_list <- performAnalyses(
#'   scenario_list      = scenarios_list,
#'   target_rates       = rep(0.5, 3),
#'   n_mcmc_iterations  = 100,
#'   n_cores            = 1L)
#'
#' ## Decision rule for more than one cohort
#' decisions_list <- getGoDecisions(
#'   analyses_list   = analyses_list,
#'   cohort_names    = c("p_1", "p_2", "p_3"),
#'   evidence_levels = c(0.5, 0.5, 0.8),
#'   boundary_rules  = quote(c(x[1] > 0.7, x[2] < 0.3, x[3] < 0.6)))
#'
#' ## Decision rule for only two of the three cohorts
#' decisions_list <- getGoDecisions(
#'   analyses_list   = analyses_list,
#'   cohort_names    = c("p_1", "p_3"),
#'   evidence_levels = c(0.5, 0.8),
#'   boundary_rules  = quote(c(x[1] > 0.7, TRUE, x[3] < 0.6)),
#'   overall_min_gos = 2L)
#'
#' ## Different decision rules for each method
#' ## This works the same way for the different evidence_levels
#' decisions_list <- getGoDecisions(
#'   analyses_list   = analyses_list,
#'   cohort_names    = c("p_1", "p_2", "p_3"),
#'   evidence_levels = c(0.5, 0.5, 0.8),
#'   boundary_rules  = list(quote(c(x[1] > 0.1, x[2] < 0.5, x[3] < 0.1)),  # "berry"
#'                          quote(c(x[1] > 0.2, x[2] < 0.4, x[3] < 0.2)),  # "exnex"
#'                          quote(c(x[1] > 0.3, x[2] < 0.3, x[3] < 0.3)),  # "exnex_adj"
#'                          quote(c(x[1] > 0.4, x[2] < 0.2, x[3] < 0.4)),  # "pooled"
#'                          quote(c(x[1] > 0.5, x[2] < 0.1, x[3] < 0.5)))) # "stratified"
#' @author Stephan Wojciekowski
#' @export
getGoDecisions <- function (

  analyses_list,

  cohort_names,
  evidence_levels,
  boundary_rules,

  overall_min_gos = 1

) {

  error_analyses_list <- simpleError(
    "Please provide an object of class analysis_list for the argument 'analyses_list'")
  error_cohort_names <- simpleError(paste(
    "Please provide a vector of strings for the argument 'cohort_names',",
    "e.g. c('p_1','p_2')"))
  error_evidence_levels <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'evidence_levels'")
  error_boundary_rules <- simpleError(paste(
    "Please provide a quote(c(...)) for the argument 'boundary_rules.'",
    "The vector c(...) inside the quote() must have the same length as the number of cohorts.",
    "See ?getGoDecisions for details"))
  error_overall_min_gos <- simpleError(
    "Please privide a positive integer for the argument 'overall_min_gos'")

  if (missing(analyses_list))           stop (error_analyses_list)
  if (missing(cohort_names))            stop (error_cohort_names)
  if (missing(evidence_levels))         stop (error_evidence_levels)
  if (missing(boundary_rules))          stop (error_boundary_rules)

  if (!is.analysis_list(analyses_list)) stop (error_analyses_list)
  if (!is.character(cohort_names))      stop (error_cohort_names)

  if (any(!cohort_names %in% colnames(analyses_list[[1]]$quantiles_list[[1]][[1]]))) stop (
    simpleError("The specified cohorts do not match the cohorts analyzed in 'analyses_list'"))

  checkEvidenceLevels <- function (evidence_levels) {

    if (!identical(length(evidence_levels), length(cohort_names))) stop(simpleError(
      "The 'evidence_levels' and the 'cohort_names' must have the same length"))

    if (is.character(evidence_levels)) {

      mean_index <- evidence_levels == "mean"

      evidence_levels_numeric <- tryCatch({
        as.numeric(evidence_levels[!mean_index])
      }, warning = function(w) w)

      if (inherits(evidence_levels_numeric, "warning")) stop(simpleError(
        "The only string allowed for the argument 'evidence_levels' is 'mean'"))


    } else {

      evidence_levels_numeric <- evidence_levels

    }

    if (!is.numeric.in.zero.one(evidence_levels_numeric)) stop (error_evidence_levels)

    available_quantiles <- round(analyses_list[[1]]$analysis_parameters$quantiles, 9)
    asked_quantiles     <- round(1 - evidence_levels_numeric, 9)
    if (any(!asked_quantiles %in% available_quantiles)) stop (simpleError(paste(
      "The 'evidence_levels' must have matches",
      "in the 'evidence_levels' provided to the call performAnalyses()",
      "that created the 'analyses_list'")))

  }
  if (is.list(evidence_levels)) {
    for (i in seq_along(evidence_levels)) {
      checkEvidenceLevels(evidence_levels[[i]])
    }
  } else {
    checkEvidenceLevels(evidence_levels)
  }

  check_boundary_rules <- tryCatch({
    x <- stats::runif(n = length(cohort_names), min = 0.001, max = 0.999)
    if (is.list(boundary_rules)) {
      for (i in seq_along(boundary_rules)) {
        if (!is.language(boundary_rules[[i]]))                                   stop ()
        if (!identical(boundary_rules[[i]][1], quote(c())))                      stop ()
        if (!identical(length(boundary_rules[[i]]) - 1L,
                       ncol(analyses_list$scenario_1$scenario_data$n_subjects))) stop ()
        eval(boundary_rules[[i]])
      }
    } else {
      if (!is.language(boundary_rules))                                        stop ()
      if (!identical(boundary_rules[1], quote(c())))                           stop ()
      ## fix number of decisions to number of cohorts
      if (!identical(length(boundary_rules) - 1L,
                     ncol(analyses_list$scenario_1$scenario_data$response_rates))) stop ()
      eval(boundary_rules)
    }
    rm(x)
  }, error = function (e) e)
  if (inherits(check_boundary_rules, "error"))     stop(error_boundary_rules)

  if (!is.positive.wholenumber(overall_min_gos))   stop (error_overall_min_gos)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  gamma_levels <- evidence_levels

  # gamma_levels <- str2expression(deparse(gamma_levels))

  ## Get method names
  method_names_matrix <- t(sapply(analyses_list,
                                  function (x) x$analysis_parameters$method_names))
  if (!all(sapply(seq_len(nrow(method_names_matrix)),
                  function (x) identical(method_names_matrix[1, ], method_names_matrix[x, ])))) {
    stop ("The scenarios where analysed with different methods")
  }
  method_names <- method_names_matrix[1, ]

  ## in case only one method has been used for all scenarios
  if (all(sapply(seq_along(method_names), function (x) {
    method_names[1] == method_names[x]
  }))) {
    method_names <- method_names[1]
  }

  ## check for input consistency
  if (!is.list(boundary_rules)) {
    boundary_rules <- list(boundary_rules)
  }
  if (!is.list(gamma_levels)) {
    gamma_levels <- list(gamma_levels)
  }

  if (length(method_names) < length(boundary_rules)) {
    stop (paste0("The lengths of 'boundary_rules' must be less than or equal to",
                 " the length of 'method_names'"))
  } else if (length(method_names) > length(boundary_rules)) {
    boundary_rules <- rep(boundary_rules, length.out = length(method_names))
  }
  if (length(method_names) < length(gamma_levels)) {
    stop (paste0("The lengths of 'evidence_levels' must be less than or equal to",
                 " the length of 'method_names'"))
  } else if (length(method_names) > length(gamma_levels)) {
    gamma_levels <- rep(gamma_levels, length.out = length(method_names))
  }

  decisions_list <- vector(mode = "list", length = length(analyses_list))
  names(decisions_list) <- names(analyses_list)

  for (s in seq_along(decisions_list)) {

    analysis_data <- analyses_list[[s]]

    methods_decisions_list <- vector(mode = "list", length = length(method_names))
    names(methods_decisions_list) <- method_names

    for (n in seq_along(method_names)) {

      go_decisions <- getGoDecisionsByCohort(
        gamma_quantiles = getPosteriorGammaQuantiles(
          method_name              = method_names[n],
          gamma_levels             = gamma_levels[[n]], # eval() ??
          quantiles                = analysis_data$analysis_parameters$quantiles,
          posterior_quantiles_list = analysis_data$quantiles_list,
          cohort_names             = cohort_names),
        decision_rule   = boundary_rules[[n]])

      ## combine new decision outcomes with previous decisions (and convert to logical)
      previous_gos <- analysis_data$scenario_data$previous_analyses$go_decisions[, -1]
      go_decisions <- go_decisions * previous_gos > 0

      ## Overall go:
      overall_go   <- apply(go_decisions, 1, function (x) sum (x) >= overall_min_gos)
      go_decisions <- cbind(overall = overall_go, go_decisions)

      ## store
      methods_decisions_list[[n]] <- go_decisions

    }

    decisions_list[[s]] <-
      list(decisions_list = methods_decisions_list,
           analysis_data  = list(
             quantiles_list      = analyses_list[[s]]$quantiles_list,
             analysis_parameters = analyses_list[[s]]$analysis_parameters),
           scenario_data  = analyses_list[[s]]$scenario_data)
  }

  names(decisions_list) <- names(analyses_list)

  class(decisions_list) <- "decision_list"

  return (decisions_list)

}

is.decision_list <- function (x) {

  if (missing(x)) stop ("Please provide an object for the argument 'x'")

  inherits(x, "decision_list")

}

#' @title getEstimates
#' @md
#' @description This function calculates the point estimates and credible intervals per cohort,
#' as well as estimates of the biases and the mean squared errors of the point estimates per cohort
#' @param analyses_list An object of class `analysis_list`,
#' as created with \code{\link[bhmbasket]{performAnalyses}}
#' @param add_parameters A vector of strings naming additional parameters
#' from the Bayesian hierarchical models, e.g. `c('mu', 'tau')`.
#' If `NULL`, no additional parameters will be evaluated,
#' Default: `NULL`
#' @param point_estimator A string indicating the type of estimator used for calculation of
#' bias and MSE. Must be one of `'median'` or `'mean'`
#' @param alpha_level A numeric in (0, 1) for the level of the credible interval.
#' Only values corresponding to quantiles saved in \code{\link[bhmbasket]{performAnalyses}}
#' will work, Default: `0.05`
#' @details Bias and MSE will only be calculated for response rate estimates of simulated trials.
#' For additional parameters, bias and MSE will not be calculated.
#'
#' Possible additional parameters are for the Bayesian hierarchical models are
#' `c('mu', 'tau')` for `'berry'`, `'exnex'`, and `'exnex_adj'`.
#' The latter two models can also access the posterior weights
#' `paste0("w_", seq_len(n_cohorts))`.
#' @return A named list of matrices of estimates of response rates and credible intervals.
#' Estimates of bias and MSE are included for response rate estimates of simulated trials.
#' @rdname getEstimates
#' @seealso
#'  \code{\link[bhmbasket]{createTrial}}
#'  \code{\link[bhmbasket]{performAnalyses}}
#' @examples
#'   scenarios_list <- simulateScenarios(
#'     n_subjects_list     = list(c(10, 20, 30)),
#'     response_rates_list = list(c(0.1, 0.2, 3)),
#'     n_trials            = 10)
#'
#'   analyses_list <- performAnalyses(
#'     scenario_list       = scenarios_list,
#'     target_rates        = c(0.1, 0.1, 0.1),
#'     calc_differences    = matrix(c(3, 2, 2, 1), ncol = 2),
#'     n_mcmc_iterations   = 100,
#'     n_cores             = 1L)
#'
#'   getEstimates(analyses_list)
#'   getEstimates(analyses_list   = analyses_list,
#'                add_parameters  = c("mu", "tau", "w_1", "w_2", "w_3"),
#'                point_estimator = "mean",
#'                alpha_level     = 0.1)
#'
#'   outcome <- createTrial(
#'     n_subjects          = c(10, 20, 30),
#'     n_responders        = c( 1,  2,  3))
#'
#'   outcome_analysis <- performAnalyses(
#'     scenario_list       = outcome,
#'     target_rates        = c(0.1, 0.1, 0.1),
#'     n_mcmc_iterations   = 100,
#'     n_cores             = 1L)
#'
#'   getEstimates(outcome_analysis)
#'   getEstimates(analyses_list  = outcome_analysis,
#'                add_parameters = c("mu", "w_1", "w_2", "w_3"))
#' @author Stephan Wojciekowski
#' @export
getEstimates <- function (

  analyses_list,
  add_parameters  = NULL,
  point_estimator = "median",
  alpha_level     = 0.05

) {

  error_analyses_list <- simpleError(
    "Please provide an object of class analysis_list for the argument 'analyses_list'")
  error_add_parameters  <- simpleError(paste(
    "Please provide a either NULL or vector of strings for the argument 'add_parameters'"))
  error_point_estimator <- simpleError(
    "Please provide either 'median' or 'mean' for the argument 'point_estimator'")
  error_alpha_level <- simpleError(
    "Please provide a numeric in (0, 1) for the argument 'alpha_level'")

  if (missing(analyses_list))                      stop (error_analyses_list)

  if (!is.analysis_list(analyses_list))            stop (error_analyses_list)

  point_estimator <- tryCatch({

    match.arg(
      point_estimator,
      choices    = c("median", "mean"),
      several.ok = FALSE)

  }, error = function (e) e)

  if (!is.null(add_parameters) &&
      !is.character(add_parameters))               stop (error_add_parameters)
  if (inherits(point_estimator, "error"))          stop (error_point_estimator)
  if (!is.single.numeric.in.zero.one(alpha_level)) stop (error_alpha_level)

  available_quantiles <- round(analyses_list[[1]]$analysis_parameters$quantiles, 9)
  asked_quantiles     <- round(c(alpha_level / 2, 1 - alpha_level / 2), 9)
  if (any(!asked_quantiles %in% available_quantiles)) stop (simpleError(paste(
    "The 'alpha_level' must be among the stored quantiles in 'analyses_list',",
    "e.g. 1 - alpha_level must be among the evidence_levels in performAnalyses()")))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  cohort_names <- getAllCohortNames(analyses_list)
  diff_indices <- grepl("diff", cohort_names)

  ## Lists to hold the results
  results_list <- vector(mode = "list", length = length(analyses_list))
  names(results_list) <- names(analyses_list)

  for (s in seq_along(analyses_list)) {

    true_rr <- analyses_list[[s]]$scenario_data$response_rates

    ## calculate the historic response rates
    hist_index <- true_rr <= 0 | true_rr >= 1
    if (any(hist_index)) {

      hist_rr <- sapply(which(hist_index), function (x) {

        true_rr[x] / analyses_list[[s]]$scenario_data$n_subjects[1, x]

      })

      true_rr[hist_index] <- hist_rr

    }

    ## calculate the differences in true rr if necessary
    if (any(diff_indices)) {

      diff_cohorts <- sapply(strsplit(substrRight(cohort_names[diff_indices], 2), ""), as.numeric)

      true_diff_rr <- t(apply(diff_cohorts, 2, function (x) {

        diff(true_rr[x])

      }))

    } else {

      true_diff_rr <- NULL

    }

    ## must be after historic response rates
    true_rr <- cbind(true_rr, true_diff_rr)

    ## prepare loop
    method_names <- analyses_list[[s]]$analysis_parameters$method_names

    results_per_method_list <- vector(mode = "list", length = length(method_names))
    names(results_per_method_list) <- method_names

    ## prepare check that additional parameters occur in at least one of the methods
    if (!is.null(add_parameters)) {
      occurences <- vector(length = length(method_names))
    }

    for (n in seq_along(method_names)) {

      if (!is.null(add_parameters)) {

        add_index <- add_parameters %in%
          colnames(analyses_list[[s]]$quantiles_list[[n]][[1]])

        if (any(add_index)) {
          occurences[n] <- TRUE
        }

        parameter_names <- c(cohort_names, add_parameters[add_index])

      } else {

        parameter_names <- cohort_names

      }

      ## Setting up gamma levels
      n_parameters <- length(parameter_names)
      gamma_levels_list <- list(rep("mean", n_parameters),
                                rep("sd", n_parameters),
                                rep(1 - alpha_level / 2, n_parameters),
                                rep(0.5, n_parameters),
                                rep(alpha_level / 2, n_parameters))

      matrix_estimates_list <- lapply(gamma_levels_list, function (gamma_levels) {

        do.call(rbind, getPosteriorGammaQuantiles(
          method_name              = method_names[n],
          gamma_levels             = gamma_levels,
          quantiles                = analyses_list[[s]]$analysis_parameters$quantiles,
          posterior_quantiles_list = analyses_list[[s]]$quantiles_list,
          cohort_names             = parameter_names))

      })

      ## Mean, SD & Quantiles
      post_quantiles <- t(do.call(rbind, lapply(matrix_estimates_list, colMeans)))
      # if (is.null(dim(post_quantiles))) {
      #   post_quantiles <- matrix(post_quantiles, byrow = TRUE,
      #                            ncol = 5, nrow = length(true_rr))
      #   rownames(post_quantiles) <- paste0("p_", seq_along(true_rr))
      # }
      colnames(post_quantiles) <- c(
        "Mean", "SD",
        paste0(c(alpha_level / 2, 0.5, 1 - alpha_level / 2) * 100, "%"))

      ## Bias and MSE
      matrix_estimates <- matrix_estimates_list[[ifelse(point_estimator == "median", 4, 1)]]

      ## if only a single trial (i.e. a trial outcome) has been evaluated
      if (identical(nrow(matrix_estimates), 1L)) {

        estimates <- post_quantiles

      } else {

        ## find all estimates for response rates
        rr_index         <- grepl("p_", colnames(matrix_estimates))
        matrix_estimates <- matrix_estimates[, rr_index]

        point_estimates  <- as.matrix(colMeans(matrix_estimates))
        var_estimates    <- as.matrix(apply(matrix_estimates, 2, stats::var))

        # if (isTRUE(all.equal(dim(point_estimates), c(1, 1)))) {
        #
        #   point_estimates <- matrix(rep(point_estimates, length(true_rr)), ncol = 1)
        #   var_estimates   <- matrix(rep(var_estimates, length(true_rr)), ncol = 1)
        #
        # }

        bias_estimates <- point_estimates - t(true_rr)
        mse_estimates  <- bias_estimates^2 + var_estimates

        colnames(bias_estimates) <- "Bias"
        colnames(mse_estimates)  <- "MSE"

        ## Combine results
        ## introduce NAs for all values that are not response rates
        na_matrix      <- matrix(NA, nrow = sum(!rr_index), ncol = 1)

        bias_estimates <- rbind(bias_estimates, na_matrix)
        mse_estimates  <- rbind(mse_estimates, na_matrix)

        estimates <- cbind(post_quantiles, bias_estimates, mse_estimates)

      }

      ## Save results
      results_per_method_list[[method_names[n]]] <- estimates

    }

    if (!is.null(add_parameters) && all(!occurences)) {
      stop (simpleError(paste(
        "The additional parameters provided in 'add_parameters' do not occur",
        "in any of the methods stored in 'analyses_list'"
      )))
    }

    results_list[[s]] <- results_per_method_list

  }

  if (length(results_list) == 1) {

    results_list <- results_list[[1]]

  } else {

    results_list <- listPerMethod(results_list)

  }

  return (results_list)

}

getRespRatesEstimates <- function (

  analyses_list,
  cohort_names = NULL,
  alpha_level  = 0.05

) {

  gamma_levels <- matrix(rep(c(1 - alpha_level / 2, 0.5, alpha_level / 2),
                             each = length(cohort_names)), nrow = 3, byrow = TRUE)

  results_list <- vector(mode = "list", length = length(analyses_list))
  names(results_list) <- names(analyses_list)

  for (s in seq_along(analyses_list)) {

    analysis_data <- analyses_list[[s]]
    method_names  <- analysis_data$analysis_parameters$method_names

    results_per_method_list <- vector(mode = "list", length = length(method_names))
    names(results_per_method_list) <- method_names

    for (method_name in method_names) {

      estimates <- apply(gamma_levels, 1, function (glevel) {

        colMeans(do.call(rbind, getPosteriorGammaQuantiles(
          method_name              = method_name,
          gamma_levels             = glevel,
          quantiles                = analysis_data$analysis_parameters$quantiles,
          posterior_quantiles_list = analysis_data$quantiles_list,
          cohort_names             = cohort_names)))

      })

      colnames(estimates) <- paste0(c(alpha_level / 2, 0.5, 1 - alpha_level / 2) * 100, "%")

      results_per_method_list[[method_name]] <- estimates

    }

    results_list[[s]] <- results_per_method_list

  }

  results_list <- listPerMethod(results_list)

  return (results_list)

}

getAllCohortNames <- function (

  analyses_list

) {

  post_names <- colnames(analyses_list[[1]]$quantiles_list[[1]][[1]])
  indices    <- grep("p_", post_names)

  return (post_names[indices])

}

