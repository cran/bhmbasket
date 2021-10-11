## R CMD check appeasement
utils::globalVariables("k")

#' @title logit
#' @description This function returns the logit of the input argument.
#' @param p A numeric in (0, 1)
#' @details This function is an alias for `stats::binomial()$linkfun`
#' @return logit of p
#' @rdname logit
#' @author Stephan Wojciekowski
#' @examples
#' logit(invLogit(0.3))
#' logit(c(0, 0.5, 1))
#' @seealso
#'  \code{\link[stats]{family}}
#' @export
logit <- function (p) {

  error_p <- simpleError("Please provide a numeric in (0, 1) for the argument 'p'")

  if (missing(p))                          stop (error_p)
  if (any(!is.numeric(p) | p < 0 | p > 1)) stop (error_p)

  stats::binomial()$linkfun(p)

}

#' @title invLogit
#' @description This function returns the inverse logit of the input argument.
#' @param theta A numeric
#' @details This function is an alias for `stats::binomial()$linkinv`
#' @return Inverse logit of theta
#' @rdname invlogit
#' @author Stephan Wojciekowski
#' @examples
#' invLogit(logit(0.3))
#' invLogit(c(-Inf, 0, Inf))
#' @seealso
#'  \code{\link[stats]{family}}
#' @export
invLogit <- function (theta) {

  error_theta <- simpleError("Please provide a numeric for the argument 'theta'")

  if (missing(theta))     stop (error_theta)
  if (!is.numeric(theta)) stop (error_theta)

  stats::binomial()$linkinv(theta)

}

firstUpper <- function (string) {

  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return (string)

}

cummulativeMovingAverage <- function (x) {
  return (cumsum(x) / seq_along(x))
}

getRowIndexOfVectorInMatrix <- function (

  vector_to_be_found,
  matrix_to_be_searched

) {

  n_col <- ncol(matrix_to_be_searched)

  if (length(vector_to_be_found) != n_col) {
    stop ("The length of the vector must be equal to the number of columns of the matrix")
  }

  index <- apply(convertVector2Matrix(sapply(seq_len(n_col), function (i) {

    matrix_to_be_searched[, i] %in% vector_to_be_found[i]

  })), 1, all)

  return (which(index))

}

substrRight <- function (x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

convertVector2Matrix <- function (vector) {

  if (is.null(dim(vector))) {
    vector <- t(as.matrix(vector))
  }

  return (vector)

}

#' @title scaleRoundList
#' @description This function applies scaling and rounding to each item of a list of numerics
#' @param list The list to which the scaling and rounding should be applied to.
#' @param scale_param A numeric for the scaling of each item of the list, Default: `1`
#' @param round_digits An integer for the number of digits.
#' If `NULL`, no rounding will be applied, Default: `NULL`
#' @return A list of scaled and rounded numerics
#' @rdname scaleRoundList
#' @author Stephan Wojciekowski
#' @examples
#' some_list <- as.list(runif(5))
#' scaleRoundList(some_list, scale_param = 100, round_digits = 2)
#'
#' scenarios_list <- simulateScenarios(
#'   n_subjects_list     = list(c(10, 20, 30)),
#'   response_rates_list = list(c(0.1, 0.2, 0.3)),
#'   n_trials            = 10)
#'
#' analyses_list <- performAnalyses(
#'   scenario_list       = scenarios_list,
#'   target_rates        = rep(0.5, 3),
#'   n_mcmc_iterations   = 100,
#'   n_cores             = 1L)
#'
#' scaleRoundList(
#'   list         = getEstimates(analyses_list),
#'   scale_param  = 100,
#'   round_digits = 2)
#' @export
scaleRoundList <- function(

  list,
  scale_param  = 1,
  round_digits = NULL

) {

  error_list  <- simpleError(
    "Please provide a list (of lists) of numerics for the argument 'list'")
  error_scale_param <- simpleError(
    "Please provide a positive numeric for the argument 'scale_param'")
  error_round_digits <- simpleError(
    "Please provide a positive integer for the argument 'round_digits'")

  if (missing(list))                                     stop (error_list)

  if (!is.list(list))                                    stop (error_list)
  if (!is.single.positive.numeric(scale_param))          stop (scale_param)
  if (!is.null(round_digits) &&
      !is.single.non.negative.wholenumber(round_digits)) stop (error_round_digits)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  list_levels <- getListLevel(list)

  out_list    <- roundList(scaleList(list, scale_param, list_levels),
                            round_digits, list_levels)

  return (out_list)

}

scaleList <- function (list, scale_param, list_levels) {

  scale_expression <- quote({
    if (is.numeric(a)) {
      a * scale_param
    } else {
      stop (simpleError("The list must contain numerics"))
    }
  })

  if (list_levels == 1) {

    return (lapply(list, function (a) eval(scale_expression)))

  } else if (list_levels == 2) {

    return (lapply(list, function (x)
      lapply(x, function (a) eval(scale_expression))))

  } else if (list_levels == 3) {

    return (lapply(list, function (x)
      lapply(x, function (y)
        lapply(y, function (a) eval(scale_expression)))))

  } else {

    stop ("lists with a nested depth greater than 3 are not supported")

  }

}

roundList <- function (list, round_digits, list_levels) {

  if (is.null(round_digits)) return (list)

  round_expression <- quote({
    if (is.numeric(a)) {
      round(a, round_digits)
    } else {
      stop (simpleError("The list must contain numerics"))
    }
  })

  if (list_levels == 1) {

    return (lapply(list, function (a) eval(round_expression)))

  } else if (list_levels == 2) {

    return (lapply(list, function (x)
      lapply(x, function (a) eval(round_expression))))

  } else if (list_levels == 3) {

    return (lapply(list, function (x)
      lapply(x, function (y)
        lapply(y, function (a) eval(round_expression)))))

  } else {

    stop ("lists with a nested depth greater than 3 are not supported")

  }

}

getListLevel <- function (

  list

) {

  if (!is.list(list)) {

    return (0)

  } else {

    list <- list[[1]]

    return (1 + getListLevel(list))

  }

}


listPerMethod <- function (

  list_per_scenario

) {

  all_method_names <- unique(as.vector(sapply(list_per_scenario, names)))

  ## Create list to hold output
  out_list <- vector(mode = "list", length = length(all_method_names))
  names(out_list) <- all_method_names

  for (n in seq_along(out_list)) {

    out_list[[n]] <- vector(mode = "list", length = length(list_per_scenario))
    names(out_list[[n]]) <- paste0("scenario_", getScenarioNumbers(list_per_scenario))

  }

  ## Copy contents to new output list
  for (name in all_method_names) {

    for (s in seq_along(list_per_scenario)) {

      out_list[[name]][[s]] <- list_per_scenario[[s]][[name]]

    }

  }

  return (out_list)

}
