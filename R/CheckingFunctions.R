### numeric ####

is.positive.numeric <- function (x) {

  return (is.numeric(x) && all(x > 0))

}

is.numeric.in.zero.one <- function (x) {

  return (is.numeric(x) && all(x > 0) && all(x < 1))

}

is.single.numeric <- function (x) {

  return (is.numeric(x) && length(x) == 1)

}

is.single.positive.numeric <- function (x) {

  return (is.single.numeric(x) && x > 0)

}

is.single.numeric.in.zero.one <- function (x) {

  return (is.single.numeric(x) && is.numeric.in.zero.one(x))

}

### whole number ####

is.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  if (is.numeric(x)) {

    return (abs(x - round(x)) < tol)

  } else {

    return (FALSE)

  }

}

is.positive.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && all(x > 0))

}

is.non.negative.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && all(x >= 0))

}

is.single.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (all(is.wholenumber(x, tol = tol)) && length(x) == 1)

}

is.single.non.negative.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (is.single.wholenumber(x, tol = tol) && all(is.non.negative.wholenumber(x, tol = tol)))

}

is.single.positive.wholenumber <- function (x, tol = .Machine$double.eps^0.5) {

  return (is.single.wholenumber(x, tol = tol) && all(is.positive.wholenumber(x, tol = tol)))

}


### evidence levels

check.evidence.levels <- function (evidence_levels,
                                   cohort_names, analyses_list, error_evidence_levels) {
  
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
