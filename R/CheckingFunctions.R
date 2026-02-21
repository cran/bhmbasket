### evidence levels
check.evidence.levels <- function (
    
  evidence_levels,
  cohort_names,
  analyses_list,
  error_evidence_levels
  
) {
  
  checkmate::assert_true(
    length(evidence_levels) == length(cohort_names),
    .var.name = "The 'evidence_levels' and the 'cohort_names' must have the same length"
  )
  
  
  if (is.character(evidence_levels)) {
    
    mean_index <- evidence_levels == "mean"
    
    evidence_levels_numeric <- tryCatch({
      as.numeric(evidence_levels[!mean_index])
    }, warning = function(w) w)
    
    checkmate::assert_false(
      inherits(evidence_levels_numeric, "warning"),
      .var.name = "The only string allowed for the argument 'evidence_levels' is 'mean'"
    )
    
    
  } else {
    
    evidence_levels_numeric <- evidence_levels
    
  }
  
  checkmate::assert_numeric(
    evidence_levels_numeric,
    lower = 0, upper = 1, any.missing = FALSE,
    .var = error_evidence_levels
  )
  
  available_quantiles <- round(analyses_list[[1]]$analysis_parameters$quantiles, 9)
  asked_quantiles     <- round(1 - evidence_levels_numeric, 9)
  
  checkmate::assert_true(
    all(asked_quantiles %in% available_quantiles),
    .var.name = paste(
      "The 'evidence_levels' must have matches",
      "in the 'evidence_levels' provided to the call performAnalyses()",
      "that created the 'analyses_list'"
    )
  )
}

### functions ####

checkForParallelBackend <- function () {
  
  "%dopar%" <- foreach::"%dopar%"
  if(!foreach::getDoParRegistered()) {
    
    message("\nCaution: No parallel backend detected for the 'foreach' framework.",
            " For execution in parallel, register a parallel backend, e.g. with:\n",
            "   doFuture::registerDoFuture()\n",
            "   future::plan(future::multisession)\n")
    
    # foreach if without parallel backend gives a warning, warning is being
    # displayed only once, should be repeated anyway, message better
    
    tt <- suppressWarnings(foreach::foreach(k = 1:2) %dopar% {k^k^k})
    rm(tt)
    
  }
  
}