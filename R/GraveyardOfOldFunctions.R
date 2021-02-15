# ## RIP 10 Dec 2019
# ##  Calculates each scenario individually
# ##  Does not take advantage of some unique trials being in more than one scenario.
# ##  Can save scenarios within function
# OldPerformAnalyses <- function (
#
#   scenario.data.list,
#   quantiles                  = c(0.025, 0.1, 0.2, 0.5, 0.975),
#   post.mean                  = FALSE,
#
#   method.names               = c("berry", "exnex", "exnex.adj", "pooled", "stratified"),
#   target.rates.list,
#   prior.parameters.berry     = NULL,
#   prior.parameters.exnex     = NULL,
#   prior.parameters.exnex.adj = NULL,
#   tau.scale                  = 1,
#
#   calc.differences           = NULL,
#
#   n.iterations               = 1e4,
#   n.cores                    = parallel::detectCores() - 1L,
#   seed                       = 123,
#   save.path                  = NULL
#
# ) {
#
#   method.names <- sort(method.names)
#
#   cat(format(Sys.time(), "%m/%d/%y"), " Perform Analyses\n", sep = "")
#
#   ## Get scenario numbers
#   scenario.numbers <- sapply(scenario.data.list, function (x) x$scenario.number)
#
#   if (length(scenario.numbers) != length(target.rates.list)){
#     stop ('"scenario.numbers", "target.rates.list" must have the same length.')
#   }
#
#   ## For each scenario
#   scenario.analyses.list <- vector(mode = "list", length = length(scenario.numbers))
#   names(scenario.analyses.list) <- paste0("Scenario.", scenario.numbers)
#
#   for (s in seq_along(scenario.numbers)) {
#
#     cat("         Analyzing Scenario ", scenario.numbers[s], " ...\n", sep = "")
#
#     ## Get scenario data
#     scenario.data <- scenario.data.list[[s]]
#
#     if (ncol(scenario.data$n.subjects.unique) != length(target.rates.list[[s]])) {
#
#       stop (paste0('The lengths of the target.rates vector must be',
#                    'equal to the number of cohorts of the scenario.'))
#
#     }
#
#     ## Create folder to save results
#     if (!is.null(save.path)) {
#
#       analysis.path <- CreateAnalysisPath(
#         save.path       = save.path,
#         scenario.number = scenario.numbers[s],
#         analysis.number = 0)
#
#     } else {
#
#       analysis.path <- NULL
#
#     }
#
#     ## Run analysis
#     method.quantiles.list            <- vector(mode = "list", length = length(method.names))
#     analysis.paramerters.list        <- vector(mode = "list", length = length(method.names))
#     names(method.quantiles.list)     <- method.names
#     names(analysis.paramerters.list) <- method.names
#
#     ## For each method
#     for (method.name in method.names) {
#
#       start.time  <- Sys.time()
#       out.message <- paste0(format(start.time, "   %H:%M", digits = 1),
#                             " - with ", FirstUpper(method.name), " ...")
#       cat(out.message, rep(".", 33 - nchar(out.message)), sep = "")
#
#       analysis.paramerters.list[[method.name]] <- PrepareAnalyses(
#         method.name                = method.name,
#         n.cohorts                  = ncol(scenario.data$n.subjects.unique),
#         target.rates               = target.rates.list[[s]],
#         prior.parameters.berry     = prior.parameters.berry,
#         prior.parameters.exnex     = prior.parameters.exnex,
#         prior.parameters.exnex.adj = prior.parameters.exnex.adj,
#         tau.scale                  = tau.scale)
#
#       method.quantiles.list[[method.name]] <- GetPostQuantiles(
#         method.name        = method.name,
#         quantiles          = quantiles,
#         post.mean          = post.mean,
#         scenario.data      = scenario.data,
#         calc.differences   = calc.differences,
#         j.parameters       = analysis.paramerters.list[[method.name]]$j.parameters,
#         j.model.file       = analysis.paramerters.list[[method.name]]$j.model.file,
#         j.data             = analysis.paramerters.list[[method.name]]$j.data,
#         n.iterations       = n.iterations,
#         seed               = seed,
#         n.cores            = n.cores,
#         save.path          = analysis.path,
#         save.trial         = ifelse(is.null(analysis.path), NULL,
#                                     which.max(scenario.data$weights.unique)))
#
#       if (method.name %in% c("berry", "exnex", "exnex.adj")) {
#         try ({PlotPostRR(method.name     = method.name,
#                          scenario.number = scenario.numbers[s],
#                          analysis.path   = analysis.path,
#                          plot.median     = "none")})
#       }
#
#       cat(" Finished after ", round(Sys.time() - start.time, 1), " ",
#           units(Sys.time() - start.time), ".\n", sep = "")
#       rm(start.time)
#
#     }
#
#     scenario.analyses.list[[s]] <- list(
#       quantiles.list      = method.quantiles.list,
#       scenario.data       = scenario.data,
#       analysis.parameters = list(
#         quantiles         = quantiles,
#         method.names      = method.names,
#         prior.parameters  = analysis.paramerters.list))
#
#     if (!is.null(save.path)) {
#
#       analysis.list <- SaveAnalyses(
#         scenario.analyses.list = list(scenario.analyses.list[[s]]),
#         save.path              = save.path,
#         analysis.numbers       = as.numeric(substr(analysis.path,
#                                                    nchar(analysis.path),
#                                                    nchar(analysis.path))))
#
#     }
#
#   }
#
#   names(scenario.analyses.list) <- paste0("scenario.", scenario.numbers)
#
#   return (scenario.analyses.list)
#
# }
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 10 Dec 2019
# ##  This function incorporates the possibility for a simple interim analysis and
# ##  returns the unique trials and their respective weights.
# OldGetUniqueTrials <- function (
#
#   n.subjects,
#   response.rates,
#
#   interim.n.subjects = NULL,
#   interim.n.min      = NULL,
#
#   n.trials = 1e4,
#
#   seed     = 123,
#   n.cores  = parallel::detectCores() - 1L
#
# ) {
#
#   "%dopar%" <- foreach::"%dopar%"
#
#   set.seed(seed)
#
#   if (is.null(dim(response.rates))) {
#     response.rates <- t(as.matrix(response.rates))
#   }
#
#   colnames(response.rates) <- paste0("rr.cohort", seq_len(ncol(response.rates)))
#
#   if (length(n.subjects) != ncol(response.rates)) {
#     stop ("n.subjects and response.rates must have same length")
#   }
#
#   if (any(response.rates < 1 & response.rates > 0)) {
#
#     ## Simulations for new cohorts
#     ## A cohort is new if it has a response rate greater than 0 and less than 1.
#     ## Response rates greater than or equal to one will be used as fixed responses.
#
#     new.cohorts <- TRUE
#     index.new   <- which(response.rates < 1 & response.rates > 0)
#
#     if (all(!is.null(interim.n.subjects), !is.null(interim.n.min))) {
#
#       ## Both interim.n.subjects and interim.n.min are specified
#
#       if (length(interim.n.subjects) != length(interim.n.min)) {
#         stop ("The length of interim.n.subjects and interim.n.min must be equal.")
#       }
#
#       if (length(interim.n.subjects) != length(n.subjects[index.new])) {
#         stop ("The length of interim.n.subjects and n.subjects must be equal.")
#       }
#
#       interim.analysis <- TRUE
#
#       ## Simulate both the data before and after interim analysis
#       pre.interim.n.responders <- GetRespondersParallel(
#         response.rates = response.rates[, index.new],
#         n.subjects     = interim.n.subjects,
#         n.trials       = n.trials,
#         n.cores        = n.cores,
#         seed           = seed)
#       post.interim.n.responders <- GetRespondersParallel(
#         response.rates = response.rates[, index.new],
#         n.subjects     = n.subjects[index.new] - interim.n.subjects,
#         n.trials       = n.trials,
#         n.cores        = n.cores,
#         seed           = seed + n.trials)
#
#       ## Check which cohorts passed the interim analysis
#       interim.passed <- GetInterimPassed(interim.n.responders = pre.interim.n.responders,
#                                          interim.n.min        = interim.n.min)
#
#       ## Combine the data before and after analysis according to
#       ## whether the cohort passed the interim analysis
#       n.responders <- pre.interim.n.responders + interim.passed * post.interim.n.responders
#
#       ## Get number of responders
#       pre.interim.n.subjets  <- matrix(interim.n.subjects,
#                                        ncol = length(interim.n.subjects),
#                                        nrow = n.trials, byrow = TRUE)
#       post.interim.n.subjets <- matrix(n.subjects[index.new] - interim.n.subjects,
#                                        ncol = length(interim.n.subjects),
#                                        nrow = n.trials, byrow = TRUE)
#       n.subjects.new <- pre.interim.n.subjets + interim.passed * post.interim.n.subjets
#
#       ## Calculate interim go probabilities
#       interim.passed   <- cbind(overall = apply(interim.passed, 1, any), interim.passed)
#       interim.go.probs <- colMeans(interim.passed)
#
#     } else if (any(!is.null(interim.n.subjects), !is.null(interim.n.min))) {
#
#       ## Only one of interim.n.subjects and interim.n.min are specified
#
#       stop (paste0("If one of interim.n.subjects or interim.n.min is not NULL, ",
#                    "the other one must also be not NULL."))
#
#     } else {
#
#       ## No interim analysis as none of interim.n.subjects and interim.n.min are specified
#
#       interim.analysis <- FALSE
#
#       ## Simulate the number of responses without interim analysis
#       n.responders <- GetRespondersParallel(
#         response.rates = response.rates[, index.new],
#         n.subjects     = n.subjects[index.new],
#         n.trials       = n.trials,
#         n.cores        = n.cores,
#         seed           = seed)
#
#       n.subjects.new <- matrix(n.subjects[index.new],
#                                ncol = length(n.subjects[index.new]),
#                                nrow = n.trials, byrow = TRUE)
#
#     }
#
#     ## Get the unique trials
#     unique.trials       <- MakeUnique(cbind(n.responders, n.subjects.new))
#     n.responders.unique <- unique.trials[[1]][, seq_len(ncol(n.responders))]
#     n.subjects.unique   <- unique.trials[[1]][, seq_len(ncol(n.subjects.new)) + ncol(n.responders)]
#     weights.unique      <- unique.trials[[2]]
#
#   } else {
#
#     ## No new cohorts, as response rates are all greater than or equal to 1
#
#     new.cohorts <- FALSE
#
#   }
#
#   if (any(response.rates >= 1 | response.rates == 0)) {
#
#     ## Historical cohorts
#
#     hist.cohorts <- TRUE
#     index.hist   <- which(response.rates >= 1 | response.rates == 0)
#
#     if (new.cohorts) {
#
#       ## In case of simulated (new) cohorts, the historic responses will be recycled to match
#       ## the number of unique trials
#
#       ## Make matrix if necessary
#       if (is.null(nrow(n.responders.unique))) {
#         n.responders.unique = as.matrix(n.responders.unique,
#                                         nrow = 1, ncol = length(n.responders.unique))
#       }
#       if (is.null(nrow(n.subjects.unique))) {
#         n.subjects.unique = as.matrix(n.subjects.unique,
#                                       nrow = 1, ncol = length(n.responders.unique))
#       }
#
#       historic.responses <- matrix(rep(response.rates[index.hist],
#                                        each = nrow(n.responders.unique)),
#                                    nrow = nrow(n.responders.unique))
#
#       historic.subjects  <- matrix(rep(n.subjects[index.hist],
#                                        each = nrow(n.subjects.unique)),
#                                    nrow = nrow(n.subjects.unique))
#
#     } else {
#
#       ## In case of no simulated cohorts, only one set of fixed responses is needed
#
#       historic.responses <- matrix(response.rates[index.hist], nrow = 1)
#       historic.subjects  <- matrix(n.subjects[index.hist], nrow = 1)
#
#     }
#
#   } else {
#
#     ## No historical cohorts, as response rates are all smaller than 1
#
#     hist.cohorts <- FALSE
#
#   }
#
#   ## Combine new and historical cohorts as appropriate
#   if (new.cohorts & hist.cohorts) {
#
#     n.responders.unique <- cbind(n.responders.unique, historic.responses)
#     n.subjects.unique   <- cbind(n.subjects.unique, historic.subjects)
#
#   } else if (hist.cohorts) {
#
#     n.responders.unique <- historic.responses
#     n.subjects.unique   <- historic.subjects
#     weights.unique      <- 1
#
#   }
#
#   ## Create list to return data
#   scenario.data <- list(n.subjects.unique   = n.subjects.unique,
#                         n.responders.unique = n.responders.unique,
#                         weights.unique      = weights.unique,
#                         response.rates      = response.rates,
#                         n.subjects          = n.subjects,
#                         n.trials            = n.trials)
#
#   if (new.cohorts) {
#     if (interim.analysis) {
#       scenario.data$interim.go.probs <- interim.go.probs
#     }
#   }
#
#   return (scenario.data)
#
# }
#
# ## RIP 10 Dec 2019
# ##  Belongs to the function above
# ##  Can save scenarios within the function
# OldSimulateScenarios <- function (
#
#   n.subjects.list,
#   response.rates.list,
#
#   interim.n.subjects.list  = NULL,
#   interim.n.min.list       = NULL,
#
#   scenario.numbers         = NULL,
#
#   n.trials                 = 1e4,
#   seed                     = 123,
#   n.cores                  = parallel::detectCores() - 1L,
#
#   save.path                = NULL
#
# ) {
#
#   if (is.null(scenario.numbers)) {
#
#     scenario.numbers <- seq_along(response.rates.list)
#
#   }
#
#   if (length(n.subjects.list)     != length(response.rates.list) |
#       length(response.rates.list) != length(scenario.numbers) |
#       length(n.subjects.list)     != length(scenario.numbers)) {
#
#     stop (paste0("n.subjects.list, response.rates.list, ",
#                  "and scenario.numbers must have same length"))
#   }
#
#   cat(format(Sys.time(), "%m/%d/%y"), " Simulating Scenarios\n", sep = "")
#
#   scenario.data.list <- vector(mode = "list", length = length(scenario.numbers))
#   for (s in seq_along(scenario.numbers)) {
#
#     start.time <- Sys.time()
#     out.message <- paste0(format(start.time, "   %H:%M"),
#                           " Simulating Scenario ", scenario.numbers[s], " ...")
#     cat(out.message, rep(".", 33 + nchar(max(scenario.numbers)) - nchar(out.message)), sep = "")
#
#     scenario.data.list[[s]] <- GetUniqueTrials(
#       n.subjects          = n.subjects.list[[s]],
#       response.rates      = response.rates.list[[s]],
#
#       interim.n.subjects  = interim.n.subjects.list[[s]],
#       interim.n.min       = interim.n.min.list[[s]],
#
#       n.trials            = n.trials,
#
#       seed                = seed,
#       n.cores             = n.cores
#     )
#
#     scenario.data.list[[s]]$scenario.number <- scenario.numbers[s]
#
#     if (!is.null(save.path)) {
#
#       SaveScenarios(
#         scenario.data.list = list(scenario.data.list[[s]]),
#         save.path          = save.path)
#
#     }
#
#     cat(" Finished after ", round(Sys.time() - start.time, 1), " ",
#         units(Sys.time() - start.time), ".\n", sep = "")
#     rm(start.time)
#
#   }
#
#   names(scenario.data.list) <- paste0("scenario.", scenario.numbers)
#
#   return (scenario.data.list)
#
# }
#
# ## RIP 10 Dec 2019
# ##  This function belongs to the functions above
# ##  It incorporates unique trials with their respective weights.
# GetGoProbabilities <- function (
#
#   scenario.analyses.list,
#
#   cohort.names,
#   decision.rules.list,
#   gamma.levels.list,
#
#   digits = NULL
#
# ) {
#
#   ## Get scenario numbers
#   scenario.numbers <- as.numeric(sub("scenario.", "", names(scenario.analyses.list)))
#
#   ## Get method names
#   method.names.matrix <- t(sapply(scenario.analyses.list,
#                                   function (x) x$analysis.parameters$method.names))
#   if (!all(sapply(1:nrow(method.names.matrix),
#                   function (x) identical(method.names.matrix[1, ], method.names.matrix[x, ])))) {
#     stop ("The scenarios where analysed with different methods.")
#   }
#   method.names <- method.names.matrix[1, ]
#
#   ## in case only one method has been used for all scenarios
#   if (all(sapply(seq_along(method.names), function (x) {
#     method.names[1] == method.names[x]
#   }))) {
#     method.names <- method.names[1]
#   }
#
#   ## check for input consistency
#   if (length(method.names) != length(decision.rules.list)) {
#     stop ("The lengths of method.names and decision.rules.list must be equal.")
#   }
#   if (length(method.names) != length(gamma.levels.list)) {
#     stop ("The lengths of method.names and gamma.levels.list must be equal.")
#   }
#
#   matrices.list        <- vector(mode = "list", length = length(method.names))
#   names(matrices.list) <- method.names
#
#   for (method.index in seq_along(method.names)) {
#
#     method.name   <- method.names[method.index]
#     decision.rule <- decision.rules.list[[method.index]]
#     gamma.levels  <- gamma.levels.list[[method.index]]
#
#     results.list        <- vector(mode = "list", length = length(scenario.numbers))
#     names(results.list) <- scenario.numbers
#
#     for (scenario.number in scenario.numbers) {
#
#       analysis.data <- scenario.analyses.list[[scenario.number]]
#
#       go.decisions <- GetGoDecisionsByCohort(
#         gamma.quantiles = GetPosteriorGammaQuantiles(
#           method.name              = method.name,
#           gamma.indices            = GetGammaIndices(
#             gamma.levels = gamma.levels,
#             quantiles    = analysis.data$analysis.parameters$quantiles),
#           posterior.quantiles.list = analysis.data$quantiles.list,
#           cohort.names             = cohort.names),
#         decision.rule   = decision.rule)
#
#       if (ncol(go.decisions) != 1) {
#
#         go.decisions <- cbind(overall = apply(go.decisions, 1, any), go.decisions)
#
#       } else {
#
#         colnames(go.decisions) <- "overall"
#
#       }
#
#       go.probs <- colSums(apply(go.decisions, 2, function(x) {
#         x * analysis.data$scenario.data$weights.unique
#       }))
#
#       response.rates <- t(mapply(function (x, y) ifelse(x < 1, x, x / y),
#                                  x = analysis.data$scenario.data$response.rates,
#                                  y = analysis.data$scenario.data$n.subjects))
#       colnames(response.rates) <- colnames(analysis.data$scenario.data$response.rates)
#
#       results.list[[scenario.number]] <- cbind(t(as.matrix(go.probs)), response.rates)
#
#     }
#
#     results.matrix <- do.call(rbind, results.list)
#     rownames(results.matrix) <- paste0("Sc.", scenario.numbers)
#
#     matrices.list[[method.name]] <- results.matrix
#
#   }
#
#   if (!is.null(digits)) {
#
#     matrices.list <- lapply(matrices.list, function (x) round(x, digits))
#
#   }
#
#   return (matrices.list)
#
# }
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 11 Dec 2019
# ##  Does the Go probabilities directly
# ##  Replaced by GetGoDecisions & an adapted GetGoProbabilities
# OldGetGoProbabilities <- function (
#
#   scenario.analyses.list,
#
#   cohort.names,
#   decision.rules.list,
#   gamma.levels.list,
#
#   digits = NULL
#
# ) {
#
#   ## Get scenario numbers
#   scenario.numbers <- as.numeric(sub("scenario.", "", names(scenario.analyses.list)))
#
#   ## Get method names
#   method.names.matrix <- t(sapply(scenario.analyses.list,
#                                   function (x) x$analysis.parameters$method.names))
#   if (!all(sapply(1:nrow(method.names.matrix),
#                   function (x) identical(method.names.matrix[1, ], method.names.matrix[x, ])))) {
#     stop ("The scenarios where analysed with different methods.")
#   }
#   method.names <- method.names.matrix[1, ]
#
#   ## in case only one method has been used for all scenarios
#   if (all(sapply(seq_along(method.names), function (x) {
#     method.names[1] == method.names[x]
#   }))) {
#     method.names <- method.names[1]
#   }
#
#   ## check for input consistency
#   gamma.levels.list <- quote(list(0.8, 0.9, 0.5, "mean"))
#
#   if (!is.list(decision.rules.list)) {
#     decision.rules.list <- list(decision.rules.list)
#   }
#   if (!is.list(gamma.levels.list)) {
#     gamma.levels.list <- list(gamma.levels.list)
#   }
#
#   if (length(method.names) < length(decision.rules.list)) {
#     stop ("The lengths of method.names must be greater than the length of decision.rules.list")
#   } else if (length(method.names) > length(decision.rules.list)) {
#     decision.rules.list <- rep(decision.rules.list, length.out = length(method.names))
#   }
#   if (length(method.names) < length(gamma.levels.list)) {
#     stop ("The lengths of method.names must be greater than the length of gamma.levels.list")
#   } else if (length(method.names) > length(gamma.levels.list)) {
#     gamma.levels.list <- rep(gamma.levels.list, length.out = length(method.names))
#   }
#
#   matrices.list        <- vector(mode = "list", length = length(method.names))
#   names(matrices.list) <- method.names
#
#   for (method.index in seq_along(method.names)) {
#
#     method.name   <- method.names[method.index]
#     decision.rule <- decision.rules.list[[method.index]]
#     gamma.levels  <- eval(gamma.levels.list[[method.index]])
#
#     results.list        <- vector(mode = "list", length = length(scenario.numbers))
#     names(results.list) <- scenario.numbers
#
#     for (scenario.number in scenario.numbers) {
#
#       analysis.data <- scenario.analyses.list[[scenario.number]]
#
#       go.decisions <- GetGoDecisionsByCohort(
#         gamma.quantiles = GetPosteriorGammaQuantiles(
#           method.name              = method.name,
#           gamma.indices            = GetGammaIndices(
#             gamma.levels = gamma.levels,
#             quantiles    = analysis.data$analysis.parameters$quantiles),
#           posterior.quantiles.list = analysis.data$quantiles.list,
#           cohort.names             = cohort.names),
#         decision.rule   = decision.rule)
#
#       if (ncol(go.decisions) != 1) {
#
#         go.decisions <- cbind(overall = apply(go.decisions, 1, any), go.decisions)
#
#       } else {
#
#         colnames(go.decisions) <- "overall"
#
#       }
#
#       go.probs <- colMeans(go.decisions)
#
#       response.rates <- analysis.data$scenario.data$response.rates
#
#       if (any(response.rates == 0 | response.rates > 1)) {
#
#         indices <- which(response.rates == 0 | response.rates > 1)
#
#         response.rates[indices] <- response.rates[indices] /
#           analysis.data$scenario.data$n.subjects[1, indices]
#
#       }
#
#       colnames(response.rates) <- colnames(analysis.data$scenario.data$response.rates)
#
#       results.list[[scenario.number]] <- cbind(t(as.matrix(go.probs)), response.rates)
#
#     }
#
#     results.matrix <- do.call(rbind, results.list)
#     rownames(results.matrix) <- paste0("Sc.", scenario.numbers)
#
#     matrices.list[[method.name]] <- results.matrix
#
#   }
#
#   if (!is.null(digits)) {
#
#     matrices.list <- lapply(matrices.list, function (x) round(x, digits))
#
#   }
#
#   return (matrices.list)
#
# }
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 12 Dec 2019
# ##  To be splitted in GetTrials and Get Scenarios
# ##  This function returns a scenario list
# ##  The new GetTrials will only return the simulated trials.
# GetTrials <- function (
#
#   n.subjects,
#   response.rates,
#
#   cohort.names = paste0("cohort.", seq_along(n.subjects)),
#
#   n.trials     = 1e4,
#
#   seed         = 123,
#   n.cores      = parallel::detectCores() - 1L
#
# ) {
#
#   "%dopar%" <- foreach::"%dopar%"
#
#   if (length(n.subjects) != length(response.rates)) {
#     stop ("n.subjects and response.rates must have same length")
#   }
#
#   set.seed(seed)
#
#   response.rates           <- ConvertVector2Matrix(response.rates)
#   colnames(response.rates) <- paste0("rr.", cohort.names)
#
#   if (any(response.rates < 1 & response.rates > 0)) {
#
#     ## Simulations for new cohorts
#     ## A cohort is new if it has a response rate greater than 0 and less than 1.
#     ## Response rates greater than or equal to one will be used as fixed responses.
#
#     new.cohorts <- TRUE
#     index.new   <- which(response.rates < 1 & response.rates > 0)
#
#     ## Simulate the number of responses without interim analysis
#     n.responders <- GetRespondersParallel(
#       response.rates = response.rates[, index.new],
#       n.subjects     = n.subjects[index.new],
#       n.trials       = n.trials,
#       n.cores        = n.cores,
#       seed           = seed)
#
#   } else {
#
#     ## No new cohorts, as response rates are all greater than or equal to 1
#
#     new.cohorts <- FALSE
#
#   }
#
#   if (any(response.rates >= 1 | response.rates == 0)) {
#
#     ## Historical cohorts
#
#     hist.cohorts <- TRUE
#     index.hist   <- which(response.rates >= 1 | response.rates == 0)
#
#     if (new.cohorts) {
#
#       ## In case of simulated (new) cohorts, the historic responses will be recycled to match
#       ## the number of unique trials
#
#       ## Make matrix if necessary
#       n.responders.hist <-  matrix(rep(response.rates[index.hist],
#                                        each = nrow(n.responders)),
#                                    nrow = nrow(n.responders))
#
#     } else {
#
#       ## In case of no simulated cohorts, only one set of fixed responses is needed
#
#       n.responders.hist <- matrix(response.rates[index.hist], nrow = 1)
#
#     }
#
#   } else {
#
#     ## No historical cohorts, as response rates are all smaller than 1
#
#     hist.cohorts <- FALSE
#
#   }
#
#   ## Combine new and historical cohorts as appropriate
#   if (new.cohorts & hist.cohorts) {
#
#     n.responders <- cbind(n.responders, n.responders.hist)
#
#   } else if (hist.cohorts) {
#
#     n.responders <- n.responders.hist
#
#   }
#
#   n.subjects <- matrix(n.subjects,  byrow = TRUE,
#                        ncol = length(n.subjects), nrow = nrow(n.responders))
#
#   colnames(n.subjects)   <- paste0("n.", cohort.names)
#   colnames(n.responders) <- paste0("r.", cohort.names)
#
#   ## Create list to return data
#   scenario.data <- list(n.subjects     = n.subjects,
#                         n.responders   = n.responders,
#                         response.rates = response.rates,
#                         n.trials       = n.trials,
#                         seed           = seed)
#
#   return (scenario.data)
#
# }
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 26 Dec 2019
# ##  This function is updated. The processing part of the scenarios
# ##  is run in parallel in the new function
# OldPerformAnalyses <- function (
#
#   scenario.data.list,
#   quantiles                  = c(0.025, 0.1, 0.2, 0.5, 0.975),
#   post.mean                  = FALSE,
#
#   method.names               = c("berry", "exnex", "exnex.adj", "pooled", "stratified"),
#   target.rates.list          = NULL,
#   prior.parameters.berry     = NULL,
#   prior.parameters.exnex     = NULL,
#   prior.parameters.exnex.adj = NULL,
#   tau.scale                  = 1,
#
#   calc.differences           = NULL,
#
#   n.iterations               = 1e4,
#   n.cores                    = parallel::detectCores() - 1L,
#   seed                       = 123
#
# ) {
#
#   ## Saving in from within this function no longer sensible
#   save.path                  = NULL
#
#   method.names <- sort(method.names)
#
#   cat(format(Sys.time(), "%m/%d/%y"), " Performing Analyses\n", sep = "")
#
#   ## Get scenario numbers
#   scenario.numbers <- sapply(scenario.data.list, function (x) x$scenario.number)
#
#   ## Exception for ExNex when prior parameters have been specified
#   if (length(method.names) != length(target.rates.list)) {
#     if (!(method.names == "exnex" && !is.null(prior.parameters.exnex))) {
#       stop ('"method.names" and "target.rates.list" must have the same length.')
#     }
#   }
#
#   ## Get unique trials for all scenarios
#   all.scenarios.n.responders <- do.call(rbind, lapply(scenario.data.list,
#                                                       function (x) x$n.responders))
#   all.scenarios.n.subjects   <- do.call(rbind, lapply(scenario.data.list,
#                                                       function (x) x$n.subjects))
#
#   n.cohorts     <- ncol(all.scenarios.n.responders)
#   trials.unique <- MakeUnique(cbind(all.scenarios.n.responders, all.scenarios.n.subjects))
#   n.responders  <- trials.unique[[1]][, seq_len(n.cohorts)]
#   n.subjects    <- trials.unique[[1]][, seq_len(n.cohorts) + n.cohorts]
#
#   scenarios.message <- paste0(as.character(scenario.numbers), sep = ", ", collapse = "")
#   cat("         Analyzing Scenarios ", substr(scenarios.message, 1, nchar(scenarios.message) - 2),
#       " ...\n", sep = "")
#
#   ## Create lists to save all results of the unique trials
#   method.quantiles.list            <- vector(mode = "list", length = length(method.names))
#   analysis.paramerters.list        <- vector(mode = "list", length = length(method.names))
#   names(method.quantiles.list)     <- method.names
#   names(analysis.paramerters.list) <- method.names
#
#   ## For each method
#   for (n in seq_along(method.names)) {
#
#     start.time  <- Sys.time()
#     out.message <- paste0(format(start.time, "   %H:%M", digits = 1),
#                           " - with ", FirstUpper(method.names[n]), " ...")
#     cat(out.message, rep(".", 33 - nchar(out.message)), sep = "")
#
#     analysis.paramerters.list[[method.names[n]]] <- PrepareAnalyses(
#       method.name                = method.names[n],
#       n.cohorts                  = n.cohorts,
#       target.rates               = target.rates.list[[n]],
#       prior.parameters.berry     = prior.parameters.berry,
#       prior.parameters.exnex     = prior.parameters.exnex,
#       prior.parameters.exnex.adj = prior.parameters.exnex.adj,
#       tau.scale                  = tau.scale)
#
#     method.quantiles.list[[method.names[n]]] <- GetPostQuantiles(
#       method.name        = method.names[n],
#       quantiles          = quantiles,
#       post.mean          = post.mean,
#       scenario.data      = list(n.subjects   = n.subjects,
#                                 n.responders = n.responders),
#       calc.differences   = calc.differences,
#       j.parameters       = analysis.paramerters.list[[method.names[n]]]$j.parameters,
#       j.model.file       = analysis.paramerters.list[[method.names[n]]]$j.model.file,
#       j.data             = analysis.paramerters.list[[method.names[n]]]$j.data,
#       n.iterations       = n.iterations,
#       seed               = seed,
#       n.cores            = n.cores,
#       save.path          = NULL,
#       save.trial         = NULL)
#
#     cat(" Finished after ", round(Sys.time() - start.time, 1), " ",
#         units(Sys.time() - start.time), ".\n", sep = "")
#     rm(start.time)
#
#   }
#
#   ## For each scenario
#   scenario.analyses.list <- vector(mode = "list", length = length(scenario.numbers))
#   names(scenario.analyses.list) <- paste0("Scenario.", scenario.numbers)
#
#   for (s in seq_along(scenario.numbers)) {
#
#     start.time  <- Sys.time()
#     out.message <- paste0("         Processing Scenario ", scenario.numbers[s], " ")
#     cat(out.message, rep(".", abs(33 - nchar(out.message))), sep = "")
#
#     ## Get scenario data
#     scenario.data <- scenario.data.list[[s]]
#
#     ## Create folder to save results
#     if (!is.null(save.path)) {
#
#       analysis.path <- CreateAnalysisPath(
#         save.path       = save.path,
#         scenario.number = scenario.numbers[s],
#         analysis.number = 0)
#
#     } else {
#
#       analysis.path <- NULL
#
#     }
#
#     ## Find the indices of the trials of a specific scenario
#     unique.trials.matrix <- cbind(n.responders, n.subjects)
#     scenario.data.matrix <- cbind(scenario.data$n.responders,
#                                   scenario.data$n.subjects)
#     trial.indices <- sapply(seq_len(nrow(scenario.data$n.responders)), function (i) {
#       GetRowIndexOfVectorInMatrix(
#         vector.to.be.found    = scenario.data.matrix[i, ],
#         matrix.to.be.searched = unique.trials.matrix)
#     })
#
#     ## Save scenario specific posterior quantiles for all methods
#     scenario.method.quantiles <- vector(mode = "list", length = length(method.names))
#     names(scenario.method.quantiles) <- method.names
#
#     for (n in seq_along(method.names)) {
#
#       # ## Target rates not required for ExNex
#       # if (ncol(scenario.data$n.subjects) != length(target.rates.list[[n]])) {
#       #
#       #   stop (paste0('The lengths of the target.rates vector must be',
#       #                'equal to the number of cohorts of the scenario.'))
#       #
#       # }
#
#       scenario.method.quantiles[[method.names[n]]] <- lapply(trial.indices, function (i) {
#         method.quantiles.list[[method.names[n]]][[i]]
#       })
#
#     }
#
#     scenario.analyses.list[[s]] <- list(
#       quantiles.list      = scenario.method.quantiles,
#       scenario.data       = scenario.data,
#       analysis.parameters = list(
#         quantiles         = quantiles,
#         method.names      = method.names,
#         prior.parameters  = analysis.paramerters.list))
#
#     if (!is.null(analysis.path)) {
#
#       analysis.list <- SaveAnalyses(
#         scenario.analyses.list = list(scenario.analyses.list[[s]]),
#         save.path              = save.path,
#         analysis.numbers       = as.numeric(substr(analysis.path,
#                                                    nchar(analysis.path),
#                                                    nchar(analysis.path))))
#
#     }
#
#     cat(" Finished after ", round(Sys.time() - start.time, 1), " ",
#         units(Sys.time() - start.time), ".\n", sep = "")
#     rm(start.time)
#
#   }
#
#   names(scenario.analyses.list) <- paste0("scenario.", scenario.numbers)
#
#   return (scenario.analyses.list)
#
# }
#
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 24 Jan 2020
# ##  This function has been updated to display the results in a nicer table
# GetRespRatesEstimates <- function (
#
#   scenario.analyses.list,
#   cohort.names,
#   alpha.level = 0.05
#
# ) {
#
#   gamma.levels <- matrix(rep(c(1 - alpha.level / 2, 0.5, alpha.level / 2),
#                              each = length(cohort.names)), nrow = 3, byrow = TRUE)
#
#   results.list <- vector(mode = "list", length = length(scenario.analyses.list))
#   names(results.list) <- names(scenario.analyses.list)
#
#   for (s in seq_along(scenario.analyses.list)) {
#
#     analysis.data <- scenario.analyses.list[[s]]
#     method.names  <- analysis.data$analysis.parameters$method.names
#
#     results.per.method.list <- vector(mode = "list", length = length(method.names))
#     names(results.per.method.list) <- method.names
#
#     for (method.name in method.names) {
#
#       estimates <- matrix(t(apply(gamma.levels, 1, function (glevel) {
#
#         colMeans(do.call(rbind, GetPosteriorGammaQuantiles(
#           method.name              = method.name,
#           gamma.levels             = glevel,
#           quantiles                = analysis.data$analysis.parameters$quantiles,
#           posterior.quantiles.list = analysis.data$quantiles.list,
#           cohort.names             = cohort.names)))
#
#       })), nrow = 1)
#
#       colnames(estimates) <- paste0(c(paste0("q", alpha.level / 2, "."),
#                                       "q0.5.",
#                                       paste0("q", 1 - alpha.level / 2, ".")),
#                                     rep(substr(cohort.names, 3, nchar(cohort.names)), each = 3))
#
#       results.per.method.list[[method.name]] <- estimates
#
#     }
#
#     results.list[[s]] <- results.per.method.list
#
#   }
#
#   results.list <- ListPerMethod(results.list)
#
#   return (results.list)
#
# }
#
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 06 Feb 2020
# ##  These functions have been updated to remove the folder structure
# createAnalysisPath <- function (
#
#   save_path,
#   scenario_number,
#   analysis_number = 0
#
# ) {
#
#   scenario_path <- file.path(save_path, paste0("Scenario ", scenario_number))
#   if (identical(analysis_number, 0)) {
#
#     analysis_number <- sum(grepl("Analysis", list.files(scenario_path))) + 1L
#
#   }
#
#   analysis_path <- file.path(scenario_path, paste0("Analysis ", analysis_number))
#   if (!dir.exists(analysis_path)) {
#
#     dir.create(analysis_path, recursive = TRUE)
#
#   }
#
#   return (analysis_path)
#
# }
#
# saveAnalyses <- function (
#
#   scenario_analyses_list,
#   save_path        = tempdir(),
#   analysis_numbers = NULL
#
# ) {
#
#   scenario_numbers <- sapply(scenario_analyses_list, function (x) x$scenario_data$scenario_number)
#
#   if (is.null(analysis_numbers)) {
#     analysis_numbers <- rep(0, length(scenario_analyses_list))
#   }
#
#   for (s in seq_along(scenario_analyses_list)) {
#
#     ## Create directory to save the results
#     analysis_path <- createAnalysisPath(
#       save_path       = save_path,
#       scenario_number = scenario_numbers[s],
#       analysis_number = analysis_numbers[s])
#
#     if (identical(analysis_numbers[s], 0)) {
#       analysis_numbers[s] <- as.numeric(substr(analysis_path,
#                                                nchar(analysis_path),
#                                                nchar(analysis_path)))
#     }
#
#     file_name <- paste0("analysis_data_", scenario_numbers[s], "_", analysis_numbers[s],".rds")
#     saveRDS(scenario_analyses_list[[s]], file = file.path(analysis_path, file_name))
#
#   }
#
#   return (list(scenario_numbers = scenario_numbers, analysis_numbers = analysis_numbers))
#
# }
# loadAnalyses <- function (
#
#   scenario_numbers,
#   analysis_numbers,
#   load_path = tempdir()
#
# ) {
#
#   if (!identical(length(scenario_numbers), length(analysis_numbers))) {
#     stop ("scenario_numbers and analysis_numbers must have equal length_")
#   }
#
#   scenario_analyses_list <- vector(mode = "list", length = length(scenario_numbers))
#
#   for (s in seq_along(scenario_numbers)) {
#
#     scenario_analyses_list[[s]] <- readRDS(
#       file.path(load_path,
#                 paste0("Scenario ", scenario_numbers[s]),
#                 paste0("Analysis ", analysis_numbers[s]),
#                 paste0("analysis_data_", scenario_numbers[s], "_", analysis_numbers[s], ".rds"))
#     )
#
#   }
#
#   names(scenario_analyses_list) <- paste0("scenario_", scenario_numbers)
#
#   return (scenario_analyses_list)
#
# }
# saveScenarios <- function (
#
#   scenario_data_list,
#   save_path = tempdir()
#
# ) {
#
#   if (!dir.exists(save_path)) {
#     dir.create(save_path)
#   }
#
#   scenario_numbers <- numeric(length(scenario_data_list))
#
#   for (s in seq_along(scenario_data_list)) {
#
#     scenario_data       <- scenario_data_list[[s]]
#     scenario_numbers[s] <- scenario_data$scenario_number
#
#     path_create <- paste0(save_path, "/Scenario ", scenario_numbers[s])
#     if (!dir.exists(path_create)) {
#       dir.create(path_create) #unlink(path_create, recursive = TRUE)
#     }
#
#     saveRDS(scenario_data, file = paste0(path_create,
#                                          "/scenario_data_",
#                                          scenario_numbers[s], ".rds"))
#
#   }
#
#   return (scenario_numbers)
#
# }
#
# loadScenarios <- function (
#
#   scenario_numbers,
#   load_path = tempdir()
#
# ) {
#
#   scenario_data_list <- vector(mode = "list", length = length(scenario_numbers))
#
#   for (s in seq_along(scenario_numbers)) {
#
#     scenario_data_list[[s]] <- readRDS(paste0(load_path, "/Scenario ", scenario_numbers[s],
#                                               "/scenario_data_", scenario_numbers[s], ".rds"))
#
#   }
#
#   return (scenario_data_list)
#
# }
#
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 06 Feb 2020
# ##  This functions have been combined
# getBiasMSE <- function (
#
#   scenario_analyses_list,
#   cohort_names    = NULL,
#   point_estimator = "median"
#
# ) {
#
#   if (is.null(cohort_names)) {
#
#     cohort_names <- getAllCohortNames(scenario_analyses_list)
#
#   }
#
#   if (point_estimator == "median") {
#     gamma_levels <- rep(0.5, times = length(cohort_names))
#   } else if (point_estimator == "mean") {
#     gamma_levels <- rep("mean", times = length(cohort_names))
#   } else {
#     stop ("point_estimator must be either 'median' or 'mean'")
#   }
#
#   results_list <- vector(mode = "list", length = length(scenario_analyses_list))
#   names(results_list) <- names(scenario_analyses_list)
#
#   for (s in seq_along(scenario_analyses_list)) {
#
#     true_rr       <- scenario_analyses_list[[s]]$scenario_data$response_rates
#     true_rr       <- true_rr[true_rr > 0 & true_rr < 1]
#     analysis_data <- scenario_analyses_list[[s]]
#     method_names  <- analysis_data$analysis_parameters$method_names
#
#     results_per_method_list <- vector(mode = "list", length = length(method_names))
#     names(results_per_method_list) <- method_names
#
#     for (method_name in method_names) {
#
#       gamma_quantiles <- getPosteriorGammaQuantiles(
#         method_name              = method_name,
#         gamma_levels             = gamma_levels,
#         quantiles                = analysis_data$analysis_parameters$quantiles,
#         posterior_quantiles_list = analysis_data$quantiles_list,
#         cohort_names             = cohort_names)
#
#       matrix_estimates <- do.call(rbind, gamma_quantiles)
#
#       point_estimates  <- as.matrix(colMeans(matrix_estimates))
#       var_estimates    <- as.matrix(apply(matrix_estimates, 2, function (x) {
#         stats::var(x)
#       }))
#
#       bias_estimates   <- point_estimates - true_rr
#       mse_estimates    <- bias_estimates^2 + var_estimates
#
#       estimates <- cbind(bias_estimates, mse_estimates)
#
#       colnames(estimates) <- c("bias", "mse")
#
#       results_per_method_list[[method_name]] <- estimates
#
#     }
#
#     results_list[[s]] <- results_per_method_list
#
#   }
#
#   if (length(results_list) == 1) {
#
#     results_list <- results_list[[1]]
#
#   } else {
#
#     results_list <- listPerMethod(results_list)
#
#   }
#
#   return (results_list)
#
# }
# getRespRatesEstimates <- function (
#
#   scenario_analyses_list,
#   cohort_names = NULL,
#   alpha_level  = 0.05
#
# ) {
#
#   if (is.null(cohort_names)) {
#
#     cohort_names <- getAllCohortNames(scenario_analyses_list)
#
#   }
#
#   gamma_levels <- matrix(rep(c(1 - alpha_level / 2, 0.5, alpha_level / 2),
#                              each = length(cohort_names)), nrow = 3, byrow = TRUE)
#
#   results_list <- vector(mode = "list", length = length(scenario_analyses_list))
#   names(results_list) <- names(scenario_analyses_list)
#
#   for (s in seq_along(scenario_analyses_list)) {
#
#     analysis_data <- scenario_analyses_list[[s]]
#     method_names  <- analysis_data$analysis_parameters$method_names
#
#     results_per_method_list <- vector(mode = "list", length = length(method_names))
#     names(results_per_method_list) <- method_names
#
#     for (method_name in method_names) {
#
#       estimates <- apply(gamma_levels, 1, function (glevel) {
#
#         colMeans(do.call(rbind, getPosteriorGammaQuantiles(
#           method_name              = method_name,
#           gamma_levels             = glevel,
#           quantiles                = analysis_data$analysis_parameters$quantiles,
#           posterior_quantiles_list = analysis_data$quantiles_list,
#           cohort_names             = cohort_names)))
#
#       })
#
#       colnames(estimates) <- paste0(c(alpha_level / 2, 0.5, 1 - alpha_level / 2) * 100, "%")
#
#       results_per_method_list[[method_name]] <- estimates
#
#     }
#
#     results_list[[s]] <- results_per_method_list
#
#   }
#
#   results_list <- listPerMethod(results_list)
#
#   return (results_list)
#
# }
#
#
# ## ------------------------------------------------------------------------------------------------
#
# ## RIP 13 Feb 2020
# ##  This function should now aim at the overall go probability
# getGoBoundaries <- function (
#
#   scenario_analysis_list,
#   cohort_names,
#   go_rates,
#   gamma_levels_list,
#   method_name
#
# ) {
#
#   boundaries <- sapply(seq_along(cohort_names), function (n) {
#
#     stats::uniroot(
#
#       f = function (x) {
#
#         go_probs_list <- getGoProbabilities(
#           getGoDecisions(
#             scenario_analyses_list = scenario_analysis_list,
#             cohort_names           = cohort_names[n],
#             boundary_rules_list    = bquote(x[1] > .(x)),
#             gamma_levels_list      = bquote(list(.(gamma_levels_list[[n]])))))
#
#         return (go_probs_list[[method_name]][[1]][1, 1] - go_rates[n])
#
#       },
#
#       interval = c(0, 1)
#
#     )$root
#
#   })
#
#   return (boundaries)
#
# }
