#'
#' #' @title PlotCumMovAvgGoProbs
#' #' @description This function plots the cummulative moving average of the go probabilities per
#' #' scenario for a given method.
#' #' @param scenario.decisions.list A list with decisions per scenario created with GetGoDecisions()
#' #' @param method.name A strings for the name of the method to be used.
#' #' @rdname PlotCumMovAvgGoProbs
#' #' @export
#' PlotCumMovAvgGoProbs <- function (
#'
#'   scenario.decisions.list,
#'   method.name
#'
#' ) {
#'
#'   n.scenarios <- length(scenario.decisions.list)
#'
#'   par(mfrow = c(n.scenarios, 1),
#'       mar   = c(2, 4, 0, 0) + 0.1)
#'
#'   for (s in seq_along(scenario.decisions.list)) {
#'
#'     if (!method.name %in%
#'         scenario.decisions.list[[s]]$analysis.data$analysis.parameters$method.names) {
#'
#'       stop ("method.name not in scenario.decisions.list")
#'
#'     }
#'
#'     decisions <- scenario.decisions.list[[s]]$decisions.list[[method.name]]
#'
#'     cumMA <- apply(decisions, 2, CummulativeMovingAverage)
#'
#'     plot(0, type = "n",
#'          xlab = "",
#'          ylab = "Go Probability",
#'          xlim = c(1, nrow(cumMA)),
#'          ylim = c(0, 1))
#'
#'     grid()
#'
#'     for (i in seq_len(ncol(cumMA))) {
#'       lines(cumMA[, i], col = i)
#'     }
#'
#'     legend("topright", bty = "n", lty = 1, cex = 1.25,
#'            col = rgb(0, 0, 0, max = 1, alpha = 0),
#'            legend = paste0(
#'              "Scenario ",
#'              scenario.decisions.list[[s]]$scenario.data$scenario.number))
#'
#'     if (s == 1) {
#'
#'       legend("topleft", legend = FirstUpper(colnames(decisions)),
#'              lty = 1, col = seq_len(ncol(cumMA)))
#'
#'     }
#'
#'   }
#'
#' }
