


#' @title getPriorParameters
#' @md
#' @description This function provides default prior parameters for the analysis methods
#' that can be used in \code{\link[bhmbasket]{performAnalyses}}.
#' @param method_names A vector of strings for the names of the methods to be used.
#' Available methods: `c("berry", "exnex", "exnex_adj", "pooled", "stratified")`
#' @param target_rates A vector of numerics in `(0, 1)` for the
#' target rate of each cohort
#' @param n_worth An integer for the number of subjects the variability of the prior should reflect
#' response rate scale, Default: `1`
#' @param tau_scale A numeric for the scale parameter of the Half-normal distribution of \eqn{\tau}
#' in the methods `"berry"`, `"exnex"`, and `"exnex_adj"`, Default: `1`
#' @param w_j A numeric in `(0, 1)` for the weight of the Ex component in the methods `"exnex"`
#' and `"exnex_adj"`, Default: `0.5`
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' Regarding the default prior parameters for `"berry"`, `"exnex"`, and `"exnex_adj"`:
#' \itemize{
#'   \item `"berry"`: The mean of \eqn{\mu} is set to `0`.
#'   Its variance is calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`.
#'   The scale parameter of \eqn{\tau} is set to `tau_scale`.
#'   \item `"exnex"`: The weight of the Ex component is set to `w_j`.
#'   For the Ex component:
#'   The target rate that results in the greatest variance is determined.
#'   The mean of \eqn{\mu} is set to that target rate.
#'   The variance of \eqn{\mu} is calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`.
#'   The scale parameter of \eqn{\tau} is set to `tau_scale`.
#'   For the Nex components:
#'   The means of \eqn{\mu_j} are set to the respective target rates.
#'   The variances of \eqn{\tau_j} are calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`, see also \code{\link[bhmbasket]{getMuVar}}.
#'   \item `"exnex_adj"`: The weight of the Ex component is set to `w_j`.
#'   For the Ex component:
#'   The target rate that results in the greatest variance is determined.
#'   The mean of \eqn{\mu} is set to `0`.
#'   The variance of \eqn{\mu} is calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`, see also \code{\link[bhmbasket]{getMuVar}}.
#'   The scale parameter of \eqn{\tau} is set to `tau_scale`.
#'   For the Nex components:
#'   The means of \eqn{\mu_j} are set to the `0`.
#'   The variances of \eqn{\tau_j} are calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`, see also \code{\link[bhmbasket]{getMuVar}}.
#'   \item `"pooled"`: The target rate that results in the greatest variance is determined.
#'   The scale parameter \eqn{\alpha} is set to that target rate times `n_worth`.
#'   The scale parameter \eqn{\beta} is set to 1 - that target rate times `n_worth`.
#'   \item `"stratified"`:
#'   The scale parameters \eqn{\alpha_j} are set to `target_rates * n_worth`.
#'   The scale parameters \eqn{\beta_j} are set to `(1 - target_rates) * n_worth`.
#' }
#' @author Stephan Wojciekowski
#' @examples
#' prior_parameters_list <- getPriorParameters(
#'   method_names = c("berry", "exnex", "exnex_adj", "pooled", "stratified"),
#'   target_rates = c(0.1, 0.2, 0.3))
#' @rdname getPriorParameters
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @references Berry, Scott M., et al. "Bayesian hierarchical modeling of patient subpopulations:
#' efficient designs of phase II oncology clinical trials."
#' \emph{Clinical Trials} 10.5 (2013): 720-734.
#' @references Neuenschwander, Beat, et al. "Robust exchangeability designs
#' for early phase clinical trials with multiple strata."
#' \emph{Pharmaceutical statistics} 15.2 (2016): 123-134.
#' @export
getPriorParameters <- function (

  method_names,
  target_rates,
  n_worth   = 1,
  tau_scale = 1,
  w_j       = 0.5

) {

  error_method_names <- simpleError(
    paste("Please provide a (vector of) strings for the argument 'method_names'\n",
          "Must be one of 'berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'"))
  error_target_rates <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'target_rates'")
  error_tau_scale    <- simpleError(
    "Please provide a positive numeric for the argument 'tau_scale'")
  error_n_worth      <- simpleError("Please provide a positive integer for the argument 'n_worth'")
  error_w_j          <- simpleError("Please provide a numeric in (0, 1) for the argument 'w_j'")

  if (missing(method_names)) stop (error_method_names)
  if (missing(target_rates)) stop (error_target_rates)

  method_names <- tryCatch({

    match.arg(
      method_names,
      choices    = c('berry', 'exnex', 'exnex_adj', 'pooled', 'stratified'),
      several.ok = TRUE)

  }, error = function (e) e)

  if (inherits(method_names, "error"))          stop (error_method_names)
  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.numeric(tau_scale))   stop (error_tau_scale)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)
  # if (!is.single.numeric.in.zero.one(w_j))      stop (error_w_j)
  if (!is.single.numeric(w_j) ||
      any(w_j < 0) ||
      any(w_j > 1))                             stop (error_w_j)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  method_names <- sort(method_names)

  prior_parameters_list <- lapply(method_names, function (method_name) {

    if (method_name == "berry") {

      getPriorParametersBerry(
        target_rates = target_rates,
        n_worth      = n_worth,
        tau_scale    = tau_scale)[[1]]

    } else if (method_name == "exnex") {

      getPriorParametersExNex(
        target_rates = target_rates,
        n_worth      = n_worth,
        tau_scale    = tau_scale,
        w_j          = w_j)[[1]]

    } else if (method_name == "exnex_adj") {

      getPriorParametersExNexAdj(
        target_rates = target_rates,
        n_worth      = n_worth,
        tau_scale    = tau_scale,
        w_j          = w_j)[[1]]

    } else if (method_name == "pooled") {

      getPriorParametersPooled(
        target_rates = target_rates,
        n_worth      = n_worth)[[1]]

    } else if (method_name == "stratified") {

      getPriorParametersStratified(
        target_rates = target_rates,
        n_worth      = n_worth)[[1]]

    }

  })

  names(prior_parameters_list) <- method_names
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title combinePriorParameters
#' @md
#' @description This function combines prior parameters from different sources and returns them
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param list_of_prior_parameters A list of items with class `prior_parameters_list`
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' This function is intended to combine the prior parameters set with the functions
#' \code{\link[bhmbasket]{setPriorParametersBerry}},
#' \code{\link[bhmbasket]{setPriorParametersExNex}},
#' \code{\link[bhmbasket]{setPriorParametersExNexAdj}},
#' \code{\link[bhmbasket]{setPriorParametersPooled}}, or
#' \code{\link[bhmbasket]{setPriorParametersStratified}},
#' in case more than one analysis method should be applied with
#' \code{\link[bhmbasket]{performAnalyses}}.
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_stratified <- setPriorParametersStratified(c(1, 2), c(3, 4))
#'  prior_parameters_berry      <- setPriorParametersBerry(0, 1, 2)
#'
#'  prior_parameters_list       <- combinePriorParameters(
#'    list(prior_parameters_berry,
#'         prior_parameters_stratified))
#' @rdname combinePriorParameters
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @export
combinePriorParameters <- function (

  list_of_prior_parameters

) {

  error_list <- simpleError("Please provide a list of of items with class 'prior_parameters_list'")

  if (missing(list_of_prior_parameters)) stop (error_list)

  if (!is.list(list_of_prior_parameters))                               stop (error_list)
  if (any(!sapply(list_of_prior_parameters, is.prior_parameters_list))) stop (error_list)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  method_names <- sapply(list_of_prior_parameters, names)

  if (!identical(length(unique(method_names)), length(method_names))) stop (simpleError(
    "Please provide only one 'prior_parameters_list' per analysis method"
  ))

  prior_parameters_list <- vector(mode = "list", length(method_names))
  names(prior_parameters_list) <- method_names

  for (n in seq_along(method_names)) {

    prior_parameters_list[[n]] <- list_of_prior_parameters[[n]][[1]]

  }

  prior_parameters_list <- prior_parameters_list[sort(names(prior_parameters_list))]
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}


is.prior_parameters_list <- function(x) {

  if (missing(x)) stop ("Please provide an object for the argument 'x'")

  inherits(x, "prior_parameters_list")

}

#' @title getMuVar
#' @md
#' @description This function returns the variance of \eqn{\mu} that is worth a certain number
#' of subjects for the distribution of the response rates.
#' @param response_rate A numeric for the response rate
#' @param tau_scale A numeric for the scale parameter of the Half-normal distribution of \eqn{\tau}
#' @param n_worth An integer for the number of subjects the variance of \eqn{\mu} should be worth
#' with regard to the variability of the distribution of the response rate, Default: `1`
#' @return Returns a numeric for the variance of \eqn{\mu}
#' @details Calculates the variance `mu_var` in
#' \deqn{logit(p) = \theta ~ N(\mu, \tau),
#' \mu ~ N(mu_mean, mu_var), \tau ~ HN(tau_scale),}
#' for `n_worth` number of observations, as in Neuenschwander et al. (2016).
#' @rdname getMuVar
#' @author Stephan Wojciekowski
#' @examples
#'   getMuVar(response_rate = 0.3,
#'            tau_scale     = 1)
#'   getMuVar(response_rate = 0.3,
#'            tau_scale     = 1,
#'            n_worth       = 2)
#' @references Neuenschwander, Beat, et al. "Robust exchangeability designs
#' for early phase clinical trials with multiple strata."
#' \emph{Pharmaceutical statistics} 15.2 (2016): 123-134.
#' @export
getMuVar <- function (

  response_rate,
  tau_scale,
  n_worth = 1

) {

  error_response_rate <- simpleError(
    "Please provide a numeric in (0, 1) for the argument 'response_rate'")
  error_tau_scale <- simpleError(
    "Please provide a positive numeric for the argument 'tau_scale'")
  error_n_worth <- simpleError(
    "Please provide a positive integer for the argument 'n_worth'")

  if (missing(response_rate)) stop (error_response_rate)
  if (missing(tau_scale))     stop (error_tau_scale)

  if (!is.numeric.in.zero.one(response_rate))   stop (error_response_rate)
  if (!is.numeric(tau_scale) ||
      length(tau_scale) > 1 ||
      tau_scale < 0)                            stop (error_tau_scale)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  mu_var <- (n_worth * response_rate * (1 - response_rate))^-1 - tau_scale^2

  return (mu_var)

}

## berry ####

getPriorParametersBerry <- function (

  target_rates,
  tau_scale = 1,
  n_worth   = 1

) {

  error_target_rates <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'target_rates'")
  error_tau_scale    <- simpleError(
    "Please provide a positive numeric for the argument 'tau_scale'")
  error_n_worth      <- simpleError(
    "Please provide a positive integer for the argument 'n_worth'")

  if (missing(target_rates))                    stop (error_target_rates)

  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.numeric(tau_scale))   stop (error_tau_scale)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]

  mu_var <- getMuVar(target_rate_max_var, tau_scale, n_worth)

  if (mu_var <= 0) stop(simpleError(paste(
    "The provided input parameters lead to a variance of mu <= 0.",
    "Consider to decrease 'tau_scale' or 'n_worth'")))

  prior_parameters <- list(
    mu_mean   = 0,
    mu_sd     = mu_var^0.5,
    tau_scale = tau_scale)

  prior_parameters_list <- list (berry = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersBerry
#' @md
#' @description This function sets prior parameters for the analysis method `"berry"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param mu_mean A numeric for the mean of \eqn{\mu}
#' @param mu_sd A positive numeric for the standard deviation of \eqn{\mu}
#' @param tau_scale A positive numeric for the scale parameter of \eqn{\tau}
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' This function sets the prior parameters for the method proposed by Berry et al. (2013).
#' Note that the implemented distribution of \eqn{\tau} is half-normal.
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_berry <- setPriorParametersBerry(0, 1, 2)
#' @rdname setPriorParametersBerry
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @references Berry, Scott M., et al. "Bayesian hierarchical modeling of patient subpopulations:
#' efficient designs of phase II oncology clinical trials."
#' \emph{Clinical Trials} 10.5 (2013): 720-734.
#' @export
setPriorParametersBerry <- function (

  mu_mean,
  mu_sd,
  tau_scale

) {

  error_mu_mean  <- simpleError(
    "Please provide a numeric for the argument 'mu_mean'")
  error_mu_sd  <- simpleError(
    "Please provide a positive numeric for the argument 'mu_sd'")
  error_tau_scale <- simpleError(
    "Please provide a positive numeric for the argument 'tau_scale'")

  if (missing(mu_mean))   stop (error_mu_mean)
  if (missing(mu_sd))     stop (error_mu_sd)
  if (missing(tau_scale)) stop (error_tau_scale)

  if (!is.single.numeric(mu_mean))            stop (error_mu_mean)
  if (!is.single.positive.numeric(mu_sd))     stop (error_mu_sd)
  if (!is.single.positive.numeric(tau_scale)) stop (error_tau_scale)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- list(
    mu_mean   = mu_mean,
    mu_sd     = mu_sd,
    tau_scale = tau_scale)

  prior_parameters_list <- list (berry = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

## exnex ####

getPriorParametersExNex <- function (

  target_rates,
  tau_scale = 1,
  n_worth   = 1,

  w_j       = 0.5

) {

  error_target_rates <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'target_rates'")
  error_tau_scale    <- simpleError(
    "Please provide a positive numeric for the argument 'tau_scale'")
  error_n_worth      <- simpleError(
    "Please provide a positive integer for the argument 'n_worth'")
  error_w_j          <- simpleError(
    "Please provide a numeric in (0, 1) for the argument 'w_j'")

  if (missing(target_rates)) stop (error_target_rates)

  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.numeric(tau_scale))   stop (error_tau_scale)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)
  if (!is.single.numeric(w_j) ||
      any(w_j < 0) ||
      any(w_j > 1))                             stop (error_w_j)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]

  mu_var <- getMuVar(target_rate_max_var, tau_scale, n_worth)

  if (mu_var <= 0) stop(simpleError(paste(
    "The provided input parameters lead to a variance of mu <= 0.",
    "Consider to decrease 'tau_scale' or 'n_worth'")))

  prior_parameters <- list(
    mu_mean   = logit(target_rate_max_var),
    mu_sd     = mu_var^0.5,
    tau_scale = tau_scale,

    mu_j  = logit(target_rates),
    tau_j = getMuVar(target_rates, 0, n_worth)^0.5,

    w_j   = w_j)

  prior_parameters_list <- list (exnex = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersExNex
#' @md
#' @description This function sets prior parameters for the analysis method `"exnex"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param mu_mean A numeric for the mean of \eqn{\mu}
#' @param mu_sd A positive numeric for the standard deviation of \eqn{\mu}
#' @param tau_scale A positive numeric for the scale parameter of \eqn{\tau}
#' @param mu_j A vector of numerics for the means \eqn{\mu_j}
#' @param tau_j A vector of positive numerics for the standard deviations \eqn{\tau_j}
#' @param w_j A numeric in `(0, 1)` for the weight of the Ex component
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' This function sets the prior parameters for the method proposed by Neuenschwander et al. (2016).
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_exnex <- setPriorParametersExNex(0, 1, 2, c(4, 5), c(6, 7), 0.8)
#' @rdname setPriorParametersExNex
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @references Neuenschwander, Beat, et al. "Robust exchangeability designs
#' for early phase clinical trials with multiple strata."
#' \emph{Pharmaceutical statistics} 15.2 (2016): 123-134.
#' @export
setPriorParametersExNex <- function (

  mu_mean,
  mu_sd,
  tau_scale,

  mu_j,
  tau_j,

  w_j

) {

  error_mu_mean   <- simpleError("Please provide a numeric for the argument 'mu_mean'")
  error_mu_sd     <- simpleError("Please provide a positive numeric for the argument 'mu_sd'")
  error_tau_scale <- simpleError("Please provide a positive numeric for the argument 'tau_scale'")
  error_mu_j      <- simpleError("Please provide a (vector of) numeric(s) for the argument 'mu_j'")
  error_tau_j     <- simpleError(
    "Please provide a (vector of) positive numeric(s) for the argument 'tau_j'")
  error_w_j   <- simpleError("Please provide a numeric in (0, 1) for the argument 'w_j'")

  if (missing(mu_mean))   stop (error_mu_mean)
  if (missing(mu_sd))     stop (error_mu_sd)
  if (missing(tau_scale)) stop (error_tau_scale)
  if (missing(mu_j))      stop (error_mu_j)
  if (missing(tau_j))     stop (error_tau_j)
  if (missing(w_j))       stop (error_w_j)

  if (!is.numeric(mu_mean))                   stop (error_mu_mean)
  if (!is.positive.numeric(mu_sd))            stop (error_mu_sd)
  if (!is.single.positive.numeric(tau_scale)) stop (error_tau_scale)
  if (!is.numeric(mu_j))                      stop (error_mu_j)
  if (!is.positive.numeric(tau_j))            stop (error_tau_j)
  if (!is.numeric(w_j) ||
      !all(w_j >= 0) ||
      !all(w_j <= 1))                         stop (error_w_j)

  if (!identical(length(mu_mean), length(mu_sd))) stop (simpleError(
    "'mu_mean' and 'mu_sd' must have save length"))

  if (!identical(length(mu_mean), 1L)) {
    if (!identical(length(w_j), length(mu_mean) + 1L)) stop (simpleError(
      "'w_j' must have length equal to length(mu_mean) + 1 if length(mu_mean) > 1"))
  }

  if (!identical(length(w_j), 1L)) {
    if (!identical(length(w_j), length(mu_mean) + 1L)) stop (simpleError(
      "'w_j' must have length 1 or 2, if length(mu_mean) = 1"))
    if (!isTRUE(all.equal(sum(w_j), 1))) stop (simpleError(
      "Sum over items in 'w_j' must equal 1 if length(w_j) > 1"))
  }

  if (!identical(length(mu_j), length(tau_j)))
    stop (simpleError("mu_j and tau_j must have the same length"))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <-  list(mu_mean   = mu_mean,
                            mu_sd     = mu_sd,
                            tau_scale = tau_scale,
                            mu_j      = mu_j,
                            tau_j     = tau_j,
                            w_j       = w_j)

  prior_parameters_list <- list (exnex = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

## exnex_adj ####

getPriorParametersExNexAdj <- function (

  target_rates,
  tau_scale = 1,
  n_worth   = 1,

  w_j       = 0.5

) {

  error_target_rates <- simpleError(
    paste("Please provide a vector of numerics in (0, 1) for the argument 'target_rates'"))
  error_tau_scale    <- simpleError("Please provide a positive numeric for the argument 'tau_scale'")
  error_n_worth      <- simpleError("Please provide a positive integer for the argument 'n_worth'")
  error_w_j          <- simpleError("Please provide a numeric in (0, 1) for the argument 'w_j'")

  if (missing(target_rates)) stop (error_target_rates)

  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.numeric(tau_scale))   stop (error_tau_scale)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)
  if (!is.single.numeric(w_j) ||
      any(w_j < 0) ||
      any(w_j > 1))                             stop (error_w_j)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- getPriorParametersExNex(target_rates, tau_scale, n_worth, w_j)[[1]]

  prior_parameters$mu_mean <- 0
  prior_parameters$mu_j    <- rep(0, length(target_rates))

  prior_parameters_list <- list (exnex_adj = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersExNexAdj
#' @md
#' @description This function sets prior parameters for the analysis method `"exnex_adj"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param mu_mean A numeric for the mean of \eqn{\mu}
#' @param mu_sd A positive numeric for the standard deviation of \eqn{\mu}
#' @param tau_scale A positive numeric for the scale parameter of \eqn{\tau}
#' @param mu_j A vector of numerics for the means \eqn{\mu_j}
#' @param tau_j A vector of positive numerics for the standard deviations \eqn{\tau_j}
#' @param w_j A numeric in `(0, 1)` for the weight of the Ex component
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' This function sets the prior parameters for the method ExNex Adjusted, which combines
#' the approach proposed by Neuenschwander et al. (2016) and the approach proposed by
#' Berry et al. (2013).
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_exnex_adj <- setPriorParametersExNexAdj(0, 1, 2, c(4, 5), c(6, 7), 0.8)
#' @rdname setPriorParametersExNexAdj
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @export
setPriorParametersExNexAdj <- function (

  mu_mean,
  mu_sd,
  tau_scale,

  mu_j,
  tau_j,

  w_j

) {

  error_mu_mean   <- simpleError("Please provide a numeric for the argument 'mu_mean'")
  error_mu_sd     <- simpleError("Please provide a positive numeric for the argument 'mu_sd'")
  error_tau_scale <- simpleError("Please provide a positive numeric for the argument 'tau_scale'")
  error_mu_j      <- simpleError("Please provide a (vector of) numeric(s) for the argument 'mu_j'")
  error_tau_j     <- simpleError(
    "Please provide a (vector of) positive numeric(s) for the argument 'tau_j'")
  error_w_j   <- simpleError("Please provide a numeric in (0, 1) for the argument 'w_j'")

  if (missing(mu_mean))   stop (error_mu_mean)
  if (missing(mu_sd))     stop (error_mu_sd)
  if (missing(tau_scale)) stop (error_tau_scale)
  if (missing(mu_j))      stop (error_mu_j)
  if (missing(tau_j))     stop (error_tau_j)
  if (missing(w_j))       stop (error_w_j)

  if (!is.numeric(mu_mean))                    stop (error_mu_mean)
  if (!is.positive.numeric(mu_sd))             stop (error_mu_sd)
  if (!is.single.positive.numeric(tau_scale))  stop (error_tau_scale)
  if (!is.numeric(mu_j))                       stop (error_mu_j)
  if (!is.positive.numeric(tau_j))             stop (error_tau_j)
  if (!is.single.numeric(w_j) ||
      any(w_j < 0) ||
      any(w_j > 1))                             stop (error_w_j)

  if (!identical(length(mu_mean), length(mu_sd))) stop (simpleError(
    "'mu_mean' and 'mu_sd' must have save length"))

  if (!identical(length(mu_mean), 1L)) {
    if (!identical(length(w_j), length(mu_mean) + 1L)) stop (simpleError(
      "'w_j' must have length equal to length(mu_mean) + 1 if length(mu_mean) > 1"))
  }

  if (!identical(length(w_j), 1L)) {
    if (!identical(length(w_j), length(mu_mean) + 1L)) stop (simpleError(
      "'w_j' must have length 1 or 2, if length(mu_mean) = 1"))
  }

  if (!identical(length(mu_j), length(tau_j)))
    stop (simpleError("mu_j and tau_j must have the same length"))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters_list <- setPriorParametersExNex(mu_mean, mu_sd, tau_scale, mu_j, tau_j, w_j)
  names(prior_parameters_list) <- "exnex_adj"

  return (prior_parameters_list)

}

## pooled ####

getPriorParametersPooled <- function (

  target_rates,
  n_worth = 1

) {

  error_target_rates <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'target_rates'")
  error_n_worth      <- simpleError(
    "Please provide a positive integer for the argument 'n_worth'")

  if (missing(target_rates)) stop (error_target_rates)

  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  target_rate_min_var <- target_rates[abs(target_rates - 0.5) == min(abs(target_rates - 0.5))][1]

  a <- target_rate_min_var * n_worth
  b <- (1 - target_rate_min_var) * n_worth

  prior_parameters <- list(a = a, b = b)

  prior_parameters_list        <- list(pooled = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersPooled
#' @md
#' @description This function sets prior parameters for the analysis method `"pooled"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param a A positive numeric for \eqn{\alpha}
#' @param b A positive numeric for \eqn{\beta}
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' The method `"pooled"` is a beta-binomial model that pools all cohorts.
#' The prior parameters are the scale parameters of the beta prior distribution.
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_pooled <- setPriorParametersPooled(1, 2)
#' @rdname setPriorParametersPooled
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @export
setPriorParametersPooled <- function (

  a,
  b

) {

  error_a <- simpleError(
    "Please provide a positive numeric for the argument 'a'")
  error_b <- simpleError(
    "Please provide a positive numeric for the argument 'b'")

  if (missing(a)) stop (error_a)
  if (missing(b)) stop (error_b)

  if (!is.single.positive.numeric(a)) stop (error_a)
  if (!is.single.positive.numeric(b)) stop (error_b)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- list(a = a, b = b)

  prior_parameters_list        <- list(pooled = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

## stratified ####

getPriorParametersStratified <- function (

  target_rates,
  n_worth = 1

) {

  error_target_rates <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'target_rates'")
  error_n_worth      <- simpleError(
    "Please provide a positive integer for the argument 'n_worth'")

  if (missing(target_rates)) stop (error_target_rates)

  if (!is.numeric.in.zero.one(target_rates))    stop (error_target_rates)
  if (!is.single.positive.wholenumber(n_worth)) stop (error_n_worth)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  a_j <- target_rates * n_worth
  b_j <- (1 - target_rates) * n_worth

  prior_parameters <- list(a_j = a_j, b_j = b_j)

  prior_parameters_list        <- list(stratified = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersStratified
#' @md
#' @description This function sets prior parameters for the analysis method `"stratified"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param a_j A vector of positive numerics for \eqn{\alpha}
#' @param b_j A vector of positive numerics for \eqn{\beta}
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' The method `"stratified"` is a beta-binomial model that assesses each cohort individually.
#' The prior parameters are the scale parameters of the beta prior distributions.
#' @author Stephan Wojciekowski
#' @examples
#'  prior_parameters_pooled <- setPriorParametersStratified(c(1, 2), c(3, 4))
#' @rdname setPriorParametersStratified
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{getMuVar}}
#' @export
setPriorParametersStratified <- function (

  a_j,
  b_j

) {

  error_a_j <- simpleError(
    "Please provide a (vector of) positive numeric(s) in for the argument 'a_j'")
  error_b_j <- simpleError(
    "Please provide a (vector of) positive numeric(s) in for the argument 'b_j'")

  if (missing(a_j)) stop (error_a_j)
  if (missing(b_j)) stop (error_b_j)

  if (!is.positive.numeric(a_j)) stop (error_a_j)
  if (!is.positive.numeric(b_j)) stop (error_b_j)

  if (!identical(length(a_j), length(b_j)))
    stop (simpleError("a_j and b_j must have the same length"))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- list(a_j = a_j, b_j = b_j)

  prior_parameters_list        <- list(stratified = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}


