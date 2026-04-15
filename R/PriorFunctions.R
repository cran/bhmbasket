#' @title getPriorParameters
#' @md
#' @description This function provides default prior parameters for the analysis methods
#' that can be used in \code{\link[bhmbasket]{performAnalyses}}.
#' @param method_names A vector of strings for the names of the methods to be used.
#' Available methods: `c("berry", "exnex", "exnex_adj", "exnex_mix", "exnex_adj_mix", "pooled", "stratified", "stratified_mix)`
#' @param target_rates A vector of numerics in `(0, 1)` for the
#' target rate of each cohort
#' @param n_worth An integer for the number of subjects the variability of the prior should reflect
#' response rate scale, Default: `1`
#' @param tau_scale A numeric for the scale parameter of the Half-normal distribution of \eqn{\tau}
#' in the methods `"berry"`, `"exnex"`, `"exnex"_mix`, `"exnex_adj"` and `"exnex_adj_mix"`, Default: `1`
#' @param w_j A numeric in `(0, 1)` for the weight of the Ex component in the methods `"exnex"`
#' and `"exnex_adj"`, Default: `0.5`
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' Regarding the default prior parameters for `"berry"`, `"exnex"`, `"exnex_mix"`, `"exnex_adj"` and `"exnex_adj_mix"`:
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
#'   
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
#'   
#'   For the Nex components:
#'   The means of \eqn{\mu_j} are set to the `0`.
#'   The variances of \eqn{\tau_j} are calculated as proposed in "Robust exchangeability designs for early
#'   phase clinical trials with multiple strata" (Neuenschwander et al. (2016))
#'   with regard to `n_worth`, see also \code{\link[bhmbasket]{getMuVar}}.
#'   \item `"exnex_mix"`: Uses the same default Ex prior construction as `"exnex"`.
#'   The Nex part default parameters are specified as a one-component mixture prior with
#'   `w_nex = 1`, `mean_nex = matrix(logit(target_rates), nrow = 1)`,
#'   and `sd_nex = matrix(sqrt(getMuVar(target_rates, 0, n_worth)), nrow = 1)`.
#'   This keeps the default mixture representation compatible with the `getPriorParameter()`'s input.
#'   \item `"exnex_adj_mix"`: Uses the same default Ex prior construction as `"exnex_adj"`.
#'   The Nex part is specified as a one-component mixture prior with
#'   `w_nex = 1`, `mean_nex = matrix(logit(target_rates), nrow = 1)`,
#'   and `sd_nex = matrix(sqrt(getMuVar(target_rates, 0, n_worth)), nrow = 1)`.
#'   The Ex component is centered as in `"exnex_adj"`.
#'   \item `"pooled"`: The target rate that results in the greatest variance is determined.
#'   The scale parameter \eqn{\alpha} is set to that target rate times `n_worth`.
#'   The scale parameter \eqn{\beta} is set to 1 - that target rate times `n_worth`.
#'   \item `"stratified"`:
#'   The scale parameters \eqn{\alpha_j} are set to `target_rates * n_worth`.
#'   The scale parameters \eqn{\beta_j} are set to `(1 - target_rates) * n_worth`.
#'   \item `"stratified_mix"`:
#'   A two-component beta mixture prior is created by default for each cohort.
#'   The first component uses `a_j = target_rates * n_worth` and
#'   `b_j = (1 - target_rates) * n_worth`.
#'   The second component is a vague prior with `a_j = 1` and `b_j = 1`.
#'   The default mixture weights are `c(0.8, 0.2)`.
#' }
#' @author Stephan Wojciekowski
#' @examples
#' prior_parameters_list <- getPriorParameters(
#'   method_names = c("berry", "exnex", "exnex_mix", "exnex_adj",
#'                    "exnex_adj_mix", "pooled", "stratified", "stratified_mix"),
#'   target_rates = c(0.1, 0.2, 0.3))
#' @rdname getPriorParameters
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{setPriorParametersBerry}}
#'  \code{\link[bhmbasket]{setPriorParametersExNex}}
#'  \code{\link[bhmbasket]{setPriorParametersExNexAdj}}
#'  \code{\link[bhmbasket]{setPriorParametersPooled}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#'  \code{\link[bhmbasket]{setPriorParametersStratifiedMix}}
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

  error_method_names <-
    paste("Providing a (vector of) strings for the argument 'method_names'\n",
          "Must be one of 'berry', 'exnex', 'exnex_mix', 'exnex_adj', 'exnex_adj_mix',",
          "'pooled', 'stratified', 'stratified_mix'")
  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_tau_scale    <-
    "Providing a positive numeric for the argument 'tau_scale'"
  error_n_worth      <- "Providing a positive integer for the argument 'n_worth'"
  error_w_j          <- "Providing a numeric in (0, 1) for the argument 'w_j'"

  checkmate::assertCharacter(method_names, any.missing = FALSE, .var.name = error_method_names)

  method_names <- tryCatch({

    match.arg(
      method_names,
      choices = c(
        "berry", "exnex", "exnex_adj", "pooled", "stratified", "stratified_mix",
        "exnex_mix", "exnex_adj_mix"
      ),
      several.ok = TRUE
    )

  }, error = function (e) {
    stop(error_method_names, call. = FALSE)
  })

  checkmate::assertNumeric(target_rates, any.missing = FALSE, .var.name = error_target_rates)
  checkmate::assertTRUE(all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates)

  checkmate::assertNumber(tau_scale, .var.name = error_tau_scale)
  checkmate::assertTRUE(tau_scale > 0, .var.name = error_tau_scale)

  checkmate::assertInt(n_worth, lower = 1, .var.name = error_n_worth)

  checkmate::assertNumeric(w_j, lower = 0, upper = 1, len = 1, .var.name = error_w_j)


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

    } else if (method_name == "exnex_mix") {

      getPriorParametersExNex(
        target_rates = target_rates,
        n_worth      = n_worth,
        tau_scale    = tau_scale,
        w_j          = w_j,
        w_nex        = 1,
        mean_nex     = matrix(logit(target_rates), nrow = 1),
        sd_nex       = matrix(sqrt(getMuVar(target_rates, 0, n_worth)), nrow = 1))[[1]]

    } else if (method_name == "exnex_adj_mix") {

      getPriorParametersExNexAdj(
        target_rates = target_rates,
        n_worth      = n_worth,
        tau_scale    = tau_scale,
        w_j          = w_j,
        w_nex        = 1,
        mean_nex     = matrix(logit(target_rates), nrow = 1),
        sd_nex       = matrix(sqrt(getMuVar(target_rates, 0, n_worth)), nrow = 1))[[1]]

    } else if (method_name == "pooled") {

      getPriorParametersPooled(
        target_rates = target_rates,
        n_worth      = n_worth)[[1]]

    } else if (method_name == "stratified") {

      getPriorParametersStratified(
        target_rates = target_rates,
        n_worth      = n_worth)[[1]]

    }

    else if (method_name == "stratified_mix") {

      getPriorParametersStratifiedMix(
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

  error_list <- "Providing a list of items with class 'prior_parameters_list'"

  checkmate::assertList(
    list_of_prior_parameters, types = "list", any.missing = FALSE, .var.name = error_list
    )

  checkmate::assertTRUE(
    all(vapply(list_of_prior_parameters, function(x) inherits(x, "prior_parameters_list"), logical(1))),
    .var.name = error_list
  )

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  method_names <- sapply(list_of_prior_parameters, names)

  checkmate::assertTRUE(
    length(unique(method_names)) == length(method_names),
  )

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

  error_response_rate <-
    "Providing a numeric in (0, 1) for the argument 'response_rate'"
  error_tau_scale <-
    "Providing a positive numeric for the argument 'tau_scale'"
  error_n_worth <-
    "Providing a positive integer for the argument 'n_worth'"

  checkmate::assertNumeric(response_rate, .var.name = error_response_rate)
  checkmate::assertTRUE(all(response_rate > 0 & response_rate < 1), .var.name = error_response_rate)

  checkmate::assertNumber(tau_scale, lower = 0, .var.name = error_tau_scale)

  checkmate::assertInt(n_worth, lower = 1, .var.name = error_n_worth)

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

  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_tau_scale    <-
    "Providing a positive numeric for the argument 'tau_scale'"
  error_n_worth      <-
    "Providing a positive integer for the argument 'n_worth'"

  checkmate::assertNumeric(target_rates, any.missing = FALSE, .var.name = error_target_rates)
  checkmate::assertTRUE(all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates)

  checkmate::assertNumber(tau_scale, .var.name = error_tau_scale)
  checkmate::assertTRUE(tau_scale > 0, .var.name = error_tau_scale)

  checkmate::assertInt(n_worth, lower = 1, .var.name = error_n_worth)

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

  error_mu_mean  <-
    "Providing a numeric for the argument 'mu_mean'"
  error_mu_sd  <-
    "Providing a positive numeric for the argument 'mu_sd'"
  error_tau_scale <-
    "Providing a positive numeric for the argument 'tau_scale'"


  checkmate::assertNumber(mu_mean, .var.name = error_mu_mean)

  checkmate::assertNumber(mu_sd, .var.name = error_mu_sd)
  checkmate::assertTRUE(mu_sd > 0, .var.name = error_mu_sd)

  checkmate::assertNumber(tau_scale, .var.name = error_tau_scale)
  checkmate::assertTRUE(tau_scale > 0, .var.name = error_tau_scale)

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
  w_j       = 0.5,

  w_nex     = NULL,
  mean_nex  = NULL,
  sd_nex    = NULL

) {

  error_target_rates <- "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_tau_scale    <- "Providing a positive numeric for the argument 'tau_scale'"
  error_n_worth      <- "Providing a positive integer for the argument 'n_worth'"
  error_w_j          <- "Providing a numeric in (0, 1) for the argument 'w_j'"
  error_w_nex        <- "Providing a numeric vector of weights in [0, 1] summing to 1 for the argument 'w_nex'"
  error_mean_nex     <- "Providing a numeric matrix for the argument 'mean_nex'"
  error_sd_nex       <- "Providing a positive numeric matrix for the argument 'sd_nex'"
  error_dim          <- "'mean_nex' and 'sd_nex' must have the same dimensions, nrow(mean_nex) must equal length(w_nex), and ncol(mean_nex) must equal length(target_rates)"

  checkmate::assertNumeric(target_rates, any.missing = FALSE, .var.name = error_target_rates)
  checkmate::assertTRUE(all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates)

  checkmate::assertNumber(tau_scale, .var.name = error_tau_scale)
  checkmate::assertTRUE(tau_scale > 0, .var.name = error_tau_scale)

  checkmate::assertInt(n_worth, lower = 1, .var.name = error_n_worth)

  checkmate::assertNumeric(w_j, lower = 0, upper = 1, len = 1, .var.name = error_w_j)

  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  mu_var <- getMuVar(target_rate_max_var, tau_scale, n_worth)

  if (mu_var <= 0) {
    stop(simpleError(paste(
      "The provided input parameters lead to a variance of mu <= 0.",
      "Consider to decrease 'tau_scale' or 'n_worth'"
    )))
  }

  prior_base <- list(
    mu_mean   = logit(target_rate_max_var),
    mu_sd     = sqrt(mu_var),
    tau_scale = tau_scale,
    w_j       = w_j
  )

  if (is.null(w_nex) && is.null(mean_nex) && is.null(sd_nex)) {

    prior_extra <- list(
      mu_j  = logit(target_rates),
      tau_j = sqrt(getMuVar(target_rates, 0, n_worth))
    )

  } else {

    checkmate::assertNumeric(w_nex, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w_nex)
    checkmate::assertTRUE(isTRUE(all.equal(sum(w_nex), 1)), .var.name = error_w_nex)

    checkmate::assertTRUE(is.matrix(mean_nex), .var.name = error_mean_nex)
    checkmate::assertTRUE(is.matrix(sd_nex), .var.name = error_sd_nex)
    checkmate::assertTRUE(all(sd_nex > 0), .var.name = error_sd_nex)
    checkmate::assertTRUE(all(dim(mean_nex) == dim(sd_nex)), .var.name = error_dim)
    checkmate::assertTRUE(nrow(mean_nex) == length(w_nex), .var.name = error_dim)
    checkmate::assertTRUE(ncol(mean_nex) == length(target_rates), .var.name = error_dim)

    prior_extra <- list(
      w_nex    = as.numeric(w_nex),
      mean_nex = mean_nex,
      sd_nex   = sd_nex
    )
  }

  prior_parameters_list <- list(exnex = c(prior_base, prior_extra))
  class(prior_parameters_list) <- "prior_parameters_list"

  return(prior_parameters_list)
}

#' @title setPriorParametersExNex
#' @md
#' @description This function sets prior parameters for the analysis method `"exnex"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#'
#' It supports two specifications for the Nex part:
#' \itemize{
#'   \item the standard ExNex specification with cohort-specific Nex priors via `mu_j` and `tau_j`
#'   \item an extended specification with a mixture prior on the Nex part via `w_nex`, `mean_nex`, and `sd_nex`
#' }
#'
#' @param mu_mean A numeric for the mean of \eqn{\mu}
#' @param mu_sd A positive numeric for the standard deviation of \eqn{\mu}
#' @param tau_scale A positive numeric for the scale parameter of \eqn{\tau}
#' @param mu_j A vector of numerics for the means \eqn{\mu_j} of the standard Nex priors.
#'   Ignored if `w_nex`, `mean_nex`, and `sd_nex` are provided.
#' @param tau_j A vector of positive numerics for the standard deviations \eqn{\tau_j}
#'   of the standard Nex priors. Ignored if `w_nex`, `mean_nex`, and `sd_nex` are provided.
#' @param w_j A numeric in `(0, 1)` for the weight of the Ex component, or a numeric vector
#'   of mixture weights summing to 1 if multiple Ex components are specified.
#' @param w_nex An optional numeric vector of mixture weights in \eqn{[0,1]} summing to 1
#'   for the Nex mixture prior.
#' @param mean_nex An optional numeric matrix of Nex mixture means with one row per Nex
#'   mixture component and one column per cohort.
#' @param sd_nex An optional positive numeric matrix of Nex mixture standard deviations with
#'   one row per Nex mixture component and one column per cohort.
#'
#' @return A list with prior parameters of class `prior_parameters_list`
#'
#' @details
#' This function sets the prior parameters for the method proposed by Neuenschwander et al. (2016).
#' If `w_nex`, `mean_nex`, and `sd_nex` are all `NULL`, the standard ExNex formulation is used.
#' Otherwise, the Nex part is specified as a finite mixture prior.
#'
#' @author Stephan Wojciekowski
#'
#' @examples
#' ## standard ExNex
#' prior_parameters_exnex <- setPriorParametersExNex(
#'   mu_mean   = 0,
#'   mu_sd     = 1,
#'   tau_scale = 2,
#'   mu_j      = c(4, 5),
#'   tau_j     = c(6, 7),
#'   w_j       = 0.8
#' )
#'
#' ## ExNex with Nex mixture prior
#' prior_parameters_exnex_mix <- setPriorParametersExNex(
#'   mu_mean   = 0,
#'   mu_sd     = 1,
#'   tau_scale = 2,
#'   w_j       = 0.8,
#'   w_nex     = c(0.7, 0.3),
#'   mean_nex  = rbind(c(4, 5), c(2, 3)),
#'   sd_nex    = rbind(c(6, 7), c(8, 9))
#' )
#'
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

  mu_j      = NULL,
  tau_j     = NULL,

  w_j,

  w_nex     = NULL,
  mean_nex  = NULL,
  sd_nex    = NULL

) {

  error_mu_mean    <- "Providing a numeric for the argument 'mu_mean'"
  error_mu_sd      <- "Providing a positive numeric for the argument 'mu_sd'"
  error_tau_scale  <- "Providing a positive numeric for the argument 'tau_scale'"
  error_mu_j       <- "Providing a (vector of) numeric(s) for the argument 'mu_j'"
  error_tau_j      <- "Providing a (vector of) positive numeric(s) for the argument 'tau_j'"
  error_w_j        <- "Providing a numeric in (0, 1) for the argument 'w_j'"
  error_w_nex      <- "Providing a numeric vector of weights in [0, 1] summing to 1 for the argument 'w_nex'"
  error_mean_nex   <- "Providing a numeric matrix for the argument 'mean_nex'"
  error_sd_nex     <- "Providing a positive numeric matrix for the argument 'sd_nex'"
  error_dim        <- "'mean_nex' and 'sd_nex' must have the same dimensions, and nrow(mean_nex) must equal length(w_nex)"
  error_mu_mean_sd <- "'mu_mean' and 'mu_sd' must have same length"
  error_mu_j_tau_j <- "'mu_j' and 'tau_j' must have the same length"
  error_w_j_long   <- "'w_j' must have length equal to length(mu_mean) + 1 if length(mu_mean) > 1"
  error_w_j_short  <- "'w_j' must have length 1 or 2 if length(mu_mean) = 1"
  error_w_j_sum    <- "Sum over items in 'w_j' must equal 1 if length(w_j) > 1"

  checkmate::assert_numeric(mu_mean, any.missing = FALSE, .var.name = error_mu_mean)
  checkmate::assert_numeric(mu_sd, lower = 0, any.missing = FALSE, .var.name = error_mu_sd)
  checkmate::assertTRUE(all(mu_sd > 0), .var.name = error_mu_sd)

  checkmate::assertNumber(tau_scale, .var.name = error_tau_scale)
  checkmate::assertTRUE(tau_scale > 0, .var.name = error_tau_scale)

  checkmate::assert_numeric(w_j, lower = 0, upper = 1, any.missing = FALSE, .var.name = error_w_j)

  checkmate::assert_true(length(mu_mean) == length(mu_sd), .var.name = error_mu_mean_sd)

  if (length(mu_mean) > 1) {
    checkmate::assert_true(length(w_j) == length(mu_mean) + 1, .var.name = error_w_j_long)
  } else {
    checkmate::assert_true(length(w_j) %in% c(1, 2), .var.name = error_w_j_short)
  }

  if (length(w_j) > 1) {
    checkmate::assert_true(isTRUE(all.equal(sum(w_j), 1)), .var.name = error_w_j_sum)
  }

  prior_base <- list(
    mu_mean   = mu_mean,
    mu_sd     = mu_sd,
    tau_scale = tau_scale,
    w_j       = w_j
  )

  if (is.null(w_nex) && is.null(mean_nex) && is.null(sd_nex)) {

    checkmate::assert_numeric(mu_j, any.missing = FALSE, .var.name = error_mu_j)
    checkmate::assert_numeric(tau_j, any.missing = FALSE, .var.name = error_tau_j)
    checkmate::assertTRUE(all(tau_j > 0), .var.name = error_tau_j)
    checkmate::assert_true(length(mu_j) == length(tau_j), .var.name = error_mu_j_tau_j)

    prior_extra <- list(
      mu_j  = mu_j,
      tau_j = tau_j
    )

  } else {

    checkmate::assertNumeric(w_nex, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w_nex)
    checkmate::assertTRUE(isTRUE(all.equal(sum(w_nex), 1)), .var.name = error_w_nex)

    checkmate::assertTRUE(is.matrix(mean_nex), .var.name = error_mean_nex)
    checkmate::assertTRUE(is.matrix(sd_nex), .var.name = error_sd_nex)
    checkmate::assertTRUE(all(sd_nex > 0), .var.name = error_sd_nex)
    checkmate::assertTRUE(all(dim(mean_nex) == dim(sd_nex)), .var.name = error_dim)
    checkmate::assertTRUE(nrow(mean_nex) == length(w_nex), .var.name = error_dim)

    prior_extra <- list(
      w_nex    = as.numeric(w_nex),
      mean_nex = mean_nex,
      sd_nex   = sd_nex
    )
  }

  prior_parameters_list <- list(exnex = c(prior_base, prior_extra))
  class(prior_parameters_list) <- "prior_parameters_list"

  return(prior_parameters_list)
}

## exnex_adj ####

getPriorParametersExNexAdj <- function (

  target_rates,
  tau_scale = 1,
  n_worth   = 1,

  w_j       = 0.5,

  w_nex     = NULL,
  mean_nex  = NULL,
  sd_nex    = NULL

) {

  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_tau_scale    <-
    "Providing a positive numeric for the argument 'tau_scale'"
  error_n_worth      <-
    "Providing a positive integer for the argument 'n_worth'"
  error_w_j          <-
    "Providing a numeric in (0, 1) for the argument 'w_j'"
  error_w_nex        <-
    "Providing a numeric vector of weights in [0, 1] summing to 1 for the argument 'w_nex'"
  error_mean_nex     <-
    "Providing a numeric matrix for the argument 'mean_nex'"
  error_sd_nex       <-
    "Providing a positive numeric matrix for the argument 'sd_nex'"
  error_dim          <-
    "'mean_nex' and 'sd_nex' must have the same dimensions, nrow(mean_nex) must equal length(w_nex), and ncol(mean_nex) must equal length(target_rates)"

  checkmate::assertNumeric(
    target_rates, any.missing = FALSE, .var.name = error_target_rates
  )
  checkmate::assertTRUE(
    all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates
  )

  checkmate::assertNumber(
    tau_scale, .var.name = error_tau_scale
  )
  checkmate::assertTRUE(
    tau_scale > 0, .var.name = error_tau_scale
  )

  checkmate::assertInt(
    n_worth, lower = 1, .var.name = error_n_worth
  )

  checkmate::assertNumeric(
    w_j, lower = 0, upper = 1, len = 1, .var.name = error_w_j
  )

  ## mixed Nex validation only if mixture inputs are supplied
  if (!(is.null(w_nex) && is.null(mean_nex) && is.null(sd_nex))) {
    checkmate::assertNumeric(
      w_nex, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w_nex
    )
    checkmate::assertTRUE(
      isTRUE(all.equal(sum(w_nex), 1)), .var.name = error_w_nex
    )

    checkmate::assertTRUE(
      is.matrix(mean_nex), .var.name = error_mean_nex
    )
    checkmate::assertTRUE(
      is.matrix(sd_nex), .var.name = error_sd_nex
    )
    checkmate::assertTRUE(
      all(sd_nex > 0), .var.name = error_sd_nex
    )
    checkmate::assertTRUE(
      all(dim(mean_nex) == dim(sd_nex)), .var.name = error_dim
    )
    checkmate::assertTRUE(
      nrow(mean_nex) == length(w_nex), .var.name = error_dim
    )
    checkmate::assertTRUE(
      ncol(mean_nex) == length(target_rates), .var.name = error_dim
    )
  }

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- getPriorParametersExNex(
    target_rates = target_rates,
    tau_scale    = tau_scale,
    n_worth      = n_worth,
    w_j          = w_j,
    w_nex        = w_nex,
    mean_nex     = mean_nex,
    sd_nex       = sd_nex
  )[[1]]

  prior_parameters$mu_mean <- 0

  if (is.null(w_nex) && is.null(mean_nex) && is.null(sd_nex)) {
    prior_parameters$mu_j <- rep(0, length(target_rates))
  }

  prior_parameters_list <- list(exnex_adj = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return(prior_parameters_list)
}

#' @title setPriorParametersExNexAdj
#' @md
#' @description This function sets prior parameters for the analysis method `"exnex_adj"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#'
#' It supports two specifications for the Nex part:
#' \itemize{
#'   \item the standard ExNex Adjusted specification with cohort-specific Nex priors via `mu_j` and `tau_j`
#'   \item an extended specification with a mixture prior on the Nex part via `w_nex`, `mean_nex`, and `sd_nex`
#' }
#'
#' @param mu_mean [numeric] Mean of \eqn{\mu}
#' @param mu_sd [numeric] Positive standard deviation of \eqn{\mu}
#' @param tau_scale [numeric] Positive scale parameter of \eqn{\tau}
#' @param mu_j [numeric] Vector of means \eqn{\mu_j} for the standard Nex priors.
#'   Ignored if `w_nex`, `mean_nex`, and `sd_nex` are provided.
#' @param tau_j [numeric] Vector of positive standard deviations \eqn{\tau_j}
#'   for the standard Nex priors. Ignored if `w_nex`, `mean_nex`, and `sd_nex` are provided.
#' @param w_j [numeric] Weight of the Ex component in `(0, 1)`, or a numeric vector
#'   of mixture weights summing to 1 if multiple Ex components are specified.
#' @param w_nex [numeric] Optional vector of mixture weights in \eqn{[0,1]} summing to 1
#'   for the Nex mixture prior.
#' @param mean_nex [numeric] Optional matrix of Nex mixture means with one row per Nex
#'   mixture component and one column per cohort.
#' @param sd_nex [numeric] Optional positive matrix of Nex mixture standard deviations with
#'   one row per Nex mixture component and one column per cohort.
#'
#' @return A list with prior parameters of class `prior_parameters_list`
#'
#' @details
#' This function sets prior parameters for the ExNex Adjusted method, which combines
#' the approach proposed by Neuenschwander et al. (2016) and the approach proposed by
#' Berry et al. (2013). If `w_nex`, `mean_nex`, and `sd_nex` are all `NULL`, the standard
#' ExNex Adjusted formulation is used. Otherwise, the Nex part is specified as a finite
#' mixture prior.
#'
#' @author Stephan Wojciekowski
#'
#' @examples
#' ## standard ExNex Adjusted
#' prior_parameters_exnex_adj <- setPriorParametersExNexAdj(
#'   mu_mean   = 0,
#'   mu_sd     = 1,
#'   tau_scale = 2,
#'   mu_j      = c(4, 5),
#'   tau_j     = c(6, 7),
#'   w_j       = 0.8
#' )
#'
#' ## ExNex Adjusted with Nex mixture prior
#' prior_parameters_exnex_adj_mix <- setPriorParametersExNexAdj(
#'   mu_mean   = 0,
#'   mu_sd     = 1,
#'   tau_scale = 2,
#'   w_j       = 0.8,
#'   w_nex     = c(0.7, 0.3),
#'   mean_nex  = rbind(c(0, 0), c(-1, -1)),
#'   sd_nex    = rbind(c(6, 7), c(8, 9))
#' )
#'
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
#' @references Neuenschwander, Beat, et al. "Robust exchangeability designs
#' for early phase clinical trials with multiple strata."
#' \emph{Pharmaceutical statistics} 15.2 (2016): 123-134.
#' @export
setPriorParametersExNexAdj <- function (

  mu_mean,
  mu_sd,
  tau_scale,

  mu_j      = NULL,
  tau_j     = NULL,

  w_j,

  w_nex     = NULL,
  mean_nex  = NULL,
  sd_nex    = NULL

) {

  error_mu_mean    <- "Providing a numeric for the argument 'mu_mean'"
  error_mu_sd      <- "Providing a positive numeric for the argument 'mu_sd'"
  error_tau_scale  <- "Providing a positive numeric for the argument 'tau_scale'"
  error_mu_j       <- "Providing a (vector of) numeric(s) for the argument 'mu_j'"
  error_tau_j      <- "Providing a (vector of) positive numeric(s) for the argument 'tau_j'"
  error_w_j        <- "Providing a numeric in (0, 1) for the argument 'w_j'"
  error_w_nex      <- "Providing a numeric vector of weights in [0, 1] summing to 1 for the argument 'w_nex'"
  error_mean_nex   <- "Providing a numeric matrix for the argument 'mean_nex'"
  error_sd_nex     <- "Providing a positive numeric matrix for the argument 'sd_nex'"
  error_dim        <- "'mean_nex' and 'sd_nex' must have the same dimensions, and nrow(mean_nex) must equal length(w_nex)"
  error_mu_mean_sd <- "'mu_mean' and 'mu_sd' must have same length"
  error_mu_j_tau_j <- "'mu_j' and 'tau_j' must have the same length"
  error_w_j_long   <- "'w_j' must have length equal to length(mu_mean) + 1 if length(mu_mean) > 1"
  error_w_j_short  <- "'w_j' must have length 1 or 2 if length(mu_mean) = 1"
  error_w_j_sum    <- "Sum over items in 'w_j' must equal 1 if length(w_j) > 1"

  checkmate::assert_numeric(
    mu_mean, any.missing = FALSE, .var.name = error_mu_mean
  )

  checkmate::assert_numeric(
    mu_sd, lower = 0, any.missing = FALSE, .var.name = error_mu_sd
  )
  checkmate::assertTRUE(
    all(mu_sd > 0), .var.name = error_mu_sd
  )

  checkmate::assertNumber(
    tau_scale, .var.name = error_tau_scale
  )
  checkmate::assertTRUE(
    tau_scale > 0, .var.name = error_tau_scale
  )

  checkmate::assert_numeric(
    w_j, lower = 0, upper = 1, any.missing = FALSE, .var.name = error_w_j
  )

  checkmate::assert_true(
    length(mu_mean) == length(mu_sd), .var.name = error_mu_mean_sd
  )

  if (length(mu_mean) > 1) {
    checkmate::assert_true(
      length(w_j) == length(mu_mean) + 1, .var.name = error_w_j_long
    )
  } else {
    checkmate::assert_true(
      length(w_j) %in% c(1, 2), .var.name = error_w_j_short
    )
  }

  if (length(w_j) > 1) {
    checkmate::assert_true(
      isTRUE(all.equal(sum(w_j), 1)), .var.name = error_w_j_sum
    )
  }

  ## standard ExNex Adjusted
  if (is.null(w_nex) && is.null(mean_nex) && is.null(sd_nex)) {

    checkmate::assert_numeric(
      mu_j, any.missing = FALSE, .var.name = error_mu_j
    )
    checkmate::assert_numeric(
      tau_j, any.missing = FALSE, .var.name = error_tau_j
    )
    checkmate::assertTRUE(
      all(tau_j > 0), .var.name = error_tau_j
    )
    checkmate::assert_true(
      length(mu_j) == length(tau_j), .var.name = error_mu_j_tau_j
    )

  } else {

    ## mixed Nex version
    checkmate::assertNumeric(
      w_nex, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w_nex
    )
    checkmate::assertTRUE(
      isTRUE(all.equal(sum(w_nex), 1)), .var.name = error_w_nex
    )

    checkmate::assertTRUE(
      is.matrix(mean_nex), .var.name = error_mean_nex
    )
    checkmate::assertTRUE(
      is.matrix(sd_nex), .var.name = error_sd_nex
    )
    checkmate::assertTRUE(
      all(sd_nex > 0), .var.name = error_sd_nex
    )
    checkmate::assertTRUE(
      all(dim(mean_nex) == dim(sd_nex)), .var.name = error_dim
    )
    checkmate::assertTRUE(
      nrow(mean_nex) == length(w_nex), .var.name = error_dim
    )
  }

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters_list <- setPriorParametersExNex(
    mu_mean   = mu_mean,
    mu_sd     = mu_sd,
    tau_scale = tau_scale,
    mu_j      = mu_j,
    tau_j     = tau_j,
    w_j       = w_j,
    w_nex     = w_nex,
    mean_nex  = mean_nex,
    sd_nex    = sd_nex
  )

  names(prior_parameters_list) <- "exnex_adj"

  return(prior_parameters_list)
}


## pooled ####

getPriorParametersPooled <- function (

  target_rates,
  n_worth = 1

) {

  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_n_worth      <-
    "Providing a positive integer for the argument 'n_worth'"

  checkmate::assertNumeric(
    target_rates, any.missing = FALSE, .var.name = error_target_rates
  )
  checkmate::assertTRUE(
    all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates
  )

  checkmate::assertInt(
    n_worth, lower = 1, .var.name = error_n_worth
  )

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

  error_a <- "Providing a positive numeric for the argument 'a'"
  error_b <- "Providing a positive numeric for the argument 'b'"

  checkmate::assertNumber(
    a, .var.name = error_a
  )
  checkmate::assertTRUE(
    a > 0, .var.name = error_a
  )

  checkmate::assertNumber(
    b, .var.name = error_b
  )
  checkmate::assertTRUE(
    b > 0, .var.name = error_b
  )

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

  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_n_worth      <-
    "Providing a positive integer for the argument 'n_worth'"

  checkmate::assertNumeric(
    target_rates, any.missing = FALSE, .var.name = error_target_rates
  )
  checkmate::assertTRUE(
    all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates
  )

  checkmate::assertInt(
    n_worth, lower = 1, .var.name = error_n_worth
  )

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

  error_a_j          <- "Providing a (vector of) positive numeric(s) for the argument 'a_j'"
  error_b_j          <- "Providing a (vector of) positive numeric(s) for the argument 'b_j'"
  error_a_j_b_j_len  <- "a_j and b_j must have the same length"

  checkmate::assertNumeric(
    a_j, any.missing = FALSE, .var.name = error_a_j
  )
  checkmate::assertTRUE(
    all(a_j > 0), .var.name = error_a_j
  )

  checkmate::assertNumeric(
    b_j, any.missing = FALSE, .var.name = error_b_j
  )
  checkmate::assertTRUE(
    all(b_j > 0), .var.name = error_b_j
  )

  checkmate::assertTRUE(
    length(a_j) == length(b_j), .var.name = error_a_j_b_j_len
  )

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  prior_parameters <- list(a_j = a_j, b_j = b_j)

  prior_parameters_list        <- list(stratified = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

## stratified_mix ####

getPriorParametersStratifiedMix <- function (

  target_rates,
  n_worth = 1,
  w       = c(0.8, 0.2)

) {

  error_target_rates <-
    "Providing a vector of numerics in (0, 1) for the argument 'target_rates'"
  error_n_worth      <-
    "Providing a positive integer for the argument 'n_worth'"
  error_w            <-
    "Providing a numeric vector of weights in [0, 1] summing to 1 for the argument 'w'"

  checkmate::assertNumeric(
    target_rates, any.missing = FALSE, .var.name = error_target_rates
  )
  checkmate::assertTRUE(
    all(target_rates > 0 & target_rates < 1), .var.name = error_target_rates
  )

  checkmate::assertInt(
    n_worth, lower = 1, .var.name = error_n_worth
  )

  checkmate::assertNumeric(
    w, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w
  )
  checkmate::assertTRUE(
    isTRUE(all.equal(sum(w), 1)), .var.name = error_w
  )

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  K <- length(w)
  J <- length(target_rates)

  if (K != 2) {
    stop(simpleError("The default implementation of 'getPriorParametersStratifiedMix()' currently expects a weight vector of length 2."))
  }

  a_j <- matrix(NA_real_, nrow = K, ncol = J)
  b_j <- matrix(NA_real_, nrow = K, ncol = J)

  ## component 1: informative stratified prior
  a_j[1, ] <- target_rates * n_worth
  b_j[1, ] <- (1 - target_rates) * n_worth

  ## component 2: vague prior
  a_j[2, ] <- rep(1, J)
  b_j[2, ] <- rep(1, J)

  prior_parameters <- list(
    w   = w,
    a_j = a_j,
    b_j = b_j
  )

  prior_parameters_list <- list(stratified_mix = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return (prior_parameters_list)

}

#' @title setPriorParametersStratifiedMix
#' @md
#' @description This function sets prior parameters for the analysis method `"stratified_mix"`
#' for use in \code{\link[bhmbasket]{performAnalyses}}.
#' @param w A numeric vector of mixture weights in \eqn{[0,1]} summing to 1.
#' @param a_j A positive numeric matrix of beta shape parameters \eqn{\alpha_j},
#' with rows per mixture component and columns per cohort.
#' @param b_j A positive numeric matrix of beta shape parameters \eqn{\beta_j},
#' with rows per mixture component and columns per cohort.
#' @return A list with prior parameters of class `prior_parameters_list`
#' @details
#' The method `"stratified_mix"` is a beta-binomial model that assesses each cohort
#' individually with a finite mixture beta prior.
#' See also the R package `RBesT`.
#' @author Stephan Wojciekowski
#' @examples
#' prior_parameters_stratified_mix <- setPriorParametersStratifiedMix(
#'   w   = c(0.8, 0.2),
#'   a_j = rbind(c(2, 3), c(1, 1)),
#'   b_j = rbind(c(8, 7), c(1, 1))
#' )
#' @rdname setPriorParametersStratifiedMix
#' @seealso
#'  \code{\link[bhmbasket]{performAnalyses}}
#'  \code{\link[bhmbasket]{getPriorParameters}}
#'  \code{\link[bhmbasket]{combinePriorParameters}}
#'  \code{\link[bhmbasket]{setPriorParametersStratified}}
#' @export
setPriorParametersStratifiedMix <- function(
    w,
    a_j,
    b_j
) {

  error_w    <- "Providing a numeric vector of weights in [0,1] summing to 1 for argument 'w'"
  error_a_j  <- "Providing a positive numeric matrix for argument 'a_j'"
  error_b_j  <- "Providing a positive numeric matrix for argument 'b_j'"
  error_dim  <- "a_j and b_j must have the same dimensions, and nrow(a_j) must equal length(w)"

  checkmate::assertNumeric(w, any.missing = FALSE, lower = 0, upper = 1, .var.name = error_w)
  checkmate::assertTRUE(isTRUE(all.equal(sum(w), 1)), .var.name = error_w)

  checkmate::assertTRUE(is.matrix(a_j), .var.name = error_a_j)
  checkmate::assertTRUE(is.matrix(b_j), .var.name = error_b_j)
  checkmate::assertTRUE(all(a_j > 0), .var.name = error_a_j)
  checkmate::assertTRUE(all(b_j > 0), .var.name = error_b_j)

  checkmate::assertTRUE(all(dim(a_j) == dim(b_j)), .var.name = error_dim)
  checkmate::assertTRUE(nrow(a_j) == length(w), .var.name = error_dim)

  prior_parameters <- list(
    w   = as.numeric(w),
    a_j = a_j,
    b_j = b_j
  )

  prior_parameters_list <- list(stratified_mix = prior_parameters)
  class(prior_parameters_list) <- "prior_parameters_list"

  return(prior_parameters_list)
}
