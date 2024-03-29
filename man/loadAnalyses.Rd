% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysisFunctions.R
\name{loadAnalyses}
\alias{loadAnalyses}
\title{loadAnalyses}
\usage{
loadAnalyses(
  scenario_numbers,
  analysis_numbers = rep(1, length(scenario_numbers)),
  load_path = tempdir()
)
}
\arguments{
\item{scenario_numbers}{A (vector of) positive integer(s) for the scenario number(s)}

\item{analysis_numbers}{A (vector of) positive integer(s) for the analysis number(s),
Default: \code{rep(1, length(scenario_numbers))}}

\item{load_path}{A string providing a path where the scenarios are being stored,
Default: \code{\link[base]{tempfile}}}
}
\value{
Returns an object of class \code{analysis_list}
}
\description{
This function loads an analysis performed with
\code{\link[bhmbasket]{performAnalyses}}
}
\examples{
  trial_data <- createTrial(
    n_subjects   = c(10, 20, 30),
    n_responders = c(1, 2, 3))

  analysis_list <- performAnalyses(
    scenario_list      = trial_data,
    target_rates       = rep(0.5, 3),
    n_mcmc_iterations  = 100)

  save_info     <- saveAnalyses(analysis_list)
  analysis_list <- loadAnalyses(scenario_numbers = save_info$scenario_numbers,
                                analysis_numbers = save_info$analysis_numbers,
                                load_path        = save_info$path)
}
\seealso{
\code{\link[bhmbasket]{performAnalyses}}
\code{\link[bhmbasket]{saveAnalyses}}
\code{\link[base]{tempfile}}
}
\author{
Stephan Wojciekowski
}
