% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PriorFunctions.R
\name{setPriorParametersPooled}
\alias{setPriorParametersPooled}
\title{setPriorParametersPooled}
\usage{
setPriorParametersPooled(a, b)
}
\arguments{
\item{a}{A positive numeric for \eqn{\alpha}}

\item{b}{A positive numeric for \eqn{\beta}}
}
\value{
A list with prior parameters of class \code{prior_parameters_list}
}
\description{
This function sets prior parameters for the analysis method \code{"pooled"}
for use in \code{\link[bhmbasket]{performAnalyses}}.
}
\details{
The method \code{"pooled"} is a beta-binomial model that pools all cohorts.
The prior parameters are the scale parameters of the beta prior distribution.
}
\examples{
 prior_parameters_pooled <- setPriorParametersPooled(1, 2)
}
\seealso{
\code{\link[bhmbasket]{performAnalyses}}
\code{\link[bhmbasket]{getPriorParameters}}
\code{\link[bhmbasket]{combinePriorParameters}}
\code{\link[bhmbasket]{setPriorParametersBerry}}
\code{\link[bhmbasket]{setPriorParametersExNex}}
\code{\link[bhmbasket]{setPriorParametersExNexAdj}}
\code{\link[bhmbasket]{setPriorParametersStratified}}
\code{\link[bhmbasket]{getMuVar}}
}
\author{
Stephan Wojciekowski
}
