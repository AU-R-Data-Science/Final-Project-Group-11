% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_functions.R
\name{metricsplot}
\alias{metricsplot}
\title{Metrics Plot}
\usage{
metricsplot(metric)
}
\arguments{
\item{metric}{A \code{character} assigning which metric to plot.}
}
\value{
A plot showing the desired matrix on a grid.
#' @author Rukesh Gusain, Michael Zirpoli, Erica Maul
}
\description{
Plot a user-assigned metric on a grid. The following metrics can be inputted as characters and, therefore, calculated: `Prevalence`, `Accuracy`, `Sensitivity`, `Specificity`, `False_Discovery_Rate`, and `Diagnostic_Odds_Ratio`.
}
