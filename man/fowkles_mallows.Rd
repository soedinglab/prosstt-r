% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{fowkles_mallows}
\alias{fowkles_mallows}
\title{Calculates the Fowkles-Mallows index.}
\usage{
fowkles_mallows(status)
}
\arguments{
\item{status}{a truth table}
}
\value{
the Fowkles-Mallows index of the input truth table
}
\description{
Calculates the Fowkles-Mallows index:
}
\details{
\deqn{\sqrt{\frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FP}} \frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FN}} }}
}
