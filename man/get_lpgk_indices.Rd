% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_method.R
\name{get_lpgk_indices}
\alias{get_lpgk_indices}
\title{Identifies the longest path in the simulation}
\usage{
get_lpgk_indices(par_loc, cell_params)
}
\arguments{
\item{par_loc}{location of the PROSSTT parameter file}

\item{cell_params}{a data frame that includes branch and pseudotime labels}
}
\value{
the indices of the cells on the longest path
}
\description{
Finds the indices of the cells on the longest path in the simulated lineage
tree.
}
