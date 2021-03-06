% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_method.R
\name{read_write_output}
\alias{read_write_output}
\title{Appends a result to a file}
\usage{
read_write_output(res, out)
}
\arguments{
\item{res}{the result matrix (produced by `evaluate_method()`)}

\item{out}{the location of the file to append to}
}
\value{
None
}
\description{
Appends an one-row matrix to the end of a file with a matrix of the same
column number (with same column names). Useful for benchmarking multiple
methods and keeping all results in one file.
}
