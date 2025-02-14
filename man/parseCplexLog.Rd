% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solvers_cplex__parse_cplex_log.R
\name{parseCplexLog}
\alias{parseCplexLog}
\title{Parses the cplex log file and reads some basic information.}
\usage{
parseCplexLog(log)
}
\arguments{
\item{log}{path of log file resulted from a carnival run OR the content
of this file read by \code{\link[readr]{read_lines}}.}
}
\value{
list variable with following fields:
- `convergence` a table that contains information on the convergence of CPLEX 
- `n_solutions` number of solutions found 
- `objective` objective function value
- `termination_reason`: reason of termination
}
\description{
Parses the cplex log file and reads some basic information.
}
\author{
Attila Gabor, 2021
}
