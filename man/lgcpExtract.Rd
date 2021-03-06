% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{lgcpExtract}
\alias{lgcpExtract}
\title{Extract samples from \code{lgcp} model}
\usage{
lgcpExtract(dirname, nchains, verbose = TRUE)
}
\arguments{
\item{dirname}{Rootname of the directory in call to \code{lgcp}.}

\item{nchains}{Integer, number of chains from call to \code{lgcp}.}

\item{verbose}{Logical indicating whether to display progress bar}
}
\value{
An array of dimension number of gridcells x number of iterations x
number of time periods.
}
\description{
Extract MCMC samples from a call to \code{lgcp}
}
