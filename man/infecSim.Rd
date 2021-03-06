% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infec_sim.R
\name{infecSim}
\alias{infecSim}
\title{Infection simulation}
\usage{
infecSim(
  region,
  t.win,
  covariates,
  mean.val,
  p = 1/mean.val,
  delta,
  rho,
  beta = c(1, 1),
  t.off = 0,
  cov.pars = c(1, 0.0075, 0.5),
  grid.size = 64^2
)
}
\arguments{
\item{region}{A \code{spatialPolygon} defining the area to simulate data for}

\item{t.win}{A vector indicating the time window to simulate data for, eg. c(1, 30)}

\item{covariates}{A \code{spatialPolygonsDataFrame} containing population density information for the area of interest. The population
density variable should be called \code{popdens}.}

\item{mean.val}{Integer, the mean number of cases per time period}

\item{p}{Probability a case generates a cluster}

\item{delta}{A vector of two values: the spatial range and temporal range parameters}

\item{rho}{Multiplicative factor of effect of a cluster}

\item{beta}{Parameters of latent field covariates}

\item{t.off}{Temporal offset parameter - number of time periods to displace the peak of the cluster intensity.}

\item{cov.pars}{Covariance parameters of latent Gaussian field covariates: partial sill, range, and nugget}

\item{grid.size}{Size of the computational grid to simulate the infectious process}
}
\value{
A list: (1) simulated data for each computational grid cell, (2) simulated case locations and time, (3) plot
of Poisson intensity, (4) plot of simulated case locations, (5) stpp object of case locations for use with stpp functions.
}
\description{
Simulate data from an infectious disease process modelled as a Poisson process
}
\details{
\enumerate{
\item Background incidence is a Poisson process with intensity lambda x (pop density)
\item A small fraction p of background events trigger a temporary increase in local
intensity to lambda x (pop density) x (1+rho) in a disc of  radius delta for the next k time-periods.
}
}
\examples{
data(square,square_pop)
infecSim(region = square,
         t.win = c(1,10),
         covariates = square_pop,
         mean.val= 100,
         p =1/100,
         delta = c(0.01,4),
         rho=3,
         t.off = 4,
         cov.pars = c(0.9,0.03,0.1),
         grid.size = 64^2)
}
