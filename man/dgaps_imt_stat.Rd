% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgaps_imt.R
\name{dgaps_imt_stat}
\alias{dgaps_imt_stat}
\title{Statistics for the \eqn{D}-gaps information matrix test}
\usage{
dgaps_imt_stat(data, theta, u, D = 1, inc_cens = TRUE)
}
\arguments{
\item{data}{A numeric vector of raw data.  Missing values are allowed, but
they should not appear between non-missing values, that is, they only be
located at the start and end of the vector.  Missing values are omitted
using \code{\link[stats]{na.omit}}.}

\item{theta}{A numeric scalar. An estimate of the extremal index
\eqn{\theta}, produced by \code{\link{dgaps}}.}

\item{u}{A numeric scalar.  Extreme value threshold applied to data.}

\item{D}{A numeric scalar.  The censoring parameter \eqn{D}. Threshold
inter-exceedances times that are not larger than \code{D} units are
left-censored, occurring with probability
\eqn{\log(1 - \theta e^{-\theta d})}{log(1 - exp(-\theta d))},
where \eqn{d = q D} and \eqn{q} is the probability with which the
threshold \eqn{u} is exceeded.}

\item{inc_cens}{A logical scalar indicating whether or not to include
contributions from right-censored inter-exceedance times, relating to the
first and last observations. See \code{\link{dgaps}}.}
}
\value{
A list
relating the quantities given on pages 18-19 of
  Suveges and Davison (2010).  All but the last component are vectors giving
  the contribution to the quantity from each \eqn{D}-gap, evaluated at the
  input value \code{theta} of \eqn{\theta}.
  \item{\code{ldj} }{the derivative of the log-likelihood with respect to
    \eqn{\theta} (the score)}
  \item{\code{Ij} }{the observed information}
  \item{\code{Jj} }{the square of the score}
  \item{\code{dj} }{\code{Jj} - \code{Ij}}
  \item{\code{Ddj} }{the derivative of \code{Jj} - \code{Ij} with respect
    to \eqn{\theta}}
  \item{\code{n_dgaps} }{the number of \eqn{D}-gaps that contribute to the
    log-likelihood.}
}
\description{
Calculates the components required to calculate the value of the information
matrix test under the \eqn{D}-gaps model, using vector data input.
Called by \code{\link{dgaps_imt}}.
}
\references{
Holesovsky, J. and Fusek, M. Estimation of the extremal index
  using censored distributions. Extremes 23, 197-213 (2020).
  \doi{10.1007/s10687-020-00374-3}
}
