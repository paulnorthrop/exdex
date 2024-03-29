% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kgaps.R
\name{kgaps}
\alias{kgaps}
\title{Maximum likelihood estimation for the \eqn{K}-gaps model}
\usage{
kgaps(data, u, k = 1, inc_cens = TRUE)
}
\arguments{
\item{data}{A numeric vector or numeric matrix of raw data.  If \code{data}
  is a matrix then the log-likelihood is constructed as the sum of
  (independent) contributions from different columns. A common situation is
  where each column relates to a different year.

  If \code{data} contains missing values then \code{\link{split_by_NAs}} is
  used to divide the data further into sequences of non-missing values,
  stored in different columns in a matrix.  Again, the log-likelihood
  is constructed as a sum of contributions from different columns.}

\item{u}{A numeric scalar.  Extreme value threshold applied to data.}

\item{k}{A non-negative numeric scalar. Run parameter \eqn{K}, as defined in
Suveges and Davison (2010).  Threshold inter-exceedances times that are
not larger than \code{k} units are assigned to the same cluster, resulting
in a \eqn{K}-gap equal to zero. Specifically, the \eqn{K}-gap \eqn{S}
corresponding to an inter-exceedance time of \eqn{T} is given by
\eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.  In practice, \eqn{k} should
be no smaller than 1, because when \eqn{k} is less than 1 the estimate
of \eqn{\theta} is always equal to 1.}

\item{inc_cens}{A logical scalar indicating whether or not to include
contributions from right-censored inter-exceedance times, relating to the
first and last observations.  It is known that these times are greater
than or equal to the time observed. See Attalides (2015) for details.
If \code{data} has multiple columns then there will be right-censored
first and last inter-exceedance times for each column.}
}
\value{
An object (a list) of class \code{c("kgaps", "exdex")} containing
    \item{\code{theta} }{The maximum likelihood estimate (MLE) of
      \eqn{\theta}.}
    \item{\code{se} }{The estimated standard error of the MLE, calculated
      using an algebraic expression for the observed information.
      If \code{k = 0} then \code{se} is returned as \code{0}.}
    \item{\code{se_exp} }{The estimated standard error of the MLE,
      calculated using an algebraic expression for the expected information.
      If the estimate of \eqn{\theta} is 0 or 1 then \code{se_exp} is
      \code{NA}.}
    \item{\code{ss} }{The list of summary statistics returned from
      \code{\link{kgaps_stat}}.}
    \item{\code{k, u, inc_cens} }{The input values of \code{k},
      \code{u} and \code{inc_cens}.}
    \item{\code{max_loglik }}{The value of the log-likelihood at the MLE.}
    \item{\code{call }}{The call to \code{kgaps}.}
}
\description{
Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
based on the \eqn{K}-gaps model for threshold inter-exceedances times of
Suveges and Davison (2010).
}
\details{
If \code{inc_cens = FALSE} then the maximum likelihood estimate of
  the extremal index \eqn{\theta} under the \eqn{K}-gaps model of
  Suveges and Davison (2010) is calculated.

  If \code{inc_cens = TRUE} then information from right-censored
  first and last inter-exceedance times is also included in the likelihood
  to be maximized, following Attalides (2015).  The form of the
  log-likelihood is given in the \strong{Details} section of
  \code{\link{kgaps_stat}}.

  It is possible that the estimate of \eqn{\theta} is equal to 1, and also
  possible that it is equal to 0. \code{\link{kgaps_stat}} explains the
  respective properties of the data that cause these events to occur.
}
\examples{
### S&P 500 index

u <- quantile(sp500, probs = 0.60)
theta <- kgaps(sp500, u)
theta
summary(theta)
coef(theta)
nobs(theta)
vcov(theta)
logLik(theta)

### Newlyn sea surges

u <- quantile(newlyn, probs = 0.60)
theta <- kgaps(newlyn, u, k = 2)
theta
summary(theta)

### Cheeseboro wind gusts

theta <- kgaps(cheeseboro, 45, k = 3)
theta
summary(theta)
}
\references{
Suveges, M. and Davison, A. C. (2010) Model
  misspecification in peaks over threshold analysis, \emph{Annals of
  Applied Statistics}, \strong{4}(1), 203-221.
  \doi{10.1214/09-AOAS292}

Attalides, N. (2015) Threshold-based extreme value modelling,
  PhD thesis, University College London.
  \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
}
\seealso{
\code{\link{kgaps_confint}} to estimate confidence intervals
  for \eqn{\theta}.

\code{\link{kgaps_methods}} for S3 methods for \code{"kgaps"}
  objects.

\code{\link{kgaps_imt}} for the information matrix test, which
  may be used to inform the choice of the pair (\code{u, k}).

\code{\link{choose_uk}} for a diagnostic plot based on
  \code{\link{kgaps_imt}}.

\code{\link{kgaps_stat}} for the calculation of sufficient
  statistics for the \eqn{K}-gaps model.

\code{\link[revdbayes]{kgaps_post}} in the
  \code{\link[revdbayes]{revdbayes}} package for Bayesian inference
  about \eqn{\theta} using the \eqn{K}-gaps model.
}
