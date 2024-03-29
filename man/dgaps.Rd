% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgaps.R
\name{dgaps}
\alias{dgaps}
\title{Maximum likelihood estimation using left-censored inter-exceedances times}
\usage{
dgaps(data, u, D = 1, inc_cens = TRUE)
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

\item{D}{A numeric scalar.  The censoring parameter \eqn{D}. Threshold
inter-exceedances times that are not larger than \code{D} units are
left-censored, occurring with probability
\eqn{\log(1 - \theta e^{-\theta d})}{log(1 - exp(-\theta d))},
where \eqn{d = q D} and \eqn{q} is the probability with which the
threshold \eqn{u} is exceeded.}

\item{inc_cens}{A logical scalar indicating whether or not to include
contributions from right-censored inter-exceedance times, relating to the
first and last observations. It is known that these times are greater
than or equal to the time observed.
If \code{data} has multiple columns then there will be right-censored
first and last inter-exceedance times for each column.}
}
\value{
An object (a list) of class \code{c("dgaps", "exdex")} containing
    \item{\code{theta} }{The maximum likelihood estimate (MLE) of
      \eqn{\theta}.}
    \item{\code{se} }{The estimated standard error of the MLE, calculated
      using an algebraic expression for the observed information.  If the
      estimate of \eqn{\theta} is 0 then \code{se} is \code{NA}.}
    \item{\code{se_exp} }{The estimated standard error of the MLE,
      calculated using an algebraic expression for the expected information.
      If the estimate of \eqn{\theta} is 0 then \code{se_exp} is \code{NA}.
      This is provided because cases may be encountered where the observed
      information is not positive.}
    \item{\code{ss} }{The list of summary statistics returned from
      \code{\link{dgaps_stat}}.}
    \item{\code{D, u, inc_cens} }{The input values of \code{D},
      \code{u} and \code{inc_cens}.}
    \item{\code{max_loglik }}{The value of the log-likelihood at the MLE.}
    \item{\code{call }}{The call to \code{dgaps}.}
}
\description{
Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
based on a model for threshold inter-exceedances times of
Holesovsky and Fusek (2020).  We refer to this as the \eqn{D}-gaps model,
because it uses a tuning parameter \eqn{D}, whereas the related \eqn{K}-gaps
model of Suveges and Davison (2010) has a tuning parameter \eqn{K}.
}
\details{
If \code{inc_cens = FALSE} then the maximum likelihood estimate of
  the extremal index \eqn{\theta} under the \eqn{D}-gaps model of
  Holesovsky and Fusek (2020) is calculated.  Under this model
  inter-exceedance times that are less than or equal to \eqn{D} are
  left-censored, as a strategy to mitigate model mis-specification resulting
  from the fact that inter-exceedance times that are equal to 0 are expected
  under the model but only positive inter-exceedance times can be observed
  in practice.

  If \code{inc_cens = TRUE} then information from the right-censored
  first and last inter-exceedance times are also included in the likelihood
  to be maximized.
  For an explanation of the idea see Attalides (2015).  The form of the
  log-likelihood is given in the \strong{Details} section of
  \code{\link{dgaps_stat}}.

  It is possible that the estimate of \eqn{\theta} is equal to 1, and also
  possible that it is equal to 0. \code{\link{dgaps_stat}} explains the
  respective properties of the data that cause these events to occur.
}
\examples{
### S&P 500 index

u <- quantile(sp500, probs = 0.60)
theta <- dgaps(sp500, u = u, D = 1)
theta
summary(theta)
coef(theta)
nobs(theta)
vcov(theta)
logLik(theta)

### Newlyn sea surges

u <- quantile(newlyn, probs = 0.60)
theta <- dgaps(newlyn, u = u, D = 2)
theta
summary(theta)

### Uccle July temperatures

# Using vector input, which merges data from different years
u <- quantile(uccle720$temp, probs = 0.9, na.rm = TRUE)
theta <- dgaps(uccle720$temp, u = u, D = 2)
theta

# Using matrix input to separate data from different years
u <- quantile(uccle720m, probs = 0.9, na.rm = TRUE)
theta <- dgaps(uccle720m, u = u, D = 2)
theta
}
\references{
Holesovsky, J. and Fusek, M. Estimation of the extremal index
  using censored distributions. Extremes 23, 197-213 (2020).
  \doi{10.1007/s10687-020-00374-3}

Suveges, M. and Davison, A. C. (2010) Model
  misspecification in peaks over threshold analysis, \emph{Annals of
  Applied Statistics}, \strong{4}(1), 203-221.
  \doi{10.1214/09-AOAS292}
}
\seealso{
\code{\link{dgaps_confint}} to estimate confidence intervals
  for \eqn{\theta}.

\code{\link{dgaps_methods}} for S3 methods for \code{"dgaps"}
  objects.

\code{\link{dgaps_imt}} for the information matrix test, which
  may be used to inform the choice of the pair (\code{u, D}).

\code{\link{choose_ud}} for a diagnostic plot based on
  \code{\link{dgaps_imt}}.

\code{\link{dgaps_stat}} for the calculation of sufficient
  statistics for the \eqn{D}-gaps model.
}
