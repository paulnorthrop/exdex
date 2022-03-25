#' exdex: Estimation of the Extremal Index
#'
#' The extremal index \eqn{\theta} is a measure of the degree of local
#' dependence in the extremes of a stationary process.  The \emph{exdex}
#' package  performs frequentist inference about \eqn{\theta} using the
#' methodologies proposed in Northrop (2015), Berghaus and Bucher (2018),
#' Suveges (2007) and Suveges and Davison (2010).
#'
#' @details Functions to implement three estimators of the extremal index
#'   are provided, namely
#' \itemize{
#'   \item{\code{\link{spm}}: semiparametric maxima estimator, using block
#'     maxima: (Northrop, 2015; Berghaus and Bucher, 2018)}
#'   \item{\code{\link{kgaps}}: \eqn{K}-gaps estimator, using threshold
#'     interexceedance times (Suveges and Davison, 2010)}
#'   \item{\code{\link{iwls}}: iterated weighted least squares estimator,
#'     using threshold interexceedance times: (Suveges, 2007)}
#' }
#' The functions \code{\link{choose_b}} and \code{\link{choose_uk}} provide
#' graphical diagnostics for choosing the tuning parameter for the
#' semiparametric estimator, the block size \eqn{b}, and the tuning parameters
#' of the \eqn{K}-gaps estimator, the threshold \eqn{u} and the run parameter
#' \eqn{K}.
#'
#' For the \eqn{K}-gaps model the `exdex` package allows missing values in the
#' data, can accommodate independent subsets of data, such as monthly or
#' seasonal time series from different years, and can incorporate information
#' from censored interexceedance times.
#'
#' See \code{vignette("exdex-vignette", package = "exdex")} for an
#' overview of the package.
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#'   maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#'   \strong{46}(5), 2307-2335. \doi{10.1214/17-AOS1621}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#'   \doi{10.1007/s10687-015-0221-5}
#' @references Suveges, M. (2007) Likelihood estimation of the extremal
#'   index. \emph{Extremes}, \strong{10}, 41-55.
#'   \doi{10.1007/s10687-007-0034-2}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{spm}}: semiparametric maxima estimator.
#' @seealso \code{\link{kgaps}}: \eqn{K}-gaps estimator.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @seealso \code{\link{choose_b}} and \code{\link{choose_uk}} for choosing
#'   tuning parameters.
#' @seealso \code{\link{newlyn}}, \code{\link{sp500}} and
#'   \code{\link{cheeseboro}} for example datasets.
#' @docType package
#' @name exdex
#' @import methods
#' @importFrom graphics plot
#' @importFrom stats coef confint nobs vcov
#' @useDynLib exdex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Newlyn sea surges
#'
#' The vector \code{newlyn} contains 2894 maximum sea-surges measured at
#' Newlyn, Cornwall, UK over the period 1971-1976. The observations are
#' the maximum hourly sea-surge heights over contiguous 15-hour time
#' periods.
#' @format A vector of length 2894.
#' @source Coles, S.G. (1991) Modelling extreme multivariate events. PhD thesis,
#'   University of Sheffield, U.K.
#' @references Fawcett, L. and Walshaw, D. (2012) Estimating return levels from
#'   serially dependent extremes. \emph{Environmetrics}, \strong{23}(3),
#'   272-283.  \doi{10.1002/env.2133}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index. \emph{Extremes}, \strong{18},
#'   585-603.  \doi{10.1007/s10687-015-0221-5}
"newlyn"

#' Daily log returns of the Standard and Poor (S&P) 500 index
#'
#' Daily log returns of the S&P 500 index, that is, the log of the ratio of
#' successive daily closing prices, from 3rd January 1990 to 9th October 2018.
#'
#' @format A vector of length 7250, created using \code{\link[zoo]{zoo}}
#'   with an "index" attribute giving the date of the corresponding negated
#'   log return.
#' @source Yahoo finance: https://finance.yahoo.com/quote/^SPX/history/
"sp500"

#' Cheeseboro hourly maximum wind gusts
#'
#' The matrix \code{cheeseboro} contains hourly maximum wind gusts recorded at
#' the Cheeseboro weather station near Thousand Oaks, Southern California, USA
#' during the month of January over the period 2000-2009. These data are
#' analysed in Reich and Shaby (2016).
#' @format A 774 by 10 numeric matrix.  Column \code{i} contains the hourly
#'   maximum wind gusts from Cheeseboro in the year 2000 + \code{i} - 1.
#'   The columns are named 2000, 2001, ..., 2009 and the rows are named
#'   day\code{j}hour\code{k}, where \code{j} is the day of the month and
#'   \code{k} the hour of the day.
#' @note There are a total of 42 missing values, located in 6 of the 10 years,
#'   namely 2000-2003 and 2005-2006.
#' @source The Remote Automated Weather Stations USA Climate Archive at
#'   \url{https://raws.dri.edu/}, more specifically the Daily Summaries of the
#'   \href{https://raws.dri.edu/cgi-bin/rawMAIN.pl?caCCHB}{Cheeseboro page}.
#' @references Reich, B. J. and Shaby, B. A. (2016). 'Time series of Extremes',
#' in Dey, D. K. and Yan, J. (eds.) Extreme Value Modeling and Risk Analysis.
#' New York: Chapman and Hall/CRC, pp. 131-151.
"cheeseboro"
