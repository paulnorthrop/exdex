#' exdex: Estimation of the Extremal Index
#'
#' The extremal index \eqn{\theta} is a measure of the degree of local
#' dependence in the extremes of a stationary process.  The \emph{exdex}
#' package  performs frequentist inference about \eqn{\theta} using the
#' methodologies proposed in Northrop (2015), Berghaus and Bucher (2018),
#' Suveges (2007), Suveges and Davison (2010) and Holesovsky and Fusek (2020).
#'
#' @details Functions to implement four estimators of the extremal index
#'   are provided, namely
#' \itemize{
#'   \item{\code{\link{spm}}: semiparametric maxima estimator, using block
#'     maxima: (Northrop, 2015; Berghaus and Bucher, 2018)}
#'   \item{\code{\link{kgaps}}: \eqn{K}-gaps estimator, using threshold
#'     inter-exceedance times (Suveges and Davison, 2010)}
#'   \item{\code{\link{dgaps}}: \eqn{D}-gaps estimator, using threshold
#'     inter-exceedance times (Holesovsky and Fusek, 2020))}
#'   \item{\code{\link{iwls}}: iterated weighted least squares estimator,
#'     using threshold inter-exceedance times: (Suveges, 2007)}
#' }
#' The functions \code{\link{choose_b}}, \code{\link{choose_uk}} and
#' \code{\link{choose_ud}} provide graphical diagnostics for choosing the
#' respective tuning parameters of the semiparametric maxima, \eqn{K}-gaps and
#' \eqn{D}-gaps estimators.
#'
#' For the \eqn{K}-gaps and \eqn{D}-gaps models the `exdex` package allows
#' missing values in the data, can accommodate independent subsets of data,
#' such as monthly or seasonal time series from different years, and can
#' incorporate information from censored inter-exceedance times.
#'
#' See \code{vignette("exdex-vignette", package = "exdex")} for an
#' overview of the package.
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#'   maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#'   \strong{46}(5), 2307-2335. \doi{10.1214/17-AOS1621}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#'   \doi{10.1007/s10687-015-0221-5}
#' @references Suveges, M. (2007) Likelihood estimation of the extremal
#'   index. \emph{Extremes}, \strong{10}, 41-55.
#'   \doi{10.1007/s10687-007-0034-2}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{spm}}: semiparametric maxima estimator.
#' @seealso \code{\link{kgaps}}: \eqn{K}-gaps estimator.
#' @seealso \code{\link{dgaps}}: \eqn{D}-gaps estimator.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @seealso \code{\link{choose_b}}, \code{\link{choose_ud}} and
#'   \code{\link{choose_ud}} for choosing tuning parameters.
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
#' The matrix \code{cheeseboro} contains hourly maximum wind gusts (in miles
#' per hour) recorded at the Cheeseboro weather station near Thousand Oaks,
#' Southern California, USA during the month of January over the period
#' 2000-2009. These data are analysed in Reich and Shaby (2016).
#' @format A 744 by 10 numeric matrix.  Column \code{i} contains the hourly
#'   maximum wind gusts (in miles per hour) from Cheeseboro in the year
#'   2000 + \code{i} - 1. The columns are named 2000, 2001, ..., 2009 and the
#'   rows are named day\code{j}hour\code{k}, where \code{j} is the day of the
#'   month and \code{k} the hour of the day.
#' @note There are 42 missing values, located in 6 of the 10 years, namely
#'   2000-2003 and 2005-2006.
#' @source The Remote Automated Weather Stations USA Climate Archive at
#'   \url{https://raws.dri.edu/}, more specifically the Daily Summaries of the
#'   \href{https://raws.dri.edu/cgi-bin/rawMAIN.pl?caCCHB}{Cheeseboro page}.
#' @references Reich, B. J. and Shaby, B. A. (2016). 'Time series of Extremes',
#' in Dey, D. K. and Yan, J. (eds.) Extreme Value Modeling and Risk Analysis.
#' New York: Chapman and Hall/CRC, pp. 131-151.
"cheeseboro"

#' Uccle maximum daily temperatures
#'
#' The dataframe \code{uccle} contains daily maximum temperatures in degrees C
#' recorded at the Uccle, Belgium from 1/1/1833 to 23/1/2011.  The Station
#' identifier in the source file is 17 and the Source identifier is 117882.
#' @format A data frame with 65036 observations on the following and 5
#'   variables.
#' \itemize{
#' \item{\code{temp:}}{ daily maximum temperature in degrees C.}
#' \item{\code{year:}}{ the year.}
#' \item{\code{month:}}{ the month of the year.}
#' \item{\code{day:}}{ day of the month.}
#' \item{\code{date:}}{ date with the \code{\link[base:Dates]{Date}} class,
#'   in the format YYYY-MM-DD.}
#' }
#' @note There are 5336 missing values.
#' @source Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of
#'  20th-century surface air temperature and precipitation series for the
#'  European Climate Assessment. \emph{Int. J. of Climatol.}, \strong{22},
#'  1441-1453 \doi{10.1002/joc.773}. Data and metadata available at
#'  \href{https://www.ecad.eu}{https://www.ecad.eu}.  The data were downloaded
#'  on 27/3/2022 using a
#'  \href{https://www.ecad.eu/dailydata/customquery.php}{Custom query (ASCII)},
#'  selecting "non-blend" for type of series.
"uccle"

#' 20th century Uccle maximum daily temperatures in July - data frame
#'
#' The dataframe \code{uccle720} contains daily maximum temperatures in degrees C
#' recorded at the Uccle, Belgium during July for the years 1901 to 1999.
#' The Station identifier in the source file is 17 and the Source identifier is
#' 117882.  These data are analysed in Holesovsky and Fusek (2020).
#' @format A data frame with 3100 observations on the following and 5
#'   variables.
#' \itemize{
#' \item{\code{temp:}}{ daily maximum temperature in degrees C.}
#' \item{\code{year:}}{ the year.}
#' \item{\code{month:}}{ the month of the year.}
#' \item{\code{day:}}{ day of the month.}
#' \item{\code{date:}}{ date with the \code{\link[base:Dates]{Date}} class,
#'   in the format YYYY-MM-DD.}
#' }
#' @note There are 6 missing values, one located in each of the years
#'   1925, 1926, 1956, 1963, 1969 and 1976.
#' @source Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of
#'  20th-century surface air temperature and precipitation series for the
#'  European Climate Assessment. \emph{Int. J. of Climatol.}, \strong{22},
#'  1441-1453 \doi{10.1002/joc.773}. Data and metadata available at
#'  \href{https://www.ecad.eu}{https://www.ecad.eu}.  The data were downloaded
#'  on 27/3/2022 using a
#'  \href{https://www.ecad.eu/dailydata/customquery.php}{Custom query (ASCII)},
#'  selecting "non-blend" for type of series.
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. \emph{Extremes}, \strong{23}, 197-213
#'   (2020). \doi{10.1007/s10687-020-00374-3}
#' @examples
#' uccle720_ts <- ts(uccle720$temp, start = c(1901, 1), frequency = 31)
#' plot(uccle720_ts, ylab = "daily maximum temperature in July / degrees C",
#'      xlab = "year")
"uccle720"

#' 20th century Uccle maximum daily temperatures in July - matrix
#'
#' The matrix \code{uccle720m} contains daily maximum temperatures in degrees C
#' recorded at the Uccle, Belgium during July for the years 1901 to 1999.
#' The Station identifier in the source file is 17 and the Source identifier is
#' 117882.  These data are analysed in Holesovsky and Fusek (2020).
#' @format A 31 by 100 numeric matrix.  Column \code{i} contains the maximum
#'   daily temperature in degrees C at Uccle in the year 1900 + \code{i} - 1.
#'   The columns are named 1900, 1901, ..., 1999 and the rows are named
#'   after the day of the month: 1, 2, .., 31.
#' @note There are 6 missing values, one located in each of the years
#'   1925, 1926, 1956, 1963, 1969 and 1976.
#' @source Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of
#'  20th-century surface air temperature and precipitation series for the
#'  European Climate Assessment. \emph{Int. J. of Climatol.}, \strong{22},
#'  1441-1453 \doi{10.1002/joc.773}. Data and metadata available at
#'  \href{https://www.ecad.eu}{https://www.ecad.eu}.  The data were downloaded
#'  on 27/3/2022 using a
#'  \href{https://www.ecad.eu/dailydata/customquery.php}{Custom query (ASCII)},
#'  selecting "non-blend" for type of series.
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. \emph{Extremes}, \strong{23}, 197-213
#'   (2020). \doi{10.1007/s10687-020-00374-3}
#' @examples
#' uccle720_ts <- ts(uccle720$temp, start = c(1901, 1), frequency = 31)
#' plot(uccle720_ts, ylab = "daily maximum temperature in July / degrees C",
#'      xlab = "year")
"uccle720m"
