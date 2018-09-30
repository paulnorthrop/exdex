#' exdex: Estimation of the Extremal Index
#'
#' Performs frequentist inference for the extremal index using the methodologies
#' described in Suveges and Davison (2010)  \url{http://dx.doi.org/10.1214/09-AOAS292} and in
#' Northrop (2015) \url{http://dx.doi.org/10.1007/s10687-015-0221-5}.
#'
#' @details Add details.
#'
#'   See \code{vignette("exdex-vignette", package = "exdex")} for an
#'   overview of the package.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{http://dx.doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @docType package
#' @name exdex
#' @import methods
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
#'   272-283.  \url{https://doi.org/10.1002/env.2133}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index. \emph{Extremes}, \strong{18},
#'   585-603.  \url{https://doi.org/10.1007/s10687-015-0221-5}
"newlyn"

#' Worldwide Terrorism Incidents
#'
#' A database of information about terrorism incidents produced by the
#' \href{https://www.rand.org/}{RAND corporation} spanning the time period
#' 1968 - 2009.
#'
#' The dataframe \code{terrorism} contains 40129 rows and the following 8
#' columns:
#'   \itemize{
#'     \item{Column 1, \code{Date}: }{The date on which the incident took
#'       place in the format YYYY-MM-DD. Has class "Date". Some days have
#'       more than one incident, others (thankfully) have none.}
#'     \item{Column 2, \code{City}: }{The city in which the incident took
#'     place.}
#'     \item{Column 3, \code{Country}: }{The country in which the incident took
#'     place.}
#'     \item{Column 4, \code{Perpetrator}: }{The group responsible for the
#'     incident (if known).}
#'     \item{Column 5, \code{Weapon}: }{The weapon used (if known).}
#'     \item{Column 6, \code{Injuries}: }{The number of people injured.}
#'     \item{Column 7, \code{Fatalities}: }{The number of people killed.}
#'     \item{Column 8, \code{Description}: }{A brief description of the
#'     incident.}
#'  }
#' A subset of these data (1968-2007 only) are analysed in Clauset and
#' Woodard (2013).
#' @format A dataframe with 40129 rows an 8 columns.
#' @source The RAND Database of Worldwide Terrorism Incidents.
#'   \url{https://www.rand.org/nsrd/projects/terrorism-incidents/download.html}
#' @references Clauset, A. and Woodard, R. (2013) Estimating the historical and
#'   future probabilities of large terrorist events, \emph{The Annals of
#'   Applied Statistics}, \strong{7}(4), 1838-1865.
#'   \url{http://dx.doi.org/10.1214/12-AOAS614}
"terrorism"
