\name{Radiation}
\alias{Radiation}
\docType{data}
\title{Global radiation data}
\description{Rank-order association between monthly mean daily global radiation, mean sunshine duration and mean temperature.}
\usage{data(Radiation)}
\format{
  A data frame with 12 observations on the following variables.
  \describe{
  	\item{Mt}{Mean temperature, \eqn{^\circ}C.}
  	\item{Gr}{Global radiation, MJ\eqn{m^2}.}
  	\item{Msd}{Mean sunshine duration, h.}
  	}}
  
\details{Montly data.}
\examples{
	data(Radiation)
	head(Radiation)
}

\source{
Sivamadhavi, V. and Selvaraj, R. S. (2012). Robust regression technique to estimate the global radiation. \emph{Indian Journal of Radio & Space Physics}, \bold{41}, 17-25.
	}
\keyword{datasets}