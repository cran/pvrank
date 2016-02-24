\name{BlaAlt}
\alias{BlaAlt}
\docType{data}
\title{Pairs of measurements of forced expiratory volume}
\description{ A problem with the use of the correlation coefficient between the first and second measurements is that there is no reason to suppose that their order is important. If the order were important the measurements would not be repeated observations of the same thing. We could reverse the order of any of the pairs and get a slightly different value of the correlation coefficient between repeated measurements.}
	
\usage{data(BlaAlt)}
\format{
A data frame with 20 observations on the following two variables.
  \describe{
  	\item{Fev1}{First FEV (litres).}
  	\item{Fev2}{Second Fev (litres).}
  	}
  }
  
\details{Pairs of measurements of forced expiratory volume (litres) a few weeks apart from 20 Scottish school children.}

\examples{
	data(BlaAlt)
	head(BlaAlt)
}

\source{
Bland, M. J. and Altman, D. G. (1996). "Measurement error and correlation coefficients". \emph{BMJ: British Medical Journal}, \bold{313}, 41--42.
}

\keyword{datasets}
