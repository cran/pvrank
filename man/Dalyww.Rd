\name{Dalyww}
\alias{Dalyww}
\docType{data}
\title{Dark contrasts: the paradox of high rates of suicide in happy places}
\description{U.S. and international data provide evidence for a paradox: the happiest places tend to have the highest suicide rates.}
\usage{data(Dalyww)}
\format{
A data frame with 51 observations on the following four variables.
  \describe{
  	\item{SRU}{Suicide rate (unadjusted).}
  	\item{ASHR}{Adjusted suicide hazard ratio.}
  	\item{RCLS}{Raw coefficients on life satisfaction.}
  	\item{ACLS}{Adjusted coefficients on life satisfaction.}
  	}
  }
  
\details{Assessing the correlation across states between their suicide rankings and their life satisfaction rankings allows us to get a better sense of the correlation between the two while still retaining all observations, including the states that are apparent outliers.}

\examples{
	data(Dalyww)
	head(Dalyww)
}

\source{
Daly, M. C. et al. S. (2011). "Dark contrasts: the paradox of high rates of suicide in happy places". \emph{Journal of Economic Behavior & Organization},
\bold{80}, 435--442.
}

\keyword{datasets}
