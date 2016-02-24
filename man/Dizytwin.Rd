\name{Dizytwin}
\alias{Dizytwin}
\docType{data}
\title{Latitude and dizygotic twinning rates}
\description{Age-standardized dizygotic twinning rates correlate positively and significantly with latitude.}
\usage{data(Dizytwin)}
\format{
A data frame with 19 observations on the following two variables.
  \describe{
  	\item{Latitude}{Geographical latitude in Europe.}
  	\item{DZT_Rate}{Dizygotic twinning rates across the countries.}
  	}
  }
  
\details{In Europe, both mean birth weight and age-standardized dizygotic twinning rates correlate positively and significantly with latitude. It seems possible that these two sorts of correlation have explanations in common via pituitary action.}

\examples{
	data(Dizytwin)
	head(Dizytwin)
}

\source{
James, W. H. (1985). "Dizygotic twinning, birth weight and latitude". \emph{Annals of Human Biology},\bold{12}, 441--447.
}

\keyword{datasets}
