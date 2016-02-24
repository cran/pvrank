\name{Franzen}
\alias{Franzen}
\docType{data}
\title{Environmental attitudes data}
\description{This data set is used to compare the concern for the natural environment between the citizens of 26 countries that participated in the 2000 International (ISSP) surveys 1993 and 2000 and to explain the differences.}
\usage{data(Franzen)}
\format{
  A data frame with 25 observations (one is missing) on the following two variables.
  \describe{
  	\item{MECISSP}{Mean Environmental concern ISSP.}
  	\item{PPP}{Purchasing power parity.}
  	}}
  
\details{Northern Ireland data are omitted because of a missing value.}
\examples{
	data(Franzen)
	head(Franzen)
}

\source{
 Franzen, A. (2003). "Environmental attitudes in international comparison: an analysis of the ISSP Surveys 1993 and 2000". \emph{Social Science Quarterly}, \bold{84}, 297--308.
	}

\keyword{datasets}
