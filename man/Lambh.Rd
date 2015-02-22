\name{Lambh}
\alias{Lambh}
\docType{data}
\title{Volcanic dust in the atmosphere}
\description{Importance to meteorology, principally the dust veils created in the atmosphere, particle sizes and distribution, heights, fall speeds and atmospheric residence times.}
\usage{data(Lambh)}
\format{
Two ime series data on the average monthly temperature for Central England and on the worldwide volcanic dust veil.
  \describe{
  	\item{DVI}{Worldwide volcanic dust veil.}
  	\item{Temp}{Average monthly temperature for Central England.}
  	}
  }
  
\details{The DVI is an attempt to place a quantitative measure on the magnitude of volcanic eruptions and allows comparisons of volcanic eruptions by a single numerical index.}

\examples{
	data(Lambh)
	head(Lambh)
}

\source{
Lamb, H. H. (1970). Volcanic dust in the atmosphere; with a chronological assessment of its meteorological significance. \emph{Philosophical Transactions of the Royal Society of London}, \bold{266}, 425-527.
}

\keyword{datasets}
