\name{FrigateB}
\alias{FrigateB}
\docType{data}
\title{Magnificent frigatebird}
\description{Males of the magnificent frigatebird (Fregata magnificens) have a large red throat pouch. They visually display this pouch and use it to make a drumming sound when seeking mates.}
\usage{data(FrigateB)}
\format{
A data frame with 18 observations on the following two variables.
  \describe{
  	\item{Vol}{Volume (\eqn{cm^3)}.}
  	\item{Frq}{Frequency (\eqn{Hz}).}
  	}
  }
  
\details{Madsen et al. (2004) wanted to know whether females, who presumably choose mates based on their pouch size, could use the pitch of the drumming sound as an indicator of pouch size. The authors estimated the volume of the pouch and the fundamental frequency of the drumming sound in 18 males.}

\examples{
	data(FrigateB)
	head(FrigateB)
}

\source{
Madsen, V. et al. (2004). "Bimodal signaling of a sexually selected trait: gular pouch drumming in the magnificent frigatebird". \emph{Condor}, \bold{106}, 156--160.
}

\keyword{datasets}
