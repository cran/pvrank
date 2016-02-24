\name{Security}
\alias{Security}
\docType{data}
\title{Accessibility, security, and accuracy in statistical databases}
\description{Effectiveness of security mechanisms for protecting numerical databases.}
\usage{data(Security)}
\format{
A data frame with 40 observations on the following three variables.
  \describe{
  	\item{OV}{Original values.}
    \item{AP}{Additively perturbed.}
  	\item{MP}{Multiplicatively perturbed.}
  	}
  }

\details{The data set include a "true" series of values drawn from a lognormal population, sorted in ascending order. Data are perturbed by using the MP and the additive AP methods together with the relative position of the perturbed value.}

\examples{
	data(Security)
	head(Security)
}

\source{
Muralidhar, K. and Batra, D. and Kirs, P. J. (1995). "Security, and accuracy in statistical databases: the case for the multiplicative fixed data perturbation approach". \emph{Management Science}, \bold{41}, 1549--1564.
}

\keyword{datasets}
