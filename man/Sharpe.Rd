\name{Sharpe}
\alias{Sharpe}
\docType{data}
\title{Performance of mutual funds, 1954-1963}
\description{Determine the extent to which differences in performance continue through time.}
\usage{data(Sharpe)}
\format{
A data frame with 34 observations on the following 2 variables.
  \describe{
  	\item{AVR}{Average annual return (per cent).}
  	\item{VAR}{Variability of annual return (per cent).}
  	}
  }
  
\details{The annual rates of return for thirty-four open-end mutual funds during the period 1954-63 are used to evaluate the predictive ability of indicators of mutual fund performance.}

\examples{
	data(Sharpe)
	head(Sharpe)
}

\source{
Sharpe, W. F. (1966). Mutual fund performance. \emph{The Journal of Business}, \bold{39}, 119-138.
}

\keyword{datasets}
