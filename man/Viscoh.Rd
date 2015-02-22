\name{Viscoh}
\alias{Viscoh}
\docType{data}
\title{Data on death measures for various conditions}
\description{The study develops measures of mortality risks that recognize the probability of death, the duration of life lost, and the role of discounting.}
\usage{data(Viscoh)}
\format{
A data frame with 29 observations on the following 3 variables.
  \describe{
  	\item{Deaths}{Deaths per 100,000.}
  	\item{LLE}{Lost life expectancy.}
  	\item{EYLL}{Expected years of life lost.}
  	}
Lost Life Expectancies and Expected Years of Life Lost calculations utilize a \eqn{3\%} annual discount rate and 10-year lags for non-immediate causes of death. 	
  }
\details{EYLL reflects the combined influence of the probability of death and the amount of lost life expectancy. These figures are strongly correlated with the probability of death rankings. Most of the shifts in the rankings tend to be relatively minor as, for example, diabetes ranks 10th in terms of the probability of death and 15th in terms of the expected years of life lost.}
\examples{
	data(Viscoh)
	head(Viscoh)
}

\source{
Viscusi, W. K. and Cogan, J. F. and Hakes, J.K. (1997). Measures of mortality risks.  \emph{Journal of Risk and Uncertainty}, \bold{14}, 213-233.
}

\keyword{datasets}
