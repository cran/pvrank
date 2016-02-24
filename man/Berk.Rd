\name{Berk}
\alias{Berk}
\docType{data}
\title{Effect of outliers}
\description{Data on the average number of births and deaths by the time of the day for a particular hospital in Brussels. The data cover a 30-year period in the nineteenth century.}
\usage{data(Berk)}
\format{
A data frame with 24 observations on the following two variables.
  \describe{
  	\item{Births}{Number of births by hour.}
  	\item{Deaths}{Number of deaths by hour.}
  	}
  }
  
\details{Twenty-two observations are clustered and show little association. Two observations (for noon and midnight) are dramatically smaller in both the y-direction and x-direction. With these two included, there is obviously a positive correlation in the data. However, a direct association between the two variables is doubtful, for if the outliers are removed then all correlations decrease and the associated p-values increase up to the point where the null hypothesis of independence cannot be rejected at any reasonable level.}

\examples{
	data(Berk)
	head(Berk)
}

\source{
Berk, R. A. (1990). "A primer on robust regression". In Fox, J. and Scott Long, J. (Eds), 292--324. \emph{Modern Methods of Data Analysis}. Sage Publications, Newbury Park, Ca, USA.
}

\keyword{datasets}
