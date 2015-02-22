\name{Zoutus}
\alias{Zoutus}
\docType{data}
\title{Data on time and dose of CT fluoroscopy-guided procedures}
\description{Two correlation/association coefficients: Pearson \eqn{r_0} and Spearman \eqn{r_1}, are reviewed and compared.}
\usage{data(Zoutus)}
\format{
  A data frame with 19 observations on the following two variables.
  \describe{
  	\item{LogTime}{Total procedure time. Minutes on a natural log scale.}
  	\item{LogDose}{Dose of CT fluoroscopy-guided procedure. Rads on a natural log scale.}
  	}}
\details{Both correlation coefficients confirm that the log of total time and the log of dose are correlated strongly and positively.}
\examples{
	data(Zoutus)
	head(Zoutus)
}

\source{
Zou, K. H. and Tuncali, K. and Silverman, S. G. (2003). Correlation and simple linear regression. \emph{Radiology}, \bold{227}, 617-628.
}

\keyword{datasets}
