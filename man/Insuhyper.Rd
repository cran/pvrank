\name{Insuhyper}
\alias{Insuhyper}
\docType{data}
\title{Role of insulin in endogenous hypertriglyceridemia}
\description{Dietary carbohydrate accentuation of endogenous triglyceride production has been studied in 33 patients. No significant correlation existed between fasting plasma triglyceride concentration and either plasma glucose or free fatty acid concentrations after the high carbohydrate diet, nor was the degree of hypertriglyceridemia related to degree of obesity.}
\usage{data(Insuhyper)}
\format{
In this dataset are described plasma glucose, insulin, FFA, and triglyceride responses to a high carbohydrate diet in 33 patients. The data frame includes 18 variables.
 \describe{
  	\item{Age}{years of age.}
  	\item{Sex}{M=male, F=female}
  	\item{Height}{cm}
  	\item{Weight}{kg}
  	\item{PI}{Ponderal Index, an estimate of obesity= Height(inches)/[weight(pounts)]^(1/3)}
  	\item{TG}{Mean plasma triglyceride concentration during the last 2 wk of the high carbohydrate dietary period}
  	\item{PG8am}{Plasma Glusose at 8 a.m. mg/100 ml}
  	\item{PG8am}{Plasma Glusose at 9 a.m. mg/100 ml}
  	\item{PG8am}{Plasma Glusose at 10 a.m.mg/100 ml}
  	\item{PG8am}{Plasma Glusose average mg/100 ml}
  	\item{PI8am}{Plasma Insulin at 8 a.m. \eqn{\mu}U/ml}
  	\item{PI8am}{Plasma Insulin at 9 a.m. \eqn{\mu}U/ml}
  	\item{PI8am}{Plasma Insulin at 10 a.m. \eqn{\mu}U/ml}
  	\item{PI8am}{Plasma Insulin average \eqn{\mu}U/ml}
  	\item{PF8am}{Plasma FFA at 8 a.m. \eqn{\mu}Eq/liter. }
  	\item{PF8am}{Plasma FFA at 9 a.m. \eqn{\mu}Eq/liter}
  	\item{PF8am}{Plasma FFA at 10 a.m. \eqn{\mu}Eq/liter}
  	\item{PF8am}{Plasma FFA average \eqn{\mu}Eq/liter}
  	}
  Row names indicate name initials of the patients. The last four variables are omitted for 13 patients.
  }
  
\details{Authors have reported that the magnitude of plasma triglyceride response in man to diets either high or low in carbohydrate was highly correlated with both the plasma glucose and insulin responses elicited by that diet. The results confirm their previous finding of a high correlation between insulin and triglyceride responses. However, with inclusion of hyperglycemic patients there was no correlation found between plasma glucose response and the magnitude of the triglyceride elevation produced by the high carbohydrate diet. Finally, there was no correlation between hypertriglyceridemia and degree of obesity, or magnitude of the FFA response.
}

\examples{
	data(Insuhyper)
	head(Insuhyper)
}

\source{
Reaven, G. M. et al.  (1967). "Role of insulin in endogenous hypertriglyceridemia". \emph{The Journal of Clinical Investigation}, \bold{46}, 1756--1767.
}

\keyword{datasets}
