\name{qrank}
\alias{qrank}
\title{Calculation of quantiles of rank correlations}
\description{For a given value of the test statistic, this routine computes approximated or exact conservative and liberal critical values under the hypothesis of no association.}
\usage{qrank(prob, n, index="spearman", approx="vggfr", print = FALSE, lower.tail = TRUE)}
\arguments{
  \item{prob}{the nominal level of significance.} 
  \item{n}{the number of ranks.}
  \item{index}{a character string that specifies the rank correlation used in the test. Acceptable values are:
  "spearman","kendall","gini","r4". Only enough of the string to be unique is required.}
  \item{approx}{a character string that specifies the type of approximation to the null distribution: "vggfr", "exact","gaussian","student".} 
  \item{print}{\code{FALSE} suppresses partial output.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probability \eqn{P [X\leq x]}{ P[X <= x]} is computed, otherwise \eqn{P[X>x]}{P[X>x]}.}  
}
\details{
This routine provides two exact quantiles corresponding to a conservative level (next smaller exact size) and a liberal level (next larger exact size). However, in the case of \eqn{n>26} (Spearman) or \eqn{n>60} (Kendall) or \eqn{n>24} (Gini) or \eqn{n>15} (r4), approximated quantiles are provided according to \code{approx}. The default option is \code{"vggfr"}. 

Exact computations use frequencies obtained by complete enumeration for Spearman's rho, Gini's gamma, and r4. A recursive formula is employed in the case of Kendall's tau.
}
\value{a list containing the following components:
\item{n}{number of ranks.}
\item{Statistic}{coefficient of rank order association}
\item{Level}{nominal level}
\item{Cq}{conservative quantile}
\item{Cv}{conservative p-value}
\item{Lq}{liberal quantile}
\item{Lv}{liberal p-value}
}
\author{Agostino Tarsitano, Ilaria Lucrezia Amerise and Marco Marozzi}
\examples{
qrank(0.10,61,"K","St")
qrank(0.01,25,"S","Ga")
qrank(0.05,27,"G","Vg",TRUE,FALSE)
qrank(0.03,11,"r4","Ex",TRUE,FALSE)
}