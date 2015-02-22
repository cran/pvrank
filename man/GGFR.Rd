\name{VGGFR}
\alias{VGGFR}
\title{Vianelli density}
\description{Plots Vianelly (generalized Gaussian) density with finite range in [-1,1].}
\usage{
VGGFR(L1, L2, add = FALSE, lwd = 2, lty = 5, col = "blue", np = 201)
}
\arguments{
  \item{L1}{positive shape parameter.}
  \item{L2}{positive shape parameter. Impacts more on the tails.}
  \item{add}{when add=TRUE the plot is superimposed to an existing graph.}
  \item{lwd}{weight of the line.}
  \item{lty}{the type of the line.}
  \item{col}{color of the curve.}
  \item{np}{number of points to be plotted.}
}
\details{The VGGFR density is given by
	\deqn{f(r;\lambda_1,\lambda_2)=\lambda_1(1-|r|^{\lambda_2})^{\lambda_1}/[2B(1/\lambda_1,\lambda_2+1)]}
	where \eqn{\lambda_1,\lambda_2>0} and \eqn{B()} is the beta function.
	}
\value{The value returned is a list contaning:
	\item{st.dev}{standard deviation}
    \item{kurt}{kurtosis}
    \item{oam}{ordinate at the mode}
    }
\note{If you want to use the output from VGGFR save the result and then select from the list the value(s) needed.}

\references{
	
	Tarsitano, A. and Amerise, I. L. (2013). Approximation of the null distribution of rank correlations. Submitted.
	
	Vianelli, S. (1983). The family of normal and lognormal distributions of order r. \emph{Metron}, \bold{41}, 3-10.}

\author{Agostino Tarsitano and Ilaria Lucrezia Amerise}

\examples{
	
# Density curve of a VGGFR model
VGGFR(2,12)
#
a<-ranktes(0.5, 28, "r4", "vg",FALSE, "two", FALSE)
b<-VGGFR(a$Lambda, add = FALSE, lwd = 2, lty = 5, col = "blue", np = 201)
#
# A family of density curves
VGGFR(1,2,col="black")
La<-seq(1,6,0.5);Lg<-seq(0,1,1/12)
for (L1 in La){
	c2<-gray(Lg, alpha= 2/6)
	for (L2 in seq(1,12,1)){
	VGGFR(L1,L2,add=TRUE,col=c2[L2])}}
# Save and use the results 
res<-VGGFR(1.5,5.5)
res$kurt-res$oam/res$st.dev
}