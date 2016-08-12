\name{VGGFR}
\alias{VGGFR}
\title{Vianelli density}
\description{Plots Vianelly (generalized Gaussian) density with finite range in [-1,1].}
\usage{
VGGFR(L1, L2, add=FALSE, lwd=2, lty=5, col="blue", ylim=c(0,1), Main="", np=201)
}

\arguments{
  \item{L1}{positive shape parameter.}
  \item{L2}{positive shape parameter. Impacts more on the tails.}
  \item{add}{when add=TRUE the plot is superimposed to an existing graph.}
  \item{lwd}{weight of the line.}
  \item{lty}{the type of the line.}
  \item{col}{color of the curve.}
   \item{ylim}{numeric vectors of length 2, giving the \eqn{y} coordinates ranges.}
   \item{Main}{a string describing the graph.}
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
\note{If you want to use the output from VGGFR save the result and then select from the list the value(s) 
	needed.}

\references{
	
Tarsitano, A. and Amerise, I. L. (2016). "Modelling of the null distribution of rank correlations". Submitted.

Vianelli, S. (1968). "Sulle curve normali di ordine $r$ per intervalli finiti delle variabili statistiche". \emph{Annali della Facolt\`a di Economia e Commercio dell'Universit\`a di Palermo}, 2.
	
Vianelli, S. (1983). "The family of normal and lognormal distributions of order r". \emph{Metron}, \bold{41}, 3-10.}

\author{Agostino Tarsitano and Ilaria Lucrezia Amerise}

\examples{
	
# Density curve of a VGGFR model
VGGFR(2, 12, add=FALSE, lwd=2, lty=5, col="darkgreen", ylim=c(0,2), Main="", np=201)
#####
#
a<-ranktes(0.5, 28, "r4", "vg",FALSE, "two", FALSE)
b<-VGGFR(a$Lambda, add = FALSE, lwd = 2, lty = 5, col = "blue", ylim=c(0,2.5),np = 201)
#####
#
# Lambert's semicircular distribution of errors (1760,1765).
# Given a probability distribution, the value with the higher probability density is 
# deemed to be more probable than the value with the lower probability density.
#
VGGFR(2,0.5,col="red",ylim=c(0,0.75),Main="Lambert's distribution of errors")
#
#####
# Pearson type II used as an approximation to the null distribution of the Fisher-Yates
# rank correlation. Fieller, E. C. and Pearson, E. S. (1961). Tests for rank correlation 
# coefficients: II. Biometrika, 48, 29-40.
n<-10
VGGFR(2, (n-4)/2, add=FALSE, lwd=2, lty=5, col="magenta2", ylim=c(0,1.1), Main="", np=201)
abline(h=0);abline(v=0,lty=2,lwd=2,col="pink2")
#####
# 
# Save and use the results 
res<-VGGFR(1.5,5.5,add = FALSE, lwd = 2, lty = 1, col = "blue", ylim=c(0,2.5),np = 201)
res$kurt-res$oam/res$st.dev
#####
# 
# A family of symmetrical beta densities
VGGFR(2,1,col="black",ylim=c(0,1.4),Main="Symmetrical beta densities")
La<-seq(1,6,0.5)
for (L1 in La){VGGFR(2,L1,add=TRUE, lwd = 1, lty = 1, col=gray(L1/6))}
#####
# 
# A family of GGFR curves
VGGFR(1,2, lwd = 1, lty = 1,col="black",ylim=c(0,5))
La<-seq(1,6,0.5);Lg<-seq(0,1,1/12)
for (L1 in La){
		c2<-gray(Lg, alpha= 2/6)
		for (L2 in seq(1,12,1)){
			VGGFR(L1,L2,add=TRUE, lwd = 1, lty = 1, col=c2[L2])
			}}	
}