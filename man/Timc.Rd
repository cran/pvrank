\name{Timc}
\alias{Timc}
\title{FORTRAN subroutine for computing VGGFR parameters}
\description{Calls a fortran subroutine which performs a controlled random search to compute the parameters of the Vianelli Generalized Gaussian with finite range distribution (VGGFR).}
\usage{Timc(n, mu2n, mu4n)}
\arguments{
\item{n }{number of ranks}
\item{mu2n}{variance of \eqn{r_{h,n}}}
\item{mu4n}{kurtosis of \eqn{r_{h,n}}}
}
\value{
  \item{Lam}{A numeric vector giving the estimates of the two parameters of the VGGFR}
   \item{EVAL}{final value of the loss function}
   }
\details{The VGGFR is a flexible density to which we can resort in the case the number of ranks is larger than the threshold for which the exact null distribution is known, but lower than the threshold for which the asymptotic Gaussian approximation becomes valid.
	
Since VGGFR has two parameters, it seems natural to involve the second and fourth moments of the null distributions of the rank correlations, which are known polynomials in \eqn{n}. The estimator of \eqn{\lambda=(\lambda_1,\lambda_2)} is obtained by solving
\deqn{C(\lambda)=\min[\max[g_2(\lambda),g_4(\lambda)]]}
with \eqn{g_2(\lambda)=\mu_2(\lambda)-\mu_{2,n}} and \eqn{g_4(\lambda)=\mu_4(\lambda)-\mu_{4,n}} where 
\eqn{\mu_{2,n}} and \eqn{\mu_{4,n}} are the second and fourth moments of the given rank correlation. Variance and kurtosis of the VGGFR are
\deqn{\mu_2(\lambda)=B(3\lambda_1^{-1},\lambda_2+1)/B(\lambda_1^{-1},\lambda_2+1)}
\deqn{\mu_4(\lambda)=B(5\lambda_1^{-1},\lambda_2+1)/(\mu_2(\lambda)B(\lambda_1^{-1},\lambda_2+1))-3}
}

\note{The controlled random search procedure is carried out on the rectangle [1,4]x[0.2n,0.5n].}

\references{
	
Brachetti, P. et al. (1997). "A new version of the Price's algorithm for global optimization". \emph{Journal of Global Optimization}, \bold{10}, 165-184.

Tarsitano, A. and Amerise, I. L. (2013). "Approximation of the null distribution of rank correlations". Submitted.

Vianelli, S. (1983). "The family of normal and lognormal distributions of order r". 
\emph{Metron}, \bold{41}, 3-10.
}

\author{Agostino Tarsitano}
