\name{Timc}
\alias{Timc}
\title{FORTRAN subroutine for computing VGGFR parameters}
\description{Calls a fortran subroutine which performs a controlled random search to compute the parameters of the Vianelli Generalized Gaussian with finite range distribution (VGGFR).}
\usage{Timc(n, mu2n, mu4n,icoef)}
\arguments{
\item{n }{number of ranks}
\item{mu2n}{variance of \eqn{r_{h,n}}}
\item{mu4n}{kurtosis of \eqn{r_{h,n}}}
\item{icoef}{an integer indicating the coefficient to be interpolated: 1=Spearman, 2=Gini, 3=Kendall, 4=r4, 5=Fisher-Yates, 6=Filliben.}
}
\value{
  \item{Lam}{A numeric vector giving the estimates of the two parameters of the VGGFR}
   \item{Eval}{Final value of the loss function}
   }
\details{
The VGGFR is a flexible density to which we can resort in the case the number of ranks is larger than the threshold for which the exact null distribution of a rank correlation is known, but lower than the threshold for which the asymptotic Gaussian approximation becomes valid.
	
Although there are several methods for estimating the parameters of the GGFR, we follow the moment-matching method as applied by Karian and Dudewicz (2000) [Ch 2.6, p. 104] in similar estimations. The second and fourth centered moments of a GGFR distribution are:

\deqn{\mu_2(\lambda)=\frac{B(3\lambda_1^{-1},\lambda_2+1)}{B(\lambda_1^{-1},\lambda_2+1)} \qquad \mu_4(\lambda)=\frac{B(5\lambda_1^{-1},\lambda_2+1)}{B(\lambda_1^{-1},\lambda_2+1)}} 

where \eqn{B(x,y)} denotes the beta function between the positive real-valued numbers \eqn{x} and \eqn{y}. The second and fourth moments of the null distributions of the rank correlations are known polynomials in \eqn{n}. 

\deqn{\mu_2(r_1)=\frac{1}{n-1} }
\deqn{\mu_2(r_2)=\frac{2(2n+5)}{9n(n-1)} }
\deqn{\mu_2(r_3)=\left[\frac{2}{3(n-1)}\right]\left[\frac{n^2+2+k_n}{n^2-k_n}\right] }	
\deqn{\mu_4(r_1)=\frac{3(25n^3-38n^2-35n+72)} {25n(n+1)(n-1)^3} }
\deqn{\mu_4(r_2)=\frac{100n^4+328n^3-127n^2-997n-372} {1350\left[0.5n(n-1)\right]^3} }	
\deqn{\mu_4(r_3)=\frac{4[35n^7-(111-35k_n)n^6+(153+29k_n)n^5-(366-59k_n)n^4+304+11k_n)n^3]} {n^{k_n}(105-2k_n)(n+k_n)^3(\!n-k_n)^4(n-3+k_n)}+}
\deqn{\frac{-[(456-114k_n)n^2-(912-492k_n)n+(1248-933k_n)]} {n^{k_n}(105-2k_n)(n+k_n)^3(\!n-k_n)^4(n-3+k_n)} }

The approximation to the null distribution of \eqn{r_4} is based on an estimation of the second moment obtained through a regression strategy.  In particular,

\deqn{\mu_2 (r_4)\approx \frac{1.00762}{(n-1)} }
\deqn{\mu_4(r_4)\approx \frac{1.0949159471}{\sqrt{n-1}}+\frac{38.7820781157}{(n-1)^2}-\frac{208.8267798530}{(n-1)^3}+\frac{396.3338168921}{(n-1)^4} }

See Tarsitano and Amerise (2016). The second and fourth moments of \eqn{r_5} are given in Fieller and Pearson (1961):

\deqn{\mu_2(r_5)=\frac{1}{(n-1)}}
\deqn{\mu_4(r_5)=\frac{1}{(n-1)^2}\Big[\frac{3(n-1)}{n+1}+\frac{(n-2)(n-3)}{n(n^2-1)}\Big(\frac{k_4}{k_2^2}\Big)^2 \Big]}

with

\deqn{k_2=\frac{\sum_{i=1}^n\xi(x_i|n)^2}{n-1} }
\deqn{k_4=\frac{n\left\{(n+1)\sum_{i=1}^n\xi(x_i|n)^4-\frac{3(n-1)}{n}\Big[\sum_{i=1}^n\xi(x_i|n)^2\Big]^2\right\}}{(n-1)(n-2)(n-3)} }	

where \eqn{\xi(x_i|n)} is the expected values of the i-th largest standardized deviate in a sample of size \eqn{n} from a Gaussian population. With regard to \eqn{r_6}, the GGFR approximation is provisionally implemented with the same structure as \eqn{r_5} using the medians in place of the means of the Gaussian order statistics. See Amerise & Tarsitano (2016).

The estimate of the parameter vector \eqn{\lambda=(\lambda_1,\lambda_2)} is obtained by solving

\deqn{C(\lambda)=\min[\max[g_2(\lambda),g_4(\lambda)]]}
with \eqn{g_2(\lambda)=\mu_2(\lambda)-\mu_{2,n}} and \eqn{g_4(\lambda)=\mu_4(\lambda)-\mu_{4,n}} where 
\eqn{\mu_{2,n}} and \eqn{\mu_{4,n}} are the second and fourth moments of the given rank correlation.

}

\note{

The GGFR density has two exponential parameters that make GGFR highly nonlinear. In addition, the presence of a beta function depending on the unknown parameters can create difficulties in numerical stability. To increase the chances of getting a global solution in reasonable computational times, TIMC executes a controlled random search (CRS) algorithm discussed, for example, in Conlon (1992) and Brachetti et al. (1997).

The CRS is a direct search technique and is purely heuristic. Given a function of two or more variables, an initial search domain is defined by specifying limits to each variable. The algorithm starts by initially storing a predetermined number of points uniformly sampled over the search space. In our version, the size of the sample of initial random points is \eqn{30}. The points are obtained  by using a Halton-sequence (Haltpon, 1964) to get a pseudo-random distribution of points within the search domain. The sample  is then gradually contracted by replacing the current worst point in it with a better point. The stop criterion is defined in terms of the worst point  \eqn{C_M} and the best point \eqn{C_m} of the array. The iterations terminate if \eqn{C_M-C_m<10^{-16}}.

The parameter space in which the search proceeds is \eqn{1.5\le \lambda_1\le 3.5}, \eqn{0.3n\le\lambda_2\le 0.6n} (Spearman's \eqn{r_1}), \eqn{0.7n\le\lambda_2\le 1.2n} (Kendall's \eqn{r_2}) and \eqn{0.3n\le \lambda_2\le 0.9n} (Gini's \eqn{r_3}), with \eqn{kn=n\ mod\ 2}. These limits were chosen via a trial-and-error process. The maximum number of function evaluations allowed is set at \eqn{100000\sqrt{Log10(n)}}.
}

\references{
	
Amerise, I. L. and Tarsitano, A. (2016). A re-edition of the probability plot correlation coefficient leading to a new rank correlation. Submitted.
	
Brachetti, P. et al. (1997). "A new version of the Price's algorithm for global optimization". \emph{Journal of Global Optimization}, \bold{10}, 165-184.

Conlon, M. (1992).  "The controlled random search procedure for function optimization". \emph{Communications in Statistics - Simulation and Computation}, \bold{21}, 912--923.

Fieller, E. C. and Pearson, E. S. (1961). Tests for rank correlation coefficients: II. \emph{Biometrika}, \bold{48}, 29--40.

Halton, j. H.  (1964). Algorithm 247: Radical-inverse quasi- random point sequence. \emph{Communications of the ACM}, \bold{7}, 701--702.

Karian, Z. A. and Dudewicz, E. J. (2000). \emph{Fitting statistical distributions. The generalized lambda distribution and generalized bootstrap methods}. Boca Raton, FL: CRC Press.

Tarsitano, A. and Amerise, I. L. (2016). "Modelling of the null distribution of rank correlations". Submitted.

Vianelli, S. (1968). "Sulle curve normali di ordine r per intervalli finiti delle variabili statistiche". \emph{Annali della Facolta` di Economia e Commercio dell'Universita` di Palermo}, 2.

Vianelli, S. (1983). "The family of normal and lognormal distributions of order r".  \emph{Metron}, \bold{41}, 3-10.

}

\author{Agostino Tarsitano}
