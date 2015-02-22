\name{comprank}
\alias{comprank}
\title{Computes various rank correlation coefficients}
\description{Computes various rank order correlations.}
\usage{comprank(x, y = NULL, index = "spearman", tiex = "woodbury", sizer = 1000000, 
	repl = 1000, print = TRUE)}
\arguments{
\item{x}{a numeric vector, matrix or data frame.}
\item{y}{NULL (default) or a vector with compatible dimensions to \eqn{x}.}
\item{index}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman", "kendall", "gini","r4". Only enough of the string to be unique is required.} 
\item{tiex}{a character string that specifies the method for breaking ties. Possible options are: "woodbury", "gh","wgh", "midrank","dubois". Only enough of the string to be unique is required. This option is ignored when ties are absent. The options: \code{ "midrank"} and \code{"dubois"} are also ignored in the calculation of \eqn{r_4}.} 
\item{sizer}{number of replications for resolving ties by randomization. Default \eqn{10^6}.}
\item{repl}{number of sampled pairs of permutations needed to apply the weighted max-min method. Default value: \eqn{10^3}.}
\item{print}{\code{FALSE} suppresses some of the output.}
}
\details{
Consider a fixed set of \eqn{n} distinct items ordered according to the different degree in which they possess two common attributes represented by \eqn{x} and \eqn{y}. Let us suppose that each attribute consist of a host of intangibles that can be ranked but not measured and that the evaluations are expressed in terms of an ordinal scale of \eqn{n} ranks: \eqn{(x_i, y_i), i=1,2, ..., n}  where \eqn{x} is re-arranged in the natural order.
 
 Gideon & Hollister (1987) pointed out that reasonable coefficients of rank correlation need to possess certain properties and gave a list of postulates for nonparametric measures of dependence based on work of R\eqn{\acute{e}}nyi (1959) and Schweizer & Wolff (1981). On this specific issue, see also Scarisini, 1984, King, T.S. & Chinchilli, V.M. (2001), and Genest & Plante (2003). In short, the properties that any index of rank correlation should satisfy are the following
 \itemize{
		\item \eqn{r(\sigma,\pi)} is defined for any pair of permutations \eqn{\sigma,\pi}.
		\item Comparability. \eqn{-1 \le r(\sigma,\pi) \le 1} with \eqn{r(\sigma,\sigma) = r(\pi,\pi) = 1} and \eqn{r(\sigma,\sigma^*) = r(\pi,\pi^*) = -1} where \eqn{\sigma^* = n+1- \sigma} and \eqn{\pi^* = n+1- \pi}.
		\item Symmetry. \eqn{r(\sigma,\pi) = r(\pi,\sigma)}.
		\item Right-invariance. \eqn{r(\sigma[\theta],\pi[\theta]) = r(\sigma,\pi)} for all \eqn{\theta \in S_n}, where \eqn{S_n} is the set of all \eqn{n!} permutations of integers \eqn{1, 2,\cdots, n}.
		\item Antisymmetry under reversal. \eqn{r(\sigma,\pi^*) = r(\sigma^*,\pi) = -r(\sigma,\pi^*)}.
		\item Zero expected value under independence. \eqn{E_{\sigma,\pi \in S_n}[r(\pi,\sigma)]=0} where all rankings are equally probable with probability \eqn{1/n!}.
	}

The \pkg{pvrank} includes four admissible rank correlations (in the sense that they have the desirable properties described above). 
 
\deqn{Spearman: r_1=1-(6/(n^3-n)S_1;\quad  S_1=\sum_{i=1}^n (x_i-y_i)^2}
\deqn{Kendall: r_2=(2/(n^2-n))S_2;\quad S_2=\sum_{i=1}^{n-1}\sum_{j=i+1}^n sign(x_i-x_j)sign(y_i-y_j)}
\deqn{Gini: r_3=(2/(n^2-k_n))S_3;\quad  S_3=\sum_{i=1}^n (|n+1-x_i-y_i|-|x_i-y_i|)}
\deqn{  r_4=\Big[\sum_{i=1}^ng_i(x,y^*)\Big]\Big[\sum_{i=1}^ng_j(x^*,y)\Big]-\Big[\sum_{i=1}^ng_i(x^*,y^*)\Big]\Big[\sum_{i=1}^ng_j(x,y)\Big]/M_n}

where the quantity \eqn{k_n} is zero if \eqn{n} is even and one if \eqn{n} is odd. The function \eqn{g} is defined as \eqn{g_h[x(i),x(j)]=max(x[i]/y[j],y[j]/x[i])}. The quantity \eqn{M_n} is given by \deqn{M_n=\Big[k_n+2\sum_{i=1}^{n/2}(n+1-i)/i\Big]^2 -n^2}

In the absence of ties \eqn{r_1, r_2, r_3} can assume \eqn{(n^3-n)/6+1}, \eqn{(n^2-n)/2+1}, \eqn{(n^2-k_n)/2+1} distinct values, respectively. The coefficient \eqn{r_4} can assume a number of different values of the order \eqn{0.25n!} more or less uniformly spaced from each other. When \eqn{n>3}, \eqn{r_1} can be zero if, and only if, \eqn{n} is not of the form \eqn{n=4m+2} where \eqn{m} is a positive integer. The coefficient \eqn{r_2} is zero or even if, and only if, \eqn{n=4m} or \eqn{n=4m+1}. The coefficient \eqn{r_2} only takes on an odd value if \eqn{n} is not in that form. For \eqn{n>3}, zero is always a value of \eqn{r_3}.  If \eqn{n\ge5}, coefficient \eqn{r_4} fails to be zero for \eqn{n=7}.

Equal values are common when rank methods are applied to rounded data or data consisting solely of small integers. A popular technique for resolving ties in rank correlation is the mid-rank method: the mean of the rankings remain unaltered, but the variance is reduced and changes according to the number and location of ties. In passing, we note that no satisfactory way was found for their use with \eqn{r_4}.

Woodbury (1940) notes that the quadratic mean has been proposed to preserve the sum of squares of the ranks: \eqn{n(n+1)(2n+1)/6}. In this case the condition on the total variance is satisfied, but not that on the total mean, which turns out to be greater than \eqn{(n+1)/2}. The Woodbury method yields a corrected values of the rank correlation which could be regarded as the average of the values of the conventional coefficient obtained by all possible rankings of tied observations. In essence, the methods "gh" (see Gideon and Hollister, 1987) and "wgh" (see Gini, 1939) consider two special orderings with the objective to determine the pair of permutations that maximizes the direct association or positive correlation and the pair of permutations that maximizes the inverse association or, in absolute terms, the negative correlation. The method "gh" yields the simple average of the bounds; the method "wgh" yields a weighted average of the bounds. See Amerise and Tarsitano (2014).
}
\value{The value of the required rank correlation.}

\references{

Amerise, I. L. and Tarsitano, A. (2014). Correction methods for ties in rank correlations. Submitted.

Genest, C. and Plante, J.-F. (2003). On Blest's measure of rank correlation. \emph{The Canadian Journal of Statistics},  \bold{31}, 35-52.

Gideon, R. and Hollister, A. (1987). A rank correlation coefficient resistant to outliers. \emph{Journal of the American Statistical Association}, \bold{82}, 656-666.

Gini, C. (1939). Sulla determinazione dell'indice di cograduazione. \emph{Metron}, \bold{13},  41-48.

King, T.S. and  Chinchilli, V.M. (2001). Robust estimators of the concordance correlation coefficient. \emph{Journal of Biopharmaceutical Statistics}, \bold{ 11}, 83-105.

R\eqn{\acute{e}}nyi, A. (1959). On measures of dependence. \emph{Acta Mathematica Hungarica}, \bold{10}, 441-451.

Scarsini, M. (1984). On measures of dependence. \emph{ Stochastica}, \bold{8}, 201-218.

Schweizer, B. and Wolff, E.F. (1981). On nonparametric measures of dependence for random variables. \emph{Annals of Statistics}, \bold{9}, 879-885.

Tarsitano, A. and Lombardo, R. (2013). A coefficient of correlation based on ratios of ranks and anti-ranks. \emph{Jahrbucher fur Nationalokonomie und Statistik}, \bold{233}, 206-224.

Woodbury, M. A. (1940). Rank correlation when there are equal variates. \emph{The Annals of Mathematical Statistics} \bold{11}, 358-362.
}
		
	
\author{Ilaria Lucrezia Amerise}
\examples{
###
data(Zoutus);attach(Zoutus)
comprank(LogDose,LogTime,"spearman","woodbury")
comprank(LogDose,LogTime,"kendall","woodbury")
detach(Zoutus)
###
data(Franzen);attach(Franzen)
plot(MECISSP,PPP, main="Environmental Attitudes in International Comparison", xlab=
"Mean Environmental concern ISSP", ylab="Purchasing power parity", pch=19, cex=0.8)
op<-par(mfrow=c(1,1))
abline(h=mean(PPP),col="black",lty=2,lwd=1)
abline(v=mean(MECISSP),col="black",lty=2,lwd=1)
par(op)
comprank(MECISSP,PPP,"kendall","gh")
comprank(MECISSP,PPP,"gini","wgh")
detach(Franzen)
###
data(Viscoh);attach(Viscoh)
Viscoh<-as.matrix(Viscoh)
a<-comprank(Viscoh,"spearman","gh",print=FALSE);a
b<-comprank(Viscoh,"r4","gh",print=FALSE);b
detach(Viscoh)
###
data(Laudaher);attach(Laudaher)
comprank(Duration,Infiltration,"gini","midrank")
comprank(Duration,Infiltration,"gini","dubois")
comprank(Duration,Infiltration,"spearman","midrank")
comprank(Duration,Infiltration,"spearman","dubois")
detach(Laudaher)
###
data(Radiation);attach(Radiation)
Mt<-c(25.32,26.80,28.81,30.92,32.91,32.35,31.21,30.35,29.99,28.41,26.67,25.40)
Gr<-c(17.61,21.11,23.07,23.55,21.94,20.60,19.10,18.95,19.70,16.17,14.19,14.65)
Msd<-c(7.86,9.15,8.77,9.11,8.24,6.77,5.45,5.35,6.39,5.39,5.37,6.35)
Radiation<-as.matrix(Radiation)
r1<-comprank(Radiation,"spearman","midrank",print=FALSE);eigen(r1)$values
r2<-comprank(Radiation,"kendall","midrank",print=FALSE);eigen(r2)$values
r3<-comprank(Radiation,"gini","midrank",print=FALSE);eigen(r3)$values
r4<-comprank(Radiation,"r4","gh",print=FALSE);eigen(r4)$values
detach(Radiation)
###
# Correlation matrix
data(Marozzi);attach(Marozzi)
Marozzi<-as.matrix(Marozzi)
cor1<-comprank(Marozzi,"spearman","midrank",print=FALSE)
rownames(cor1)<-colnames(Marozzi);rownames(cor1)<-colnames(Marozzi)
print(round(cor1,3))
cor1 <-comprank(Marozzi,"kendall","midrank",print=FALSE)
print(round(cor1,3))
cor1<-comprank(Marozzi,"gini","midrank",print=FALSE)
print(round(cor1,3))
cor1<-comprank(Marozzi,"r4","wgh",print=FALSE)
print(round(cor1,3))
detach(Marozzi)
}