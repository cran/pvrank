\name{comprank}
\alias{comprank}
\title{Computes various rank correlation coefficients}
\description{Computes various rank order correlations with and withou ties present.}
\usage{comprank(p, q = NULL, indr, tiex = "woodbury", sizer = 100000,
	 repgin = 1000, print = TRUE)}
\arguments{
\item{p}{a numeric vector, matrix or data frame.}
\item{q}{NULL (default) or a vector with compatible dimensions to \eqn{p}.}
\item{indr}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman", "kendall", "gini","r4", "fy1" (Fisher-Yates - means), "fy2" (Fisher-Yates - medians), "sbz" (Symmetrical Borroni-Zenga). Only enough of the string to be unique is required.}
\item{tiex}{a character string that specifies the method for breaking ties. Possible options are: "woodbury", "gh","wgh", "midrank","dubois". Only enough of the string to be unique is required. This option is ignored when ties are absent. The options: \code{"midrank"} and \code{"dubois"} are also ignored in the calculation of \eqn{r_4}, Fisher-Yates means and Fisher-Yates medians.}
\item{sizer}{number of replications for resolving ties by randomization. Default value: \eqn{10^5}.}
\item{repgin}{number of sampled pairs of permutations needed to apply the weighted max-min method (wgh). Default value: \eqn{10^3}.}
\item{print}{if \code{FALSE} suppresses some of the output.}
}
\details{
Consider a fixed set of \eqn{n} distinct items ordered according to the different degree in which they possess two common attributes represented by \eqn{p} and \eqn{q}. Let us suppose that each attribute consist of a host of intangibles that can be ranked but not measured and that the evaluations are expressed in terms of an ordinal scale of \eqn{n} ranks: \eqn{(p_i, q_i), i=1,2, ..., n}  where \eqn{p} is re-arranged in the natural order.

 Gideon & Hollister (1987) pointed out that reasonable coefficients of rank correlation need to possess certain properties and gave a list of postulates for nonparametric measures of dependence based on works of R\eqn{\acute{e}}nyi (1959) and Schweizer & Wolff (1981). On this specific issue, see also Scarisini, (1984), King & Chinchilli (2001), and Genest & Plante (2003). In short, the properties that any index of rank correlation should satisfy are the following
 \itemize{
		\item \eqn{r(p,q)} is defined for any pair of permutations \eqn{p,q}.
		\item Comparability. \eqn{-1 \le r(p,q) \le 1} with \eqn{r(p,p) = r(q,q) = 1} and \eqn{r(p,p^*) = r(q,q^*) = -1} where \eqn{p^* = n+1- p} and \eqn{q^* = n+1- q}. It is also
required that increasing levels (decreasing levels) of agreement are reflected by increasing positive (negative) values.
		\item Symmetry. \eqn{r(p,q) = r(q,p)}.
		\item Right-invariance. \eqn{r(p[\theta],q[\theta]) = r(p,q)} for all \eqn{\theta, p, q \in S_n}, where \eqn{S_n} is the set of all \eqn{n!} permutations of integers \eqn{1, 2,\cdots, n}.
		\item Antisymmetry under reversal. \eqn{r(p,q^*) = r(p^*,q) = -r(p,q^*)}.
		\item Zero expected value under independence. \eqn{E_{p,q \in S_n}[r(q,p)]=0} when all rankings are equally probable with probability \eqn{1/n!}
	}

 \pkg{pvrank} includes six admissible rank correlations (in the sense that they have the desirable properties described above).

\deqn{Spearman: r_1=1-(6/(n^3-n)S_1;\quad  S_1=\sum_{i=1}^n (p_i-q_i)^2}
\deqn{Kendall: r_2=(2/(n^2-n))S_2;\quad S_2=\sum_{i=1}^{n-1}\sum_{j=i+1}^n sign(p_i-q_j)sign(p_i-q_j)}
\deqn{Gini: r_3=(2/(n^2-k_n))S_3;\quad  S_3=\sum_{i=1}^n (|n+1-p_i-q_i|-|p_i-q_i|)}
\deqn{  r_4=\Big[\sum_{i=1}^ng_i(p,q^*)\Big]\Big[\sum_{i=1}^ng_j(p^*,q)\Big]-\Big[\sum_{i=1}^ng_i(p^*,q^*)\Big]\Big[\sum_{i=1}^ng_j(p,q)\Big]/M_n}
where \eqn{k_n} is zero if \eqn{n} is even and one if \eqn{n} is odd. The function \eqn{g\left(.\right)} is defined as \eqn{g_h(p_i,q_j)=max(p_i/q_j,q_j/p_i)}. The quantity \eqn{M_n} is given by
\deqn{M_n=\Big[k_n+2\sum_{i=1}^{n/2}(n+1-i)/i\Big]^2 -n^2}

Three other admissible coefficients are:

\deqn{FY - means: r_5=\frac{\sum_{i=1}^n \xi(p_i|n)\xi(q_i|n)}{\sum_{i=1}^n \xi(p_i|n)^2}}
\deqn{FY - medians: r_6=\frac{\sum_{i=1}^n \zeta(p_i|n)\zeta(q_i|n)}{\sum_{i=1}^n \zeta(p_i|n)^2}}
\deqn{Symmetrical BZ: r_7=1.5\frac{\sum_{i=1}^n\sum_{j=1}^n|(i+p_i)-(j+p_j|-|(i-p_i)-|(j-p_j)|)}{n(n^2-1)}}

The expression \eqn{\xi(p_i|n)} denotes the expected values of the i-th largest standardized deviate in a sample of size \eqn{n} from a standad Gaussian population. See, for example, function evNormOrdStats in package \code{EnvStats}, Millard (2013). An implementation of the Gaussian rank correlation coefficient is given in the \code{rococo} package. See Bodenhofer et al. (2013). Finally, \eqn{\zeta(p_i|n)} is the i-th order statistic median from a standard Gaussian distribution. These medians \eqn{\zeta(p_i|n)} are exactly related to the order statistic medians \eqn{m_i} from a uniform distribution on \eqn{[0,1]}, that is,  \eqn{\zeta(p_i|n)=\Phi^{-1}(m_i)} where \eqn{\Phi^{-1}(x)} is the quantile function of the standard Gaussian cumulative distribution function. The medians of the order statistics from a unit uniform distribution are given by \eqn{m_i=B^{-1}(0.5,i,n+1-i)} where \eqn{B^{-1}(x,a,b)} denotes the inverse Beta distribution with parametes \eqn{a,b} evaluated at the percentage \eqn{x}. Fortunately, \eqn{B^{-1}(,)} is a standard function: \code{qbeta} in package \code{stats}. See Tarsitano & Amerise (2016) and Amerise & Tarsitano (2016).

	In the absence of ties \eqn{r_1, r_2, r_3} can assume \eqn{(n^3-n)/6+1}, \eqn{(n^2-n)/2+1}, \eqn{(n^2-k_n)/2+1} distinct values, respectively. The coefficient \eqn{r_4} can assume a number of different values of the order \eqn{0.25n!} more or less uniformly spaced from each other. When \eqn{n>3}, \eqn{r_1} can be zero if, and only if, \eqn{n} is not of the form \eqn{n=4m+2} where \eqn{m} is a positive integer. The quantity \eqn{S_2} is even if, and only if, \eqn{n=4m} or \eqn{n=4m+1}; \eqn{S_2} only takes on an odd value if \eqn{n} is not in that form. For \eqn{n>3}, zero is always a value of \eqn{r_3}.  If \eqn{n\ge5}, coefficient \eqn{r_4} fails to be zero for \eqn{n=7}. The number of distinct values of \eqn{r_5} and \eqn{r_6} is only limited by the number of bits used in the representation of float-point values in a computer.

	Equal values are common when rank methods are applied to rounded data or data consisting solely of small integers. A popular technique for resolving ties in rank correlation is the mid-rank method: the mean of the rankings remain unaltered, but the variance is reduced and changes according to the number and location of ties.

	Woodbury (1940) notes that the quadratic mean has been proposed to preserve the sum of squares of the ranks: \eqn{n(n+1)(2n+1)/6}. In this case the condition on the total variance is satisfied, but not that on the total mean, which turns out to be greater than \eqn{(n+1)/2}. The Woodbury method yields a corrected values of the rank correlation which could be regarded as the average of the values of the conventional coefficient obtained by all possible rankings of tied observations. In essence, the methods "gh" (see Gideon and Hollister, 1987) and "wgh" (see Gini, 1939) consider two special orderings with the objective to determine the pair of permutations that maximizes the direct association or positive correlation and the pair of permutations that maximizes the inverse association or, in absolute terms, the negative correlation. The method "gh" yields the simple average of the bounds; the method "wgh" yields a weighted average of the bounds. See Amerise and Tarsitano (2014).
}
\value{ a list with the following components :

\item{r}{The value of the required rank correlation.}

\item{ities}{Method for dealing with ties.}	

}

\references{

Amerise, I. L. and Tarsitano, A. (2014). Correction methods for ties in rank correlations. \emph{Journal of Applied statistics}, \bold{42}, 1--13.

Amerise, I. L. and Tarsitano, A. (2016). A new symmetrical test of bivariate independence. Submitted.
  
Borroni, G. C. (2013). "A new rank correlation measure". \emph{Statistical Papers}, \bold{54}, 255--270.

Bodenhofer, U. et al. (2013). "Testing noisy numerical data for monotonic association". \emph{information Sciences}, \bold{245}, 21--37.

Genest, C. and Plante, J.-F. (2003). "On Blest's measure of rank correlation". \emph{The Canadian Journal of Statistics},  \bold{31}, 35--52.

Gideon, R. and Hollister, A. (1987). "A rank correlation coefficient resistant to outliers". \emph{Journal of the American Statistical Association}, \bold{82}, 656--666.

Gini, C. (1939). "Sulla determinazione dell'indice di cograduazione". \emph{Metron}, \bold{13},  41--48.

King, T.S. and  Chinchilli, V.M. (2001). "Robust estimators of the concordance correlation coefficient". \emph{Journal of Biopharmaceutical Statistics}, \bold{11}, 83--105.

Millard, S.P. (2013). EnvStats: An R Package for Environmental Statistics. Springer, New York.

R\eqn{\acute{e}}nyi, A. (1959). "On measures of dependence". \emph{Acta Mathematica Hungarica}, \bold{10}, 441--451.

Scarsini, M. (1984). "On measures of dependence". \emph{ Stochastica}, \bold{8}, 201--218.

Schweizer, B. and Wolff, E.F. (1981). "On nonparametric measures of dependence for random variables". \emph{Annals of Statistics}, \bold{9}, 879--885.

Tarsitano, A. and Lombardo, R. (2013). "A coefficient of correlation based on ratios of ranks and anti-ranks". \emph{Jahrbucher fur Nationalokonomie und Statistik}, \bold{233}, 206--224.

Tarsitano, A. and Amerise, I. L. (2016). "Effectiveness of rank correlation statistics in non-linear relationships". Submitted.

Woodbury, M. A. (1940). "Rank correlation when there are equal variates". \emph{The Annals of Mathematical Statistics} \bold{11}, 358--362.

}

\author{Ilaria Lucrezia Amerise}

\examples{

data(Zoutus);attach(Zoutus);print(Zoutus)
a<-comprank(LogDose,LogTime,"spearman","woodbury")
cat(a$r,a$ities,"\n")
a<-comprank(LogDose,LogTime,"kendall","woodbury")
cat(a$r,a$ities,"\n")
detach(Zoutus)
#####
#
# Yotopoulos, P. A.  Nugent, J. B. (1973). A balanced-growth version of the
# linkage hypothesis: a test. The Quarterly Journal of Economics, 87, 157-171.
x<-1:18
y<-c(7,1,4,8,3,9,2,5,10,6,17,13,14,12,11,16,15,18)
a<-comprank(x,y,"gini");cat(a$r,a$ities,"\n")
#####
#
data(Franzen);attach(Franzen)
op<-par(mfrow=c(1,1))
plot(MECISSP,PPP, main="Environmental Attitudes in International Comparison",
xlab="Mean Environmental concern ISSP", ylab="Purchasing power parity",
pch=19, cex=0.8,col="salmon3")
abline(h=mean(PPP),col="darkred",lty=2,lwd=1)
abline(v=mean(MECISSP),col="darkred",lty=2,lwd=1)
par(op)
a<-comprank(MECISSP,PPP,"kendall","gh");cat(a$r,a$ities,"\n")
a<-comprank(MECISSP,PPP,"gini","wgh");cat(a$r,a$ities,"\n")
detach(Franzen)
#####
#
data(Viscoh);attach(Viscoh)
Viscoh<-as.matrix(Viscoh)
a<-comprank(Viscoh,"spearman","gh",print=FALSE)
print(a$r);cat(" method:", a$ities,"\n")
b<-comprank(Viscoh,"r4","gh",print=FALSE)
print(b$r);cat(" method:", b$ities,"\n")
c<-comprank(Viscoh,"fy1","wgh",print=FALSE)
print(c$r);cat("method:", c$ities,"\n")
d<-comprank(Viscoh,"fy2","wgh",print=FALSE)
print(d$r);cat(" method:", d$ities,"\n")
d<-comprank(Viscoh,"sbz","wgh",print=FALSE)
print(d$r);cat(" method:", d$ities,"\n")
detach(Viscoh)
#####
#
data(Laudaher);attach(Laudaher)
a1<-comprank(Duration,Infiltration,"gini","midrank")
a2<-comprank(Duration,Infiltration,"gini","dubois")
a3<-comprank(Duration,Infiltration,"spearman","midrank")
a4<-comprank(Duration,Infiltration,"spearman","dubois")
cat("Coefficient","method","\n",a1$r,a1$ities,"\n",a2$r,a2$ities,"\n",a3$r,a3$ities,
"\n",a4$r,a4$ities,"\n")
detach(Laudaher)
#####
# Asymptotic confidence intervals.
r.cofint <- function(r, n, index, Level=.95) {
	asd<-rep(0,4)
	asd[1]<-1/sqrt(n-1) # Spearman
	asd[2]<-sqrt((4*n+10)/(9*n*(n-1))) # Kendall
	asd[3]<-1/sqrt(1.5*n) # Gini
	asd[4]<-1/sqrt(1.00762*(n-1)) # r4
	# Fisher-Yates-means, Fisher-Yates-medians, symmetrical BZ
	if (index<=4) {zse<-r*asd[index]} else {zse<-atan(r*(1-0.6/(n+8)))/sqrt(n-3)} 
	rlow <- r - zse * qnorm((1-Level)/2,lower.tail=FALSE);rlow<-max(-1,rlow)
 	rupp <- r + zse * qnorm((1-Level)/2,lower.tail=FALSE);rupp<-min(1,rupp)
 	out<-list(Lower.r=rlow, Upper.r=rupp)
 	return(out)
}
#
# Rajalingam S. and Zeya O. "Summative assessments by Spearman`s correlation
# coefficient: a case study. in enhancing learning: teaching & learning international 
# conference. Nov 24-26 2011. Miri, Sarawak: Curtin University

FA=c(39.09,39.77,35.62,34.69,34.42,35.57,35.94,38.68,38.41,36.4,37.95,37.03,38.03,
35.5,38.55,30.68,26.3,36.28)
SA=c(27.5,34,24,24,17,17.5,16.5,26,25.5,28.5,26.5,13,12.5,9.5,28.5,23,24.5,22)
n<-length(FA)
op<-par(mfrow=c(1,1))
plot(FA,SA, main="Academic progresses of students",
xlab="Formative assessment", ylab="summative assessment",
pch=19,col="steelblue4")
text(FA,SA,label=1:n,cex=0.8,pos=2)
abline(h=mean(SA),col="darkblue",lty=2,lwd=1)
abline(v=mean(FA),col="darkblue",lty=2,lwd=1)
par(op)
rct<-c("spearman","kendall","gini","r4", "fy1","fy2","sbz")
for (index in 1:7){
	r<-comprank(FA,SA,rct[index])$r
	cir<-r.cofint(r,n,index,Level=0.99)
	cat(rct[index],"Low:",cir$Lower,"Value:",r," Upp: ",cir$Upper.r,"\n")
}
#####
#
# Daniel, C.  Wood, F. S. Fitting Equations to Data. New York: John Wiley,
# 1971, p. 45. Pilot-plant data
# The response variable (y) corresponds to the acid content determined by
# titration and the explanatory variable (x) is the organic acid content determined 
# by extraction and weighting
y<-c(76, 70, 55, 71, 55, 48, 50, 66, 41, 43, 82, 68, 88, 58, 64, 88, 89, 88,
84, 88)
x<-c(123, 109, 62, 104, 57, 37, 44, 100, 16, 28, 138, 105, 159, 75, 88, 164,
169, 167, 149, 167)
out1<-comprank(x,y,"sp","woodbury")
out2<-comprank(x,y,"ke","woodbury")
out3<-comprank(x,y,"gi","woodbury")
out4<-comprank(x,y,"r4","woodbury")
out5<-comprank(x,y,"fy1","woodbury")
out6<-comprank(x,y,"fy2","woodbury")
out7<-comprank(x,y,"sbz","woodbury")
ind<-c(out1$r,out2$r,out3$r,out4$r,out5$r,out6$r,out7$r)
cat(out1$ities,"\n")
print(round(ind,5))
#####
#
data(DietFish);attach(DietFish)
op<-par(mfrow=c(1,1))
plot(occ_D,occ_H, main="Comparison of the diets of banded killfish",xlab="occ F_diaphanus",
ylab="occ F_heteroclitus", pch=19, cex=0.8,col="darkcyan")
text(occ_D,occ_H, labels = rownames(DietFish), cex=0.5, pos=3)
abline(h=mean(occ_D),col="darkred",lty=2,lwd=1)
abline(v=mean(occ_H),col="darkred",lty=2,lwd=1)
par(op)
a<-comprank(occ_D,occ_H,"sp","woodbury");cat(a$r,a$ities,"\n")
a<-comprank(occ_D,occ_H,"sp","gh");cat(a$r,a$ities,"\n")
a<-comprank(occ_D,occ_H,"sp","wgh");cat(a$r,a$ities,"\n")
a<-comprank(occ_D,occ_H,"sp","midrank");cat(a$r,a$ities,"\n")
a<-comprank(occ_D,occ_H,"sp","dubois");cat(a$r,a$ities,"\n")
detach(DietFish)
#####
#
data(Radiation);attach(Radiation);Radiation<-as.matrix(Radiation)
r1<-comprank(Radiation,"spearman","midrank",print=FALSE);eigen(r1$r)
r2<-comprank(Radiation,"kendall","midrank",print=FALSE);eigen(r2$r)
r3<-comprank(Radiation,"gini","midrank",print=FALSE);eigen(r3$r)
r4<-comprank(Radiation,"r4","gh",print=FALSE);eigen(r4$r)
r5<-comprank(Radiation,"fy1","gh",print=FALSE);eigen(r5$r)
r6<-comprank(Radiation,"fy2","gh",print=FALSE);eigen(r6$r)
r7<-comprank(Radiation,"sbz","gh",print=FALSE);eigen(r7$r)
detach(Radiation)
#####
#
# Correlation matrix
data(Marozzi);attach(Marozzi)
Marozzi<-as.matrix(Marozzi)
cor1<-comprank(Marozzi,"spearman","midrank",print=FALSE)
rownames(cor1$r)<-colnames(Marozzi);rownames(cor1$r)<-colnames(Marozzi)
print(round(cor1$r,3))
cor1 <-comprank(Marozzi,"kendall","midrank",print=FALSE)
print(round(cor1$r,3))
cor1<-comprank(Marozzi,"gini","midrank",print=FALSE)
print(round(cor1$r,3))
cor1<-comprank(Marozzi,"r4","wgh",print=FALSE)
print(round(cor1$r,3))
cor1<-comprank(Marozzi,"fy1","wgh",print=FALSE)
print(round(cor1$r,3))
cor1<-comprank(Marozzi,"fy2","wgh",print=FALSE)
print(round(cor1$r,3))
cor1<-comprank(Marozzi,"sbz","wgh",print=FALSE)
print(round(cor1$r,3))
detach(Marozzi)
}
