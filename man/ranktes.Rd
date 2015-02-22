\name{ranktes}
\alias{ranktes}
\title{Tests of independence between two rankings}
\description{Performs various independence tests based on rank correlations.}
\usage{ranktes(r, n, index = "spearman", approx = "exact", CC = FALSE, 
	type = "two-sided", print = TRUE)}
\arguments{
  \item{r}{the value of the test statistic.} 
  \item{n}{the number of ranks.}
  \item{index}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman","kendall","gini","r4". Only enough of the string to be unique is required.}  
  \item{approx}{a character string that specifies the type of approximation to the null distribution of the statistic required in \code{index}: "vggfr", "exact","gaussian","student". Only enough of the string to be unique is required.}  
 \item{CC}{if true, a continuity correction is applied. Ignored if \code{approx= "exact"} or if \code{index="r4"}.}
 \item{type}{type of alternative hypothesis. The options are \code{"two-sided"} (default), \code{"greater"} or \code{"less"}. Only enough of the string to be unique is required; \code{"greater"} corresponds to positive association, \code{"less"} to negative association.}
\item{print}{\code{FALSE} suppresses some of the output.}
}

\value{A list containing the following components:
\item{n}{number of ranks}
\item{Statistic}{coefficient of monotone association}
\item{r}{observed value}
\item{approx}{type of approximation}
\item{tails}{type of alternative hypothesis}
\item{Cpv}{Conservative \eqn{p}-value}
\item{Lpv}{Liberal \eqn{p}-value}
\item{Lambda}{If \code{approx}=\code{"vggfr"} returns the vector containing the two shape parameters of the VGGFR density that best fits the null distribution of the required coefficients. It is \code{NULL} otherwise}
}
\details{

It is important to note that, as correctly observe Iman and Conover (1978), the discreteness of rank correlations often leads into situations where no critical region has exactly the size \eqn{\alpha}. If \code{approx}\code{="exact"}  the routine provides the next smaller exact size called conservative \eqn{p}-value or the next larger exact size called liberal \eqn{p}-value. A test can be considered conclusive if both the conservative and the liberal significance levels lie on the same side with respect of \eqn{\alpha}.

If \code{approx}\code{="exact"} then \code{rankes} uses precomputed exact null distributions for \eqn{n\le 26} (Spearman or \eqn{r_1}; see Gustafson, 2009). \eqn{n\le 24} (Gini or  \eqn{r_3}; see Girone et al. 2010). \eqn{n\le60} (Kendall or \eqn{r_2}; here we apply the recursive formula given by Panneton and Robillard, 1972).  \eqn{n\le 15} (\eqn{r_4}); see Tarsitano and Amerise (2015).

 In the case \code{approx}\eqn{\ne} \code{"exact"}, this routine computes critical values and significance levels by applying the density required in \code{approx}. Now, of course,  conservative and liberal p-values coincide. 
 
The statistics involved in the t-Student approximation are

\deqn{r_h^+=r_h\sqrt{\frac{m_h}{1-r_h^2}} \sim t_{m_h} \ h=1,\cdots,4}.

where \eqn{m_1=n-2}, \eqn{m_2=[9n(n-1)/(4n+10)-1]}, \eqn{m_3=[3(n-1)(n^2-k_n)/2(n^2+2+ k_n)]}, \eqn{m_4=(n-2.01524)/1.00762} with \eqn{k_n=n\ mod \ 2}.

The statistics involved in the Gaussian distribution are
\deqn{
r_1^*=r_1\sqrt{n-1}, \ r_2^*=r_2\sqrt{\frac{9n(n-1)}{4n+10}}, \ r_3^*=r_3\sqrt{1.5n}, \ r_4^*=r_4\sqrt{1.00762(n-1)}.
}

The t-Student approximation of Spearman's \eqn{r_1} and the Gaussian approximation of Kendall's \eqn{r_2} are well known. 

Vittadini (1996) proposes the t-Student approximation for Kendall's tau.  Landenna et. al. (1989) suggest the symmetrical t-Student for the Gini cograduation coefficient \eqn{r_3}.

Cifarelli and Regazzini (1977) give the Gaussian approximation to the null distribution of Gini coefficient. See also Genest et al. (2010).

The Gaussian approximation to the null distribution of \eqn{r_4} is obtained in Tarsitano and Amerise (2015). The same authors propose a t-Student approximation.

 If \code{approx="vggfr"} the Vianelli generalized Gaussian distribution with finite range is fitted by using the method of moments. The \code{approx="vggfr"}approximation to the null distribution of \eqn{r_4} seems to be very poor. Presumably this is an effect the reduced number of ranks for which the exact frequencies of \eqn{r_4} are fully known, which prevent us from having better possibilities for the fitting. The combination  \code{approx="vggfr"} and  \code{index="r4"} is not implemented at the moment.
}

\author{Agostino Tarsitano, Ilaria Lucrezia Amerise, Marco Marozzi}

\note{In the case of a two-sided test, the \eqn{p}-value is defined as two times the one-sided \eqn{p}-value, bounded above by 1.}

\references{

Cifarelli, D. M. and Regazzini, E. (1977). On a distribution-free test of independence based on Gini's rank association coefficient. \emph{Recent Developments in Statistics (Proceedings of the European Meeting of Statisticians, Grenoble, 1976)}, Amsterdam, North-Holland, 375--385.

Genest, N. B. and Neslehova, J. and Ben Ghorbal, N.  (2010). "Spearman's footrule and Gini's gamma: a review with complements" \emph{Journal of Nonparametric Statistics}, \bold{22}, 937--954.

Girone, G. et al. (2010). La distribuzione campionaria dell'indice di cograduazione di Gini per dimensioni campionarie fino a 24. \emph{Annali del Dipartimento di Scienze Statistiche "Carlo Cecchi" - Universita` di Bari}, \bold{24}, 246--271.

Gustafson, L. (2009). Spearman rho null distribution. Available at \url{http://www.luke-g.com/math/spearman/index.html}

Iman, L. and Conover, W. J. (1978). Approximations of the critical region for Spearman's rho with and without ties present. \emph{Communication in Statistics - Simulation and Computation}, \bold{7}, 269-282.

Landenna, G. and Scagni, A. and Boldrini, M. (1989). An approximated distribution of the Gini's rank association coefficient. \emph{Communications in Statistics. Theory and Methods}, \bold{18}, 2017-2026.

Tarsitano, A. and Amerise, I. L. (2013). "Approximation of the null distribution of rank correlations". Submitted.

Tarsitano, A. and Amerise, I. L. (2015). "On a measure of rank order association". Submitted.

Vittadini, G. (1996). "Una famiglia di distribuzione per i test di associazione". In \emph{Atti della XXXVIII riunione scientifica della S.I.S.}, Rimini 9-13 Aprile 1996, 2, 521--528.
}
\examples{
data(Atar);attach(Atar)
op<-par(mfrow=c(1,1))
plot(TBL,TFL,main="",xlab="Backward Linkage Index",ylab="Forward Linkage Index",pch=19,
cex=0.9,col="olivedrab")
abline(h=mean(TFL),col="black",lty=2,lwd=1)
abline(v=mean(TBL),col="black",lty=2,lwd=1)
par(op)
r<-comprank(TBL,TFL,"spearman","wgh")
ranktes(r, length(TBL), "spearman", "st",FALSE, "two", TRUE)
detach(Atar)
###
data(Sharpe);attach(Sharpe)
op<-par(mfrow=c(1,1))
plot(AVR,VAR, type = "p",pch=19,cex=1.1,col="violetred",main="Mutual fund performance")
text(AVR,VAR, labels = rownames(Sharpe), cex=0.5, pos=3)
abline(h=mean(AVR),col="black",lty=2,lwd=1)
abline(v=mean(VAR),col="black",lty=2,lwd=1)
par(op)
r<-comprank(AVR,VAR,"gini","wgh")
ranktes(r, length(AVR), "gini", "st",FALSE, "two", TRUE)
detach(Sharpe)
###
# Sun,J.-G. and Jurisicova, A. and Casper, R.F. (1997). "Detection of Deoxyribonucleic 
# Acid Fragmentation in Human Sperm: Correlation with Fertilization In Vitro".
# Biology of Reproduction, 56, 602-607.
n<-c(222,298,143,143,291,148);
r<-c(-0.18,-0.12,-0.16,-0.20,-0.06,-0.003)
App<-c("Ga","St","Vg")
N<-length(n)
Ta<-matrix(NA,N,5)
for (i in 1:length(n)){
		Ta[i,1]<-r[i];Ta[i,2]<-n[i]
		for (j in 1:3){
				app<-App[j]
				a<-ranktes(r[i],n[i],"S",app,FALSE,"t",FALSE);Ta[i,2+j]<-a$Cpv
		}}
Df<-matrix(Ta,6,5)
rownames(Df)<-c("Conc. sperm/mL","Motility", "Fertilization rate", "Cleavage rate",  
"Male age","Abstinence days")
colnames(Df)<-c("Spearman's rho","n of samples","Appr. Gaussian","Appr. t-Student", 
"Appr. GGFR")
Df<-as.data.frame(Df)
print(round(Df,5))
###
data(Starshi);attach(Starshi)
op<-par(mfrow=c(1,1))
plot(Sm15F,Sm15M, type = "p",pch=19,cex=0.9,col="violetred",main="Smokers ")
text(Sm15F,Sm15M,labels = rownames(Starshi),cex=0.6,pos=c(1,rep(2,10),3,2))
abline(h=mean(Sm15M),col="black",lty=2,lwd=1)
abline(v=mean(Sm15F),col="black",lty=2,lwd=1)
par(op)
r<-comprank(Sm15F,Sm15M,"r4","wgh")
a<-ranktes(r, length(Sm15F), "r4", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"s","wgh")
a<-ranktes(r, length(Sm15F), "s", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"k","wgh")
a<-ranktes(r, length(Sm15F), "k", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"g","wgh")
a<-ranktes(r, length(Sm15F), "g", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
detach(Starshi)
}