\name{ranktes}
\alias{ranktes}
\title{Tests of independence between two rankings}
\description{Performs various independence tests based on rank correlations.}
\usage{ranktes(r, n, index = "", approx = "exact", CC = FALSE, type = "two-sided",
	 print = TRUE)}
\arguments{
  \item{r}{the value of the test statistic.} 
  \item{n}{the number of ranks.}
  \item{index}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman","kendall","gini","r4", "fy" (that is, Fisher-Yates correlation) and "fa" (Filliben-Amerise rank correlation). Only enough of the string to be unique is required.}  
  \item{approx}{a character string that specifies the type of approximation to the null distribution of the statistic required in \code{index}: "vggfr", "exact","gaussian","student". Only enough of the string to be unique is required.}  
 \item{CC}{if true, a continuity correction is applied. Ignored if \code{approx= "exact"} or if \code{index="r4"} or if \code{index="fy"} or if \code{index="fa"}.}
 \item{type}{type of alternative hypothesis. The options are \code{"two-sided"} (independence), 
   \code{"greater"} (concordance) or \code{"less"} (discordance). Only enough of the string to be unique is required.}
\item{print}{\code{FALSE} suppresses some of the output.}
}

\value{A list containing the following components:
\item{n}{number of ranks}
\item{statistic}{type of rank correlation}
\item{r}{observed value of the coefficient}
\item{approx}{type of approximation}
\item{tails}{type of alternative hypothesis}
\item{Cpv}{conservative \eqn{p}-value}
\item{Lpv}{liberal \eqn{p}-value}
\item{Lambda}{if \code{approx}=\code{"vggfr"} returns the vector containing the two shape parameters of the VGGFR density that best fits the null distribution of the required coefficients. It is \code{NULL} otherwise}
}
\details{
Upon computing \eqn{r_h}, it is common practice to determine whether the value  is large enough, in absolute terms, to lead to the conclusion that the rank correlation coefficient that would be obtained from the entire set of permutations, say it \eqn{\rho_h} is different from zero. To this end, we consider the test \eqn{H_0:\rho_h=0} against:
  
\eqn{H_1:\rho_h>0} (concordance). Only a large positive \eqn{r_h} can be considered in line with this alternative hypothesis. To be considered significant, \eqn{r_h} must be greater and its value must be equal to or larger than the critical values corresponding to the prespecified \eqn{\alpha} level of significance.

\eqn{H_1:\rho_h<0} (discordance). Only a large negative \eqn{r_h} value will provide support for this alternative hypothesis. More specifically, \eqn{r_h} must be less and \eqn{|r_h|} must be equal to or larger than the critical values corresponding to the prespecified \eqn{\alpha} level of significance.

\eqn{H_1:\rho_h\ne 0} (independence). Either a large negative or a large positive value of \eqn{r_h} are conform to this alternative hypothesis. A significant value of  \eqn{r_h} is obtained if the prespecified \eqn{\alpha} level of significance is equal to or greater than 2\eqn{Prob(-|r_h|)} where \eqn{Prob(.)} is the probability density used to approximate the null distribution of \eqn{r_h}.


	It is important to note that, as correctly observe Iman and Conover (1978), the discreteness of rank correlations often leads into situations where no critical region has exactly the size \eqn{\alpha}. If \code{approx}\code{="exact"}  the routine provides the next smaller exact size called conservative \eqn{p}-value or the next larger exact size called liberal \eqn{p}-value. A test can be considered conclusive if both the conservative and the liberal significance levels lie on the same side with respect of \eqn{\alpha}. Compared to an exact \eqn{p}-value, a conservative \eqn{p}-value tends to understate the evidence against \eqn{H_0}, whereas a liberal \eqn{p}-value tends to overstate it.

If \code{approx}\code{="exact"} then \code{ranktes} uses precomputed exact null distributions. In particular, \eqn{n\le 26} for \eqn{r_1} (see Gustafson, 2009); \eqn{n\le 24} for \eqn{r_3} (see Girone et al. 2010). We have obtained the null distribution of \eqn{r_4}, \eqn{r_5} and \eqn{r_6} for \eqn{n\le 15} by calculating rank correlations for all the \eqn{n!} permutations of the integers \eqn{1,2, \cdots, n}. Kendall's \eqn{r_2} benefits from a recurrence relationship (Panneton and Robillard, 1972) that can handle \eqn{r_2} for \eqn{n} up \eqn{60}. In this regard, it is very useful the package \code{Mpfr}.
 
In the case \code{approx} \eqn{\ne} \code{"exact"}, this routine computes critical values and significance levels by applying the density function indicated in \code{approx}. Now, of course,  conservative and liberal p-values coincide. 
 
	The statistics involved in t-Student approximations are:

\deqn{r_h^+=r_h\sqrt{\frac{m_{h,a}}{1-r_h^2}} \sim t_{m_{h,b}}, h=1,\cdots,4}

where \eqn{m_{1,a}=n-2}, \eqn{m_{2,a}=[9n(n-1)/(4n+10)-1]}, \eqn{m_{3,a}=[3(n-1)(n^2-k_n)/2(n^2+2+ k_n)]} with \eqn{k_n=n\ mod \ 2} and \eqn{m_{4,a}=2(n-2.01524)/1.00762}. Furthermore,  \eqn{m_{1,b}=n-2}, \eqn{m_{2,b}=\lfloor m_{2,a} \rfloor}, \eqn{m_{3,b}=\lfloor m_{3,b}+0.5 \rfloor}, \eqn{m_{4,b}=\lfloor m_{4,a} \rfloor}.

The t-Student approximation of \eqn{r_1} and the Gaussian approximation of Kendall's \eqn{r_2} are well known. Vittadini (1996) proposes the t-Student approximation to Kendall's \eqn{r_2}.  Landenna et. al. (1989) suggest the t-Student for the Gini cograduation coefficient \eqn{r_3}. Tarsitano and Amerise (2015) developed  a t-Student approximation to the null distribution of \eqn{r_4}. Terry (1952) has pointed out that the t-distribution with \eqn{n-2} degrees of freedom provides a good approximation to the null distribution of \eqn{r_5^+}, We conjecture that a similar result holds with Filliben-Amerise rank correlation \eqn{r_6^+}.

The statistics involved in the Gaussian distribution are:

\deqn{r_1^*=\frac{r_1}{\sqrt{n-1}}, \ r_2^*=r_2\sqrt{\frac{4n+10}{9n(n-1)}}, \ r_3^*=\frac{r_3}{\sqrt{1.5n}}, \  r_4^*=\frac{r_4}{\sqrt{1.00762(n-1)}}}

\deqn{r_5^*=\frac{r_5}{\sqrt{n-1}}, \ r_6^*=\frac{r_6}{\sqrt{n-1}}}

Cifarelli and Regazzini (1977) give the Gaussian approximation to the null distribution of Gini coefficient. See also Genest et al. (2010). The Gaussian approximation to the null distribution of \eqn{r_4} is obtained in Tarsitano and Amerise (2015). The Gaussian approximations to \eqn{r_5} and  \eqn{r_6} are asymptotic results.

If \code{approx="vggfr"} the Vianelli generalized Gaussian distribution with finite range is fitted by using the method of moments.
}

\author{Agostino Tarsitano, Ilaria Lucrezia Amerise, Marco Marozzi}

\note{In the case of a two-sided test, the \eqn{p}-value is defined as two times the one-sided \eqn{p}-value, bounded above by 1.}

\references{

Cifarelli, D. M. and Regazzini, E. (1977). "On a distribution-free test of independence based on Gini's rank association coefficient". \emph{Recent Developments in Statistics (Proceedings of the European Meeting of Statisticians, Grenoble, 1976)}, Amsterdam, North-Holland, 375--385.

Genest, N. B. and Neslehova, J. and Ben Ghorbal, N.  (2010). "'Spearman's footrule and Gini's gamma: a review with complements" \emph{Journal of Nonparametric Statistics}, \bold{22}, 937--954.

Girone, G. et al. (2010). "La distribuzione campionaria dell'indice di cograduazione di Gini per dimensioni campionarie fino a 24". \emph{Annali del Dipartimento di Scienze Statistiche "Carlo Cecchi" - Universita` di Bari}, \bold{24}, 246--271.

Gustafson, L. (2009).  rho null distribution. Available at \url{http://www.luke-g.com/}

Iman, L. and Conover, W. J. (1978). "Approximations of the critical region for 's rho with and without ties presen"t. \emph{Communication in Statistics - Simulation and Computation}, \bold{7}, 269--282.

Landenna, G. and Scagni, A. and Boldrini, M. (1989). "An approximated distribution of the Gini's rank association coefficient". \emph{Communications in Statistics. Theory and Methods}, \bold{18}, 2017--2026.

Tarsitano, A. and Amerise, I. L. (2015). "On a measure of rank-order association".  \emph{Journal of Statistical and Econometric Methods}, \bold{4}, 83--105.

Tarsitano, A. and Amerise, I. L. (2016). "Modelling the null distribution of rank correlations". Submitted.

Terry, M. E. (1952). "Some rank order tests which are most powerful against specific parameteric alternatives". \emph{Annals of Mathematical Statistics},  \bold{23}, 346--366.

Vittadini, G. (1996). "Una famiglia di distribuzione per i test di associazione". In \emph{Atti della XXXVIII riunione scientifica della S.I.S.}, Rimini 9-13 Aprile 1996, 2, 521--528.

}

\examples{
	
# G. P. Watkins (1933). An Ordinal Index of Correlation, Journal of the 
# American Statistical Association, 28:182, 139-151.
#  20-item series for area and density have been made up to cover the
# original 13 states and the four others earliest admitted to the Union.
State<-c("Georgia","North_Carolina","New_York","Luisiana","Pennsylvania",
"Virginia","Tennessee","Ohio","Kentucky","Maine","South_Carolina",
"West_Virginia","Maryland","Vermont","New_Hampshire","Massachusetts",
"New_Jersey","Connecticut","Delaware","Rhode_Island")
Area<-c(1,2,	3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
Density<-c(17,13,5,18,6,14,12,7,11,20,15,10,8,19,16,2,3,4,9,1)
op<-par(mfrow=c(1,1))
plot(Area,Density,main="",xlab="Area",ylab="Density",pch=19,cex=0.9,
col="darkgreen" )
abline(h=mean(Area),col="black",lty=2,lwd=1)
abline(v=mean(Density),col="darkblue",lty=2,lwd=1)
par(op)
r<-comprank(Area,Density,"fa","wgh")$r
ranktes(r, length(Area), "fa", "ga",FALSE, "two", TRUE)
#####
#
\dontrun{
data(Atar);attach(Atar)
op<-par(mfrow=c(1,1))
plot(TBL,TFL,main="",xlab="Backward Linkage Index",ylab=
"Forward Linkage Index",pch=19, cex=0.9,col="magenta")
abline(h=mean(TFL),col="black",lty=2,lwd=1)
abline(v=mean(TBL),col="black",lty=2,lwd=1)
par(op)
r<-comprank(TBL,TFL,"fy","wgh")$r
ranktes(r, length(TBL), "fy", "vggfr",FALSE, "two", TRUE)
detach(Atar)
}
#####
#
data(Sharpe);attach(Sharpe)
op<-par(mfrow=c(1,1))
plot(AVR,VAR, type = "p",pch=19,cex=1.1,col="tomato",main="Mutual fund 
performance")
text(AVR,VAR, labels = rownames(Sharpe), cex=0.5, pos=3)
abline(h=mean(AVR),col="black",lty=2,lwd=1)
abline(v=mean(VAR),col="black",lty=2,lwd=1)
par(op)
r<-comprank(AVR,VAR,"gini","wgh")$r
ranktes(r, length(AVR), "gini", "st",FALSE, "greater", TRUE)
detach(Sharpe)
#####
#
\dontrun{
# Sun,J.-G. and Jurisicova, A. and Casper, R.F. (1997). "Detection of
# Deoxyribonucleic Acid Fragmentation in Human Sperm: Correlation
# with Fertilization In Vitro".  Biology of Reproduction, 56, 602-607.
n<-c(222,298,143,143,291,148); r<-c(-0.18,-0.12,-0.16,-0.20,-0.06,-0.003)
App<-c("Ga","St","Vg")
N<-length(n);Ta<-matrix(NA,N,5)
for (i in 1:length(n)){Ta[i,1]<-r[i];Ta[i,2]<-n[i]
		for (j in 1:3){
				app<-App[j]
				a<-ranktes(r[i],n[i],"S",app,FALSE,"t",FALSE);Ta[i,2+j]<-a$Cpv
		}}
Df<-matrix(Ta,6,5)
rownames(Df)<-c("Conc. sperm/mL","Motility", "Fertilization rate",
"Cleavage rate", "Male age","Abstinence days")
colnames(Df)<-c("Spearman","n of samples","Appr. Gaussian",
"Appr. t-Student", "Appr. GGFR")
Df<-as.data.frame(Df); print(round(Df,5))
}
#####
#
\dontrun{
data(Starshi);attach(Starshi)
op<-par(mfrow=c(1,1))
plot(Sm15F,Sm15M, type = "p",pch=19,cex=0.9,col="darkorange",
main="Smokers ")
text(Sm15F,Sm15M,labels = rownames(Starshi),cex=0.6,pos=
c(1,rep(2,10),3,2))
abline(h=mean(Sm15M),col="black",lty=2,lwd=1)
abline(v=mean(Sm15F),col="black",lty=2,lwd=1)
par(op)
r<-comprank(Sm15F,Sm15M,"r4","wgh")$r
a<-ranktes(r, length(Sm15F), "r4", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"sp","wgh")$r
a<-ranktes(r, length(Sm15F), "sp", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"ken","wgh")$r
a<-ranktes(r, length(Sm15F), "ken", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"g","wgh")$r
a<-ranktes(r, length(Sm15F), "g", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"fy","wgh")$r
a<-ranktes(r, length(Sm15F), "fy", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
r<-comprank(Sm15F,Sm15M,"fa","wgh")$r
a<-ranktes(r, length(Sm15F), "fa", "ex",TRUE, "two", FALSE)
cat(a$Value,a$Cpv,a$Lpv,"\n")
detach(Starshi)
}
#####
#
\dontrun{
All.App<-function(r,n,index,type){
# Computes p-values of an observed rank correlation statistic
	A<-rep(r,9)
	names(A)<-encodeString(c(index,"t-Student, CC=F","Gaussian, CC=F", "VGGFR,
	CC=F","t-Student, CC=T","Gaussian, CC=T", "VGGFR, CC=T",
	"Exact p Conservative","Exact p Liberal"),justify="right")
		a<-ranktes(r,n,index,"St",FALSE,type,FALSE);A[2]<-a$Cpv
		a<-ranktes(r,n,index,"Ga",FALSE,type,FALSE);A[3]<-a$Cpv
		a<-ranktes(r,n,index,"Vg",FALSE,type,FALSE);A[4]<-a$Cpv
		a<-ranktes(r,n,index,"St",TRUE,type,FALSE);A[5]<-a$Cpv
		a<-ranktes(r,n,index,"Ga",TRUE,type,FALSE);A[6]<-a$Cpv
		a<-ranktes(r,n,index,"Vg",TRUE,type,FALSE);A[7]<-a$Cpv
		a<-ranktes(r,n,index,"Ex",FALSE,type,FALSE);A[8]<-a$Cpv;A[9]<-a$Lpv
		A<-as.matrix(A)
return(A)}
data(Gabbs);attach(Gabbs)
B<-matrix(0,9,6)
colnames(B)<-colnames(Gabbs[1:6])
rownames(B)<-encodeString(c("index","t-Student, CC=F","Gaussian, CC=F",
"VGGFR, CC=F","t-Student, CC=T", "Gaussian, CC=T", "VGGFR, CC=T",
"Exact p Conservative","Exact p Liberal"), justify="right")
index<-"spearman"
for (i in 1:6){r<-comprank(Gabbs[,i],Gabbs[,7],index, print=FALSE)$r
					B[,i]<-All.App(r,19,index,"less")}
print(round(B,5))
detach(Gabbs)
}
#####
#
\dontrun{
data(Dalyww);attach(Dalyww)
op<-par(mfrow=c(1,1))
plot(ACLS,ASHR,main="The paradox of high rates of suicide in happy places",
xlab="Adjusted Life Satisfaction", ylab="Adjusted Suicide Risk",pch=19,
cex=0.8,col="steelblue")
text(ACLS,ASHR,labels=rownames(Dalyww),cex=0.7,pos=2)
abline(h=mean(ASHR),col="black",lty=2,lwd=1)
abline(v=mean(ACLS),col="black",lty=2,lwd=1)
par(op)
r<-comprank(ACLS,ASHR,"spearman")$r;n<-length(ASHR)
out<-ranktes(r,n,"s","ga",FALSE,"greater",FALSE)
cat(round(out$Value,3),round(out$Cpv,5),round(out$Lpv,5),"\n")
out<-ranktes(r,n,"s","st",FALSE,"greater",FALSE)
cat(round(out$Value,3),round(out$Cpv,5),round(out$Lpv,5),"\n")
out<-ranktes(r,n,"s","vg",FALSE,"greater",FALSE)
cat(round(out$Value,3),round(out$Cpv,5),round(out$Lpv,5),"\n")
#
r<-comprank(ACLS,ASHR,"kendall")$r
out<-ranktes(r,n,"kendall","st",FALSE,"greater",FALSE)
cat(round(out$Value,3),round(out$Cpv,5),round(out$Lpv,5),"\n")
#
r<-comprank(ACLS,ASHR,"r4")$r
out<-ranktes(r,n,"r4","st",FALSE,"greater",FALSE)
cat(round(out$Value,3),round(out$Cpv,5),round(out$Lpv,5),"\n")
detach(Dalyww)
}
}