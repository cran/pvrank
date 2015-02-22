\name{rankor}
\alias{rankor}
\title{Computes values and p-values of some rank correlations}
\description{Computes rank correlations and evaluates exact and estimated \eqn{p}-values under the null hypothesis of no association. If there are no ties, all the tests are exact, at least within the limitations imposed by the complete enumeration of all permutations. When there are ties, the function offers two solutions as described below.}
\usage{rankor(x, y, index = "spearman", approx = "exact", tiex = "woodbury", CC = FALSE, 
	type = "two-sided", print = TRUE, sizer = 1000000, repl = 1000)}
\arguments{
 \item{x}{a numeric vector.}  
 \item{y}{a numeric vector with compatible dimension to \eqn{x}.}  
 \item{index}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman","kendall","gini" and "r4". Only enough of the string to be unique is required.}  
 \item{approx}{a character string that specifies the type of approximation to the null distribution: "vggfr", "exact","gaussian" and "student". Only enough of the string to be unique is required. Exact computations use frequencies precomputed and stored in specific look-up tables.} 
 \item{tiex}{a character string that specifies the method for breaking ties. Possible options are: "woodbury", "gh","wgh", "midrank","dubois". Only enough of the string to be unique is required. Ignored if there are no equal scores. The options: \code{ "midrank"} and \code{"dubois"} are also ignored in the calculation of \eqn{r_4}.}
 \item{CC}{if true, a continuity correction is applied. Ignored if \code{approx}= \code{"exact"} and if \code{index}=\code{"r4"}.}
 \item{type}{type of alternative hypothesis. The options are \code{"two-sided"} (default), \code{"greater"} or \code{"less"}. Only enough of the string to be unique is required; \code{"greater"} corresponds to positive association, \code{"less"} to negative association.}
\item{print}{\code{FALSE} suppresses some of the output.}
 \item{sizer}{number of replications for resolving ties by randomization. Default value: \eqn{10^6}. Use of this specification is recommended for applying the Woodury method to the Gini index.}
\item{repl}{number of sampled pairs of permutations needed to apply the weighted max-min method. Default value: \eqn{10^3}.}
}

\details{
	
The presence of ties can be dealt with two approaches. (i) Use the fact that for any given number of ties of each multiplicity the mean and the variance of the statistics will tend towards that of the system without any ties as \eqn{n} increases. Consequently, the standard Gaussian distribution can be used after resolving ties by the method of mid-rank or the method of Dubois. (ii) Determine a synthesis of the values of the given statistic over all possible ways in which ranks could have been caused by truly differing observations. This solution enables the execution of rank correlation tests by using the usual procedure applied in the absence of ties. See the description of \code{comprank}.

It can be noted that each one of \eqn{r_h, h=1,2,3} includes a quantity \eqn{S_h}, which is always an integer. It is possible to add a continuity correction for \eqn{r_h, h=1,2,3} to adjust the fact that we are approximating a discrete distribution with a continuous one. Kendall and Dickinson Gibbons (1990) [p. 65] improved the approximation of \eqn{r_2} by subtracting one from the observed \eqn{S_2} if it is positive or adding one if it is negative. Also, Kendall and Dickinson Gibbons (1990) [p. 70] suggested subtracting one to \eqn{S_1} if it is less than \eqn{(n^3-n)/6} and adding one if it exceeds \eqn{(n^3-n)/6}. The correction to \eqn{S_3} can follow the same scheme as that of \eqn{S_1}. The formulation of coefficients leads to the recommendation of not using the correction for continuity with \eqn{r_4}.

If \eqn{p<0.5} the use of the continuity correction decreases the \eqn{p}-value, but without the correction the calculations are slanted in favor of not rejecting the null hypothesis. If \eqn{n} is sufficiently large, the effect of  the continuity correction is minimal. In this regard, Haber (1982) notes that, while it is agreed that continuity correction usually improves the approximation to the cumulative distribution function of a discrete random variable the use of continuity correction when performing a statistical test is still an issue of major controversy.  See, for example, Kruskal & Wallis(1952). 

}
\value{A list containing the following components:
\item{n}{number of ranks}
\item{Statistic}{coefficient of monotone association}
\item{r}{observed value}
\item{approx}{type of approximation}
\item{tails}{the type of alternative}
\item{Cpv}{Approximate or exact conservative p-value}
\item{Lpv}{Approximate or exact or liberal p-value}
}
\author{Agostino Tarsitano, Ilaria Lucrezia Amerise, Marco Marozzi}
\note{
	
The approximations to the exact null distribution of the rank correlations statistics considered in the \pkg{pvrank} are discussed for the command \code{"ranktes"}.
}

\references{
	
Haber, M. (1982). The continuity correction and statistical testing. \emph{International Statistical Review}, \bold{50}, 135-144.

Kendall, M. and Dickinson Gibbons, J. (1990). \emph{Rank Correlation Methods, 5th ed.}, Oxford University Press,  New York. 

Kruskal, W.H. and Wallis, A. (1952). Use of ranks in one-criterion variance analysis. \emph{Journal of the American Statistical Association}, \bold{47}, 583-621
}

\examples{
data(Lambh);attach(Lambh)
rankor(DVI,Temp,"spearman","ga","gh",FALSE,"two",TRUE)
rankor(DVI,Temp,"r4","ga","gh",FALSE,"two",TRUE)
detach(Lambh)	
###
data(Berk);attach(Berk)
op<-par(mfrow=c(1,2), mgp=c(1.8,.5,0), mar=c(2.8,2.7,2,1),oma=c(0,0,0,0))
plot(Births,Deaths,main="",pch=19,cex=0.99,cex.lab=0.99,cex.axis=0.8)
abline(h=mean(Deaths),col="black",lty=2,lwd=1)
abline(v=mean(Births),col="black",lty=2,lwd=1)
text(Births[12],Deaths[12],labels="noon",pos=3,cex=0.7)
text(Births[24],Deaths[24],labels="midnight",pos=4,cex=0.7)
W<-Births[-c(12,24)];Z<-Deaths[-c(12,24)]
plot(W,Z,main="",xlab="Births",ylab="Deaths",pch=19,cex=0.99,cex.lab=0.99,
cex.axis=0.8)
abline(h=mean(Z),col="black",lty=2,lwd=1)
abline(v=mean(W),col="black",lty=2,lwd=1)
A<-matrix(NA,10,5);Series<-c("Complete","Clean")
colnames(A)<-c("Data set","Coeff.","Value", "Cons. Two-tail p","Lib. Two-tail p")
a0<-cor.test(Births,Deaths, method = "pearson", alternative = "t")
k<-1;A[k,1]<-Series[1];A[k,2]<-"pearson";A[k,3]<-round(a0$estimate,4)
A[k,4]<-round(a0$p.value,5);A[k,5]<-round(a0$p.value,5)
for (j in c("spearman","kendall","gini","r4")){k<-k+1
		a<-rankor(Births,Deaths,j,"st","gh",FALSE,"two",FALSE)
		A[k,1]<-Series[1];A[k,2]<-j;A[k,3]<-round(a$Value,4);A[k,4]<-round(a$Cpv,5)
		A[k,5]<-round(a$Lpv,5)}
a1<-cor.test(W,Z, method = "pearson", alternative = "t")		
k<-k+1;A[k,1]<-Series[2];A[k,2]<-"pearson";A[k,3]<-round(a1$estimate,4)
A[k,4]<-round(a1$p.value,5);A[k,5]<-round(a1$p.value,5)
for (j in c("spearman","kendall","gini","r4")){k<-k+1
		a<-rankor(W,Z,j,"st","wgh",FALSE,"two",FALSE)
		A[k,1]<-Series[2];A[k,2]<-j;A[k,3]<-round(a$Value,4);A[k,4]<-round(a$Cpv,5)
		A[k,5]<-round(a$Lpv,5)}
A<-as.data.frame(A)
print(A,print.gap=4,right=FALSE)
detach(Berk)
###
data(Locre);attach(Locre)
op<-par(mfrow=c(1,1))
plot(Males,Females,main="Fer cryin' out loud - there is a sex difference",xlab="Females",
       ylab="Males",pch=19,cex=0.8,col="steelblue")
text(Males,Females,labels=1:length(Females),cex=0.7,pos=2)
abline(h=mean(Females),col="black",lty=2,lwd=1)
abline(v=mean(Males),col="black",lty=2,lwd=1)
par(op)
out<-rankor(Males,Females,"g","vg","gh",FALSE,"two")
cat(out$Value,out$Cpv,"\n")
cor.test(Males,Females, alternative="two.sided",method="pearson", continuity= FALSE)
detach(Locre)
###
# Daniel, C.  Wood, F. S. Fitting Equations to Data. New York: John Wiley, 1971, p. 45
# Pilot-plant data
# The response variable (y) corresponds to the acid content determined by titration and
# the explanatory variable (x) is the organic acid content determined by extraction and
# weighting
y<-c(76, 70, 55, 71, 55, 48, 50, 66, 41, 43, 82, 68, 88, 58, 64, 88, 89, 88, 84, 88)
x<-c(123, 109, 62, 104, 57, 37, 44, 100, 16, 28, 138, 105, 159, 75, 88, 164, 169, 167, 
149, 167)
out<-rankor(x,y,"s","ex","dubois");out
cat(out$Value,out$Cpv,"\n")
}
