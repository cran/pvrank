\name{rankor}
\alias{rankor}
\title{Critical values and p-values of some rank correlations}
\description{Computes rank correlations and evaluates exact and estimated \eqn{p}-values under the null hypothesis of no association. If there are no ties, all the tests are exact, at least within the limitations imposed by the complete enumeration of all permutations. When there are ties, the function offers two solutions as described below.}
\usage{rankor(p, q, index = "spearman", approx = "exact", tiex = "woodbury", CC = FALSE, 
	type = "two-sided", print = TRUE, sizer = 100000, repgin = 1000)}
\arguments{
 \item{p}{a numeric vector.}  
 \item{q}{a numeric vector with compatible dimension to \eqn{p}.}  
 \item{index}{a character string that specifies the rank correlation that is used in the test statistic. Acceptable values are: "spearman","kendall","gini", "r4", "fy" (Fisher-Yates correlation), and "fa" (Filliben-Amerise rank correlation). Only enough of the string to be unique is required.}  
 \item{approx}{a character string that specifies the type of approximation to the null distribution. Acceptable values are: "vggfr", "exact","gaussian" and "student". Only enough of the string to be unique is required. Exact computations use frequencies precomputed and stored in specific look-up tables.} 
 \item{tiex}{a character string that specifies the method for breaking ties. Possible options are: "woodbury", "gh","wgh", "midrank", "dubois". Only enough of the string to be unique is required. Ignored if there are no equal scores.}
 \item{CC}{if true, a continuity correction is applied.}
 \item{type}{type of alternative hypothesis. The options are \code{"two-sided"} (default), \code{"greater"} or \code{"less"}. Only enough of the string to be unique is required; \code{"greater"} corresponds to positive association, \code{"less"} to negative association.}
\item{print}{\code{FALSE} suppresses some of the output.}
\item{sizer}{number of replications for resolving ties by randomization. Default value: \eqn{10^5}. Use of this specification is recommended for applying the Woodury method to the Gini index.}
\item{repgin}{number of sampled pairs of permutations needed to apply the weighted max-min method (wgh). Default value: \eqn{10^3}.}
}

\details{
	
The presence of ties is assessed through two approaches. (i) Use the fact that for any given number of ties of each multiplicity the mean and the variance of the statistics will tend towards that of data without any ties as \eqn{n} increases. Consequently, the standard Gaussian distribution can be used after resolving ties by the method of mid-rank or the method of Dubois. (ii) Determine a synthesis of the values of the given statistic over all possible ways in which ranks could have been caused by truly differing observations. This solution enables the execution of rank correlation tests by using the usual procedure applied in the absence of ties. See the description of \code{comprank}.

It can be noted that each one of \eqn{r_h, h=1,2,3} includes a quantity \eqn{S_h}, which is always an integer. It is possible to add a continuity correction for \eqn{r_h, h=1,2,3} to adjust the fact that we are approximating a discrete distribution with a continuous one. Kendall and Dickinson Gibbons (1990) [p. 65] improved the approximation of \eqn{r_2} by subtracting one from the observed \eqn{S_2} if it is positive or adding one if it is negative. Also, Kendall and Dickinson Gibbons (1990) [p. 70] suggested subtracting one to \eqn{S_1} if it is less than \eqn{(n^3-n)/6} and adding one if it exceeds \eqn{(n^3-n)/6}. The correction to \eqn{S_3} can follow the same scheme as that of \eqn{S_1}. The formulation of coefficients leads to the recommendation of not using the correction for continuity with \eqn{r_4}, \eqn{r_5} and \eqn{r_6}. Furthermore, Iman and Conover (1978) point out that the continuity correction should not te considered with ties present.

If \eqn{p<0.5} the use of the continuity correction decreases the \eqn{p}-value, but without the correction the calculations are slanted in favor of not rejecting the null hypothesis. If \eqn{n} is sufficiently large, the effect of  the continuity correction is minimal. In this regard, Haber (1982) notes that, while it is agreed that continuity correction usually improves the approximation to the cumulative distribution function of a discrete random variable the use of continuity correction when performing a statistical test is still an issue of major controversy.  See, for example, Kruskal & Wallis (1952). 

}
\value{A list containing the following components:
\item{n}{number of ranks}
\item{statistic}{type of coefficient of rank correlation}
\item{r}{observed value of the coefficient}
\item{approx}{type of approximation}
\item{tails}{the type of alternative}
\item{Cpv}{approximate or exact conservative p-value}
\item{Lpv}{approximate or exact or liberal p-value}
}
\author{Agostino Tarsitano, Ilaria Lucrezia Amerise, Marco Marozzi}

\note{
The approximations to the exact null distribution of therank correlations statistics considered in \pkg{pvrank} are discussed with regard the function \code{"ranktes"}.
}

\references{
	
Haber, M. (1982). "The continuity correction and statistical testing". \emph{International Statistical Review}, \bold{50}, 135--144.

Iman, R. L. and Conover, W. J. (1978). "Approximations of the critical region for spearman`s rho with and without ties present". \emph{Communications in Statistics-Simulation and Computation}, \bold{7}, 269-282.

Kendall, M. and Dickinson Gibbons, J. (1990). \emph{Rank Correlation Methods, 5th ed.}, Oxford University Press,  New York. 

Kruskal, W.H. and Wallis, A. (1952). "Use of ranks in one-criterion variance analysis". \emph{Journal of the American Statistical Association}, \bold{47}, 583--621.

}

\examples{
	
data(Lambh);attach(Lambh)
rankor(DVI,Temp,"spearman","ga","gh",FALSE,"less",TRUE)
rankor(DVI,Temp,"r4","ga","gh",FALSE,"less",TRUE)
rankor(DVI,Temp,"fy","ga","gh",FALSE,"less",TRUE)
detach(Lambh)
#####
#
data(Security);attach(Security)
rankor(AP,OV,"sp","ga","gh",FALSE,"greater",TRUE)
rankor(MP,AP,"fy","ga","gh",FALSE,"greater",TRUE)
rankor(OV,MP,"fa","ga","gh",FALSE,"greater",TRUE)
detach(Security)
#####
#
data(Dizytwin);attach(Dizytwin)
op<-par(mfrow=c(1,1), mgp=c(1.8,.5,0), mar=c(2.8,2.7,2,1),oma=c(0,0,0,0))
plot(Latitude,DZT_Rate,main="Latitude and dizygotic twinning rates",xlab="Latitude",
ylab="DZT_Rate", pch=19,cex=0.9, col= "tan4")
text(Latitude,DZT_Rate,labels=rownames(Dizytwin),cex=0.6,pos=c(rep(3,10),1,3,1,
rep(3,4),1.3))
abline(v=mean(Latitude),col="black",lty=2,lwd=1)
abline(h=mean(DZT_Rate),col="darkblue",lty=2,lwd=1)
par(op)
rankor(Latitude,DZT_Rate,"r4","vg","gh",FALSE,"two",TRUE)
rankor(Latitude,DZT_Rate,"sp","ex","gh",FALSE,"two",TRUE)
rankor(Latitude,DZT_Rate,"ken","ex","gh",FALSE,"two",TRUE)
detach(Dizytwin)	
#####
#
# Bland, J.M. and Mutoka, C. and Hutt, M.S.R. ``Kaposi's sarcoma in Tanzania". East African
# Journal of Medical Research, 4, 47--53 (1977).
# Data from a study of the geographical distribution of a tumor, Kaposi's sarcoma, in
# mainland Tanzania.
Region<-c("Tabora","Iringa","Coast","Mara","Ruvuma","Mbeya","Shinyanga","Kigoma",
"Singida","Dodoma","Morogoro","Westlake","Arusha","Mtwara","Mwanza","Tanga",
"Kilimanjara")
Popw10kHF<-c(1.8, 2.6, 4.0, 4.4, 6.6, 6.7, 9.0, 9.2, 10.8, 11.1, 11.7, 12.5, 13.7, 14.8,
20.7, 23.0, 57.3)
CasePMY<-c(2.37, 8.46, 1.28, 4.29, 7.21, 2.06, 1.66, 4.22, 6.17, 2.6, 6.33, 6.6, 2.46,
6.4, 8.54, 4.54, 6.65)
op<-par(mfrow=c(1,1), mgp=c(1.8,.5,0), mar=c(2.8,2.7,2,1),oma=c(0,0,0,0))
plot(Popw10kHF,CasePMY,main="",pch=19,cex=0.9,cex.lab=0.9,cex.axis=0.8,col="navy")
text(Popw10kHF,CasePMY,labels=Region,pos=3,cex=0.7)
abline(v=mean(Popw10kHF),col="black",lty=2,lwd=1)
abline(h=mean(CasePMY),col="darkblue",lty=2,lwd=1)
par(op)
rankor(CasePMY,Popw10kHF,"sp","ga","dubois",FALSE,"greater",TRUE)
#####
#
\dontrun{
data(Berk);attach(Berk)
op<-par(mfrow=c(1,1), mgp=c(1.8,.5,0), mar=c(2.8,2.7,2,1),oma=c(0,0,0,0))
plot(Births,Deaths,main="",pch=19,cex=0.9,cex.lab=0.9,cex.axis=0.8)
abline(h=mean(Deaths),col="black",lty=2,lwd=1)
abline(v=mean(Births),col="black",lty=2,lwd=1)
text(Births[12],Deaths[12],labels="noon",pos=3,cex=0.7)
text(Births[24],Deaths[24],labels="midnight",pos=4,cex=0.7)
W<-Births[-c(12,24)];Z<-Deaths[-c(12,24)]
plot(W,Z,main="",xlab="Births",ylab="Deaths",pch=19,cex=0.99,cex.lab=0.99,
cex.axis=0.8)
text(W,Z,labels=paste("h",1:24,sep=""),cex=0.6, pos=
c(1,3,4,1,1,1,4, 2,1,1,1,1,1,1,1,1,3,1,1,1,2,1,4,2) )
abline(h=mean(Z),col="black",lty=2,lwd=1)
abline(v=mean(W),col="black",lty=2,lwd=1)
par(op)
A<-matrix(NA,10,5);Series<-c("Complete","Clean")
colnames(A)<-c("Data set","Coeff.","Value", "C. Two-tail p","L. Two-tail p")
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
}
#####
#
\dontrun{
data(BlaAlt);attach(BlaAlt)
op<-par(mfrow=c(1,1))
plot(Fev1,Fev2,main="",xlab="First measurement of Fev",ylab="Second measurement of Fev",
pch=19,cex=0.9, col= "gold2" )
abline(v=mean(Fev1),col="black",lty=2,lwd=1)
abline(h=mean(Fev2),col="darkblue",lty=2,lwd=1)
par(op)
rankor(Fev1,Fev2,"sp","ga","woodbury",FALSE,"two",TRUE)
rankor(Fev1,Fev2,"sp","ga","midrank",FALSE,"two",TRUE)
rankor(Fev1,Fev2,"sp","ga","dubois",FALSE,"two",TRUE)
rankor(Fev1,Fev2,"sp","ga","gh",FALSE,"two",TRUE)
rankor(Fev1,Fev2,"sp","ga","wgh",FALSE,"two",TRUE)
detach(BlaAlt)
}
#####
#
data(Locre);attach(Locre)
op<-par(mfrow=c(1,1))
plot(Males,Females,main="Fer cryin' out loud - there is a sex difference",
xlab="Females",ylab="Males",pch=19,cex=0.8,col="steelblue")
text(Males,Females,labels=1:length(Females),cex=0.7,pos=2)
abline(h=mean(Females),col="black",lty=2,lwd=1)
abline(v=mean(Males),col="black",lty=2,lwd=1)
par(op)
out<-rankor(Males,Females,"g","vg","gh",FALSE,"greater")
cat(out$Value,out$Cpv,"\n")
detach(Locre)
#####
#
\dontrun{
# Relationship between the size of caudolateral curvilinear osteophyte# of the canine femoral neck and the radiographic view
data(Femurs);attach(Femurs)
cos<-ifelse(Gender==1,"steelblue","pink4")
txs<-ifelse(Gender==1,"F","M")
op<-par(mfrow=c(1,2))
	plot(Age,CCO,main="Plot_1",
	xlab="Age",ylab="CCO",pch=19,cex=0.7,col=cos)
	text(Age,CCO,labels=txs,cex=0.6,pos=2)
	abline(h=mean(CCO),col="black",lty=2,lwd=1)
	abline(v=mean(Age),col="black",lty=2,lwd=1)
	plot(BW,CCO,main="Plot_2",
	xlab="Body weight",ylab="CCO",pch=19,cex=0.7,col=cos)
	text(BW,CCO,labels=txs,cex=0.6,pos=2)
	abline(h=mean(CCO),col="black",lty=2,lwd=1)
	abline(v=mean(BW),col="black",lty=2,lwd=1)
par(op)
out<-rankor(Age,CCO,"g","st","gh",FALSE,"two")
cat(out$Value,out$Cpv,"\n")
out<-rankor(BW,CCO,"fa","st","gh",FALSE,"two")
detach(Femurs)
}
#####
#
\dontrun{
data(FrigateB);attach(FrigateB)
op<-par(mfrow=c(1,1))
plot(Vol,Frq,pch=19,col="darkgreen")
abline(h=mean(Frq),col="black",lty=2,lwd=1)
abline(v=mean(Vol),col="black",lty=2,lwd=1)
par(op)
Evar<-function(A){
	index<-c("s","k","g","r4","fy","fa");approx<-c("ex","ga","st","vg")
	for (ind in index){for(app in approx){rankor(A[,1],A[,2],ind,app,"wgh",FALSE,
			                   "less",TRUE)}}}
Evar(FrigateB)
detach(FrigateB)
}
#####
#
\dontrun{
All.index<-function(X,type){
# Computes the observed rank correlation statistics and their p-values.
	Index<-c("spearman","kendall","gini","r4","fy","fa")
	n<-nrow(X);A<-matrix(0,6,6)
	colnames(A)<-c(" Observed","  t-Student","  Gaussian","  VGGFR",
	"Exact Cons."," Exact Lib.")
	rownames(A)<-c("Spearman","Kendall ","Gini ","   r4", "Fisher-Yates"," fa")
	for (i in 1:6){
		a<-rankor(X[,1],X[,2],Index[i],"St","woodbury",FALSE,type,FALSE)
		r<-a$Value;A[i,2]<-a$Cpv
		a<-ranktes(r,n,Index[i],"Ga",FALSE,type,FALSE);A[i,3]<-a$Cpv
		a<-ranktes(r,n,Index[i],"Vg",FALSE,type,FALSE);A[i,4]<-a$Cpv
		a<-ranktes(r,n,Index[i],"Ex",FALSE,type,FALSE);A[i,5]<-a$Cpv;A[i,6]<-a$Lpv
		A[i,1]<-r}
return(A)}
# Data for the sample run is from Sokal and Rohlf (Box 15.6, 1981; 
# or Box 15.7, 1995): Computation of rank correlation coefficient between
# the total length of 15 aphid stem mothers and the mean thorax length of 
# their parthenogenetic offspring.
X<-matrix(c(8.7,5.95,8.5,5.65,9.4,6,10,5.7,6.3,4.7,7.8,5.53,11.9,6.4,6.5,4.18,6.6,
6.15,10.6, 5.93,10.2,5.7,7.2,5.68, 8.6,6.13,11.1,6.3, 11.6,6.03), 15, 2, byrow=TRUE) 
A<-All.index(X,"two-sided") # type of alternative
print(round(A,4))
}
}
