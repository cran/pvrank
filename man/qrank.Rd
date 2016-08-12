\name{qrank}
\alias{qrank}
\title{Quantiles of exact and approximated null distributions of rank correlations}
\description{For a given level of significance, this routine computes approximated or exact conservative 
	and/or liberal critical values under the hypothesis of no association.}
\usage{qrank(prob, n, index="spearman", approx="vggfr", print = FALSE,
	lower.tail = TRUE)}
\arguments{
  \item{prob}{the nominal level of significance.} 
  \item{n}{the number of ranks.}
  \item{index}{a character string that specifies the rank correlation used in the test. Acceptable values are:
  "spearman","kendall","gini","r4","fy","fa". Only enough of the string to be unique is required.}
  \item{approx}{a character string that specifies the type of approximation to the null distribution: "vggfr", "exact","gaussian","student".} 
  \item{print}{\code{FALSE} suppresses partial output.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probability \eqn{P [X\leq x]}{ P[X <= x]} is computed, otherwise \eqn{P[X>x]}{P[X>x]}. In brief, lower tailed tests are used to test for negative correlation and upper tailed tests are used to test for positive correlation.}  
}
\details{
This routine provides two exact quantiles corresponding to a conservative level (next smaller exact size) and a liberal level (next larger exact size). However, in the case of \eqn{n>26} (Spearman) or \eqn{n>60} (Kendall) or \eqn{n>24} (Gini) or \eqn{n>15} (\eqn{r_4}, Fisher-Yates coefficient and Filliben-Amerise rank correlation), an approximated, but unique quantile is produced according to \code{approx}. The default option is \code{"vggfr"} in the case of Spearman and \eqn{r_4}; \code{"gaussian"} for Kendall, "fy", and "fa"; \code{"student"} for Gini's cograduation. 

Exact computations use frequencies obtained by complete enumeration for Spearman's \eqn{r_1}, Gini's \eqn{r_2}, \eqn{r_4}, \eqn{r_5} (Fisher-Yates coefficient), and \eqn{r_6} (Filliben-Amerise rank correlation). A recursive formula is employed in the case of Kendall's \eqn{r_3}.

}
\value{a list containing the following components:
\item{n}{number of ranks.}
\item{Statistic}{coefficient of rank order association}
\item{Level}{nominal level}
\item{Cq}{conservative quantile}
\item{Cv}{conservative p-value}
\item{Lq}{liberal quantile}
\item{Lv}{liberal p-value}
}

\note{ 
The quantile function \eqn{Q(.)} of a symmetrical distribution satisfies:
\deqn{Q(0.5)-Q(p)=Q(1-p)-Q(0.5)\quad for \ 0<p< 0.5}
}


\author{Agostino Tarsitano, Ilaria Lucrezia Amerise and Marco Marozzi}

\examples{ 
\dontrun{
data(Insuhyper); attach(Insuhyper)
op<-par(mfrow=c(1,1), mgp=c(1.8,.5,0), mar=c(2.8,2.7,2,1),oma=c(0,0,0,0))
plot(PI,TG,main="Rank correlation between obesity and triglyceride response",
xlab="Ponderal Index", ylab="Plasma triglyceride concentration", pch=19,
cex=0.9, col= "rosybrown4")
text(PI,TG,labels=rownames(Insuhyper),cex=0.6,pos=c(rep(3,10),1,3,1,rep(3,4),1.3))
abline(v=mean(PI),col="black",lty=2,lwd=1)
abline(h=mean(TG),col="darkblue",lty=2,lwd=1)
par(op)
r<-comprank(PI,TG,"spearman","gh")$r
a1<-qrank(0.025, length(PI), "sp", "vggfr", print = FALSE,lower.tail = TRUE)$Cq
a2<-qrank(0.975, length(PI), "sp", "vggfr", print = FALSE,lower.tail = TRUE)$Cq
cat(round(a1,4),round(r,4),round(a2,4))
r<-comprank(PI,TG,"kendall")$r
b1<-qrank(0.95, length(PI), "ke", "ex", lower.tail = TRUE)$Cq # p = .05, one-tailed (upper)
b2<-qrank(0.05, length(PI), "ke", "ex", lower.tail = TRUE)$Cq # p = .95, one-tailed (upper)
cat(round(b2,4),round(r,4),round(b1,4))
detach(Insuhyper)
}
#####
#
a<-qrank(0.10,61,"K","St")
a<-qrank(0.01,25,"S","Ga",FALSE,FALSE);a$Cq
a<-qrank(0.03,11,"fy","Ga",FALSE,FALSE);a$Cq
a<-qrank(0.001,15,"r4","Ex",FALSE,FALSE);cat(a$Cq,a$Lq,"\n")
a<-qrank(0.01,14,"fil","Ex",FALSE,FALSE);cat(a$Cq,a$Lq,"\n")
#####
#
\dontrun{
a<-qrank(0.05,27,"G","Vg",FALSE,FALSE);a$Cq
}
}