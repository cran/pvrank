\name{DealwT}
\alias{DealwT}
\title{Treatment of ties}
\description{Calls a set of FORTRAN subroutines to deal with tied scores.}
\usage{DealwT(n, x, y, indr="spearman", ties = "woodbury", sizer = 1000000,
	repl= 1000, print = TRUE)}
\arguments{
	\item{n}{number of ranks.}
	\item{x}{a numeric vector.}  
	\item{y}{a numeric vector with compatible dimensions to x.}
	\item{indr}{rank order association coefficient. Possible options are: 1="spearman",
		
		2="kendall", 3="gini", 4="r4".}
	\item{ties}{method for breaking ties. Possible options are: 1="woodbury", 2="gh",
		
		3="wgh", 4= "midrank", 5="dubois", 6="no ties".}
	\item{sizer}{number of replications for resolving ties by randomization. Default \eqn{10^6}.}
	\item{repl}{number of sampled pairs of permutations needed to apply the weighted 
		
		max-min method. Default \eqn{10^3}.}
	\item{print}{\code{FALSE} suppresses some of the output.}
	}
\value{
\item{n}{number of ranks}
\item{x}{a numeric vector}  
\item{y}{a numeric vector with the same dimension as x}
\item{indr}{rank correlation coefficient}
\item{ties}{method for breaking ties }
\item{sizer}{number of replications for resolving ties by randomization}
\item{repl}{number of sampled pairs of permutations needed to apply the weighted max-min method}
\item{r}{the observed value of the reuired coefficient} 
\item{Hilo}{a numeric vector containing the bounds obtained with the
	Gideon-Hollister method, if demanded} 
}
\note{The routines are written in FORTRAN 90. When a sequence of more than nine tied scores is present in one or both rankings, the execution is halted.}
\author{Agostino Tarsitano}
