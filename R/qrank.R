qrank<-function(prob,n,index="spearman",approx="vggfr",print=FALSE,lower.tail=TRUE){
	if (!is.numeric(prob)| !is.numeric(n)) {stop("Non-numeric argument to mathematical function")}
	if (!(prob>0 & prob<1)) {stop("The nominal significance level must satisfy 0<prob<1)");print(prob)}
	nn<-trunc(n)
	if(n<0) {stop("A negative number of ranks is indicated")}
	if(!(nn==n)) {warning("A non integer number of ranks is indicated. A truncation is necessary.");n<-nn}
	if (n<5) {stop("The number of ranks must be at least 5");print(n)}
	prob1<-prob;Lowt<-"(lower tail)."
	if (lower.tail==FALSE) {prob1<-1-prob;Lowt<-"(upper tail)."}
	cifer<-9;precis<-53
	nu<-mpfr(n,300);nm1<-nu-1;nm2<-nu-2
	C32<-mpfr(1.00762,300);C33<-mpfr(2.01524,300)
	ain<-c("spearman","kendall","gini","r4")
	index<-tolower(index)
	index<-match.arg(index, ain, several.ok = TRUE)
	apx<-c("gaussian","student","vggfr","exact");approx<-tolower(approx)
	approx<-match.arg(approx, apx, several.ok = TRUE)	
	ksw<-FALSE
	if ((index=="spearman" & n>26) | (index=="kendall" & n>60) | (index=="gini" & n>24) | (index=="r4" & n>15)) {ksw<-TRUE}
	if (approx=="exact" & ksw){
		    cat("The exact p-value is not available for the indicated number of ranks \n  An approximate p-value is computed \n")
			if (index=="r4") {approx<-"student"} else {approx<-"vggfr"}
		}
	if (approx=="vggfr"){eappr<-"VGGFR"
		    a<-ranktes(0.0,n,index,"vggfr",F,"l",F)
			L1<-a$Lambda[1];L2<-a$Lambda[2]
			dVian<-function(x,L1,L2){L1*(1-abs(x)^L1)^L2/(2*beta(1/L1,1+L2))}
			Crit<-function(x,prob,L1,L2){abs(prob-integrate(dVian,L1=L1,L2=L2,-1,x)$value)}
			rc<-optimize(Crit,c(-1,1),prob,L1,L2,tol=0.0000001)$minimum}
	if ("student"==approx){eappr<-"t-Student"
			 if (index=="spearman"){Lamx<-nm2
			 									  Tx<-qt(prob, as.numeric(Lamx), lower.tail = TRUE)
			 	 								   a<-Tx^2/Lamx;rc<-sqrt(a/(1+a))}
			 if (index=="gini"){kn<-n%%2;z1<-3*nm1*(nu^2-kn);z2<-2*(nu^2+2+kn);Lm<-0.5*(z1/z2 -1)
		  	    							  Lamx<-2*trunc(Lm+0.5)
		  	    							  Tx<-qt(prob,as.numeric(Lamx), lower.tail = TRUE)
		  	    							   a<-Tx^2/Lamx;rc<-sqrt(a/(1+a))}
			 if (index=="kendall"){z1<-(4*nu+10)/(9*nu*nm1);Lm<-0.5*(1/z1-1);Lamx<-2*trunc(Lm+0.5)
			 								 Tx<-qt(prob, as.numeric(Lamx), lower.tail = TRUE)
			 								 a<-Tx^2/Lamx;rc<-sqrt(a/(1+a))}
		     if (index=="r4") {Lamx<-trunc((nu-C33)/C32)
		     							     Tx<-qt(prob, as.numeric(Lamx), lower.tail = TRUE)
		     							     a<-Tx^2/Lamx;rc<-sqrt(a/(1+a))} 
		    if (prob<0.5){rc<- -abs(rc)}  					  			 						  			 						  
			 }
	if ("gaussian" == approx){eappr<-"Gaussian"
			 zx<-qnorm(prob, mean = 0, sd = 1)
			 if (index=="spearman"){rc<-zx/sqrt(nm1)}
			 if (index=="gini"){rc<-zx/sqrt(1.5*nu)}
			 if (index=="kendall"){rc<-zx/sqrt((9*nu*nm1)/(2*(2*nu+5)))}
			 if (index=="r4"){rc<-sqrt(C32)*zx/sqrt(nm1)}
			 if (prob<0.5){rc<- -abs(rc)}
			 		}
	if(!(approx=="exact")){
			if (print){
					cat("Statistic: ",index," n:",n, " Nominal significance level:",noquote(sprintf("%.5f",as.numeric(prob1))),Lowt,"\n")
					cat("Approximate quantile:",noquote(sprintf("%.5f",as.numeric(rc))),"\n")
					if (!(approx=="vggfr")){cat("Approx.: ",eappr,"\n")}
					else {cat("Approx. : VGGFR with (",noquote(sprintf("%.6f",as.numeric(L1))),",",noquote(sprintf("%.6f",as.numeric(L2))),")","\n")}
					}
		outrank<-list(n=n, Statistic=index, Level=prob, Cq=as.numeric(rc))
		return(outrank)}
	eappr<-"Exact"
	(opr <- mpfr_default_prec())
	stopifnot(opr == (oprec <- mpfr_default_prec(300)), 300  == mpfr_default_prec())
	mpfr_default_prec(opr)
	C1<-mpfr(25,300);C2<-mpfr(38,300);C3<-mpfr(35,300);C4<-mpfr(72,300)
	C5<-mpfr(5,300);C6<-mpfr(9,300);C7<-mpfr(100,300);C8<-mpfr(328,300);C9<-mpfr(127,300)
	C0<-mpfr(997,300);C11<-mpfr(372,300);C12<-mpfr(1350,300);C13<-mpfr(111,300)
	C14<-mpfr(153,300);C15<-mpfr(366,300);C16<-mpfr(301,300);C17<-mpfr(456,300)
	C18<-mpfr(912,300);C19<-mpfr(1248,300);C20<-mpfr(107,300);C21<-mpfr(4,300)
	C22<-mpfr(76,300);C23<-mpfr(182,300);C24<-mpfr(307,300);C25<-mpfr(315,300)
	C26<-mpfr(342,300);C27<-mpfr(420,300);C28<-mpfr(315,300);C29<-mpfr(105,300)
	nu<-mpfr(n,300);nm1<-nu-1;nm2<-nu-2;nm3<-nu-3;np1=nu+1
	if  ("kendall" == index) {Cs1<-2/(n^2-n);Cm<-(n^2-n)/2+1
		 	B<-mpfr(matrix(0,Cm,2),300)			 				
		 	fname<-paste("Kend",n,".txt",sep="")
		 	fpath <- system.file("extdata", fname, package="pvrank")
		 	A<-read.table(fpath,sep="",colClasses = "character")
		    B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		    B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)
		     j1<-which(B[,2]<=prob)[length(which(B[,2]<=prob))]
		     j2<-which(B[,2]>=prob)[1]		       					
		     Cpv<-formatMpfr(B[j1,1],digits=cifer);Lpv<-formatMpfr(B[j2,1],digits=cifer)
		     Cqv<-formatMpfr(B[j1,2],digits=cifer);Lqv<-formatMpfr(B[j2,2],digits=cifer)
		     }
	if  ("gini" == index) {kn<-n%%2;Cs1<-2/(n^2-kn);Cm<-(n^2-kn)/2+1
		 	B<-mpfr(matrix(0,Cm,2),300)	
		 	fname<-paste("Gini",n,".txt",sep="")
		 	fpath <- system.file("extdata", fname, package="pvrank")
		 	A<-read.table(fpath,sep="",colClasses = "character")
		    B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		    B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)
		     j1<-which(B[,2]<=prob)[length(which(B[,2]<=prob))]
		     j2<-which(B[,2]>=prob)[1]
			Cpv<-formatMpfr(B[j1,1],digits=cifer);Lpv<-formatMpfr(B[j2,1],digits=cifer)
			Cqv<-formatMpfr(B[j1,2],digits=cifer);Lqv<-formatMpfr(B[j2,2],digits=cifer)
			}
	  if ("spearman" == index) {Cs1<-12/(n^3-n);Cs2<-3*np1/nm1;Cm<-(n^3-n)/6+1	
		 	B<-mpfr(matrix(0,Cm,2),300)
		 	fname<-paste("Spear",n,".txt",sep="")
		 	fpath <- system.file("extdata", fname, package="pvrank")
		 	A<-read.table(fpath,sep="",colClasses = "character")
		    B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		    B[,1]<-Cs1*B[,1]-Cs2;B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)
		    j1<-which(B[,2]<=prob)[length(which(B[,2]<=prob))];j2<-which(B[,2]>=prob)[1]
		    Cpv<-formatMpfr(B[j1,1],digits=cifer);Lpv<-formatMpfr(B[j2,1],digits=cifer)
		    Cqv<-formatMpfr(B[j1,2],digits=cifer);Lqv<-formatMpfr(B[j2,2],digits=cifer)
		    }
	  if ("r4" == index) {
		 	fname<-paste("R4",n,".csv",sep="")
		 	fpath <- system.file("extdata", fname, package="pvrank")
			B<-read.csv(fpath,header=F,sep=" ")
			B<-as.matrix(noquote(B))
		    B[,1]<-B[,1]/10000
		    ws1<-sum(as.numeric(B[,2]));B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])
		    j1<-which(B[,2]<=prob)[length(which(B[,2]<=prob))]
		    j2<-which(B[,2]>=prob)[1]
		    Cpv<-B[j1,1];Lpv<-B[j2,1];Cqv<-B[j1,2];Lqv<-B[j2,2]}
		    if(!lower.tail){Cqv<-1-Cqv;Lqv<-1-Lqv}
		    if(print){
		    	cat("Statistic:",index,"n:",n, "Nominal significance level:",noquote(sprintf("%.6f",as.numeric(prob1))),Lowt,"Approx.:",eappr,"\n")
		   		cat("Conservative quantile:",noquote(sprintf("%.6f",as.numeric(Cpv)))," Real level:",noquote(sprintf("%.6f",as.numeric(Cqv))),"\n")
		    	cat("Liberal quantile          :",noquote(sprintf("%.6f",as.numeric(Lpv))), " Real level:",noquote(sprintf("%.6f",as.numeric(Lqv))),"\n")
		    }		    
		    outrank<-list(n=n,Statistic=index,Level=prob,Cq=as.numeric(Cpv),Cv=as.numeric(Cqv),Lq=as.numeric(Lpv),Lv=as.numeric(Lqv))
	return(outrank)}