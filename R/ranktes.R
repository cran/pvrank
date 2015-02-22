ranktes<-function(r, n, index="spearman", approx="exact", CC=FALSE, type="two-sided", print=TRUE){
# Computes the p-value of an observed rank correlation
	 if (!is.numeric(r)| !is.numeric(n)) {stop("Non-numeric argument to mathematical function")}
	 if (!is.logical(CC)){type<-CC;CC<-FALSE}
	 if (!(approx=="exact") | index=="r4")	{CC<-FALSE}
	 if (n<5){stop("For n<5 rank correlation tests of independence are not reliable.")}
	 cifer<-9;precis<-53
	 ain<-c("spearman","kendall","gini","r4")
	 apx<-c("gaussian","student","vggfr","exact")
	 alter<-c("two-sided", "less", "greater")
	 index<-tolower(index);approx<-tolower(approx);type=tolower(type)
	 index<-match.arg(index, ain, several.ok = TRUE) 
	 approx<-match.arg(approx, apx, several.ok = TRUE)
	 type<-match.arg(type, alter, several.ok = TRUE)
	 (opr <- mpfr_default_prec())
	 stopifnot(opr == (oprec <- mpfr_default_prec(300)), 300  == mpfr_default_prec())
	 mpfr_default_prec(opr)
	 ksw<-FALSE
	 if ((index=="spearman" & n>26) | (index=="kendall" & n>60) | (index=="gini" & n>24) | (index=="r4" & n>15)) {ksw<-TRUE}
	 if (approx=="exact" & ksw){
		    cat("The exact p-value is not available for the indicated number of ranks \n  An approximate p-value is computed \n")
			if (index=="r4") {approx<-"student"} else {approx<-"vggfr"}
			}
	C1<-mpfr(25,300);C2<-mpfr(38,300);C3<-mpfr(35,300);C4<-mpfr(72,300)
	C5<-mpfr(5,300);C6<-mpfr(9,300);C7<-mpfr(100,300);C8<-mpfr(328,300);C9<-mpfr(127,300)
	C0<-mpfr(997,300);C11<-mpfr(372,300);C12<-mpfr(1350,300);C13<-mpfr(111,300)
	C14<-mpfr(153,300);C15<-mpfr(366,300);C16<-mpfr(301,300);C17<-mpfr(456,300)
	C18<-mpfr(912,300);C19<-mpfr(1248,300);C20<-mpfr(107,300);C21<-mpfr(4,300)
	C22<-mpfr(76,300);C23<-mpfr(182,300);C24<-mpfr(307,300);C25<-mpfr(315,300)
	C26<-mpfr(342,300);C27<-mpfr(420,300);C28<-mpfr(315,300);C29<-mpfr(105,300)
	C32<-mpfr(1.00762,300);C33<-mpfr(2.01524,300)
	nu<-mpfr(n,300);ccf<-mpfr(0,300)
	nm1<-nu-1;nm2<-nu-2;nm3<-nu-3;np1=nu+1
	kn<-n%%2;kkn<-mpfr(kn,300)
	nm1<-nu-1;nm2<-nu-2;nm3<-nu-3;np1<-nu+1;np2<-nu+2;ccf<-0
	if (CC & index=="spearman") {
			S1<-trunc((nu*nu*nu-nu)/6)-(1-r)*(nu*nu*nu-nu)/6
			if (S1<0){sig<-1} else {sig<- -1}
			ccf<- sig*6/(nu*nu*nu-nu)} 
	if (CC & index=="gini") {kkn<-n%%2
			S1<-trunc(0.25*r*(nu*nu-kkn))
			if (S1<0){sig<- 1} else {sig<- -1}
			ccf<- sig*2/(nu^2-kkn)} 
	 if (CC & index=="kendall") {
			S1<-trunc(0.5*r*nu*nm1)
			if (S1<0){sig<- 1} else {sig<- -1}
			ccf<- sig*2/(nu*nm1)}
	 if (CC & index=="r4") {ccf<-0} 
	 rc<-mpfr(r,300)-ccf
	 if (approx=="exact") {rc<- -abs(rc)}
#
		Vian<-function(x,L1,L2){
			La1<-lgamma(1/L1);Lb1<-lgamma(1+L2)
			Lc1<-lgamma(1/L1+1+L2);d1<-exp(La1+Lb1-Lc1)
			L1*(1-abs(x)^L1)^L2/(2*d1)}
#
	if(!(approx=="exact")){
			if (approx=="gaussian"){eappr<-"Gaussian"
					if ("spearman" == index) {zx<-rc*sqrt(nm1);estat<-"Spearman's rho"}
		 			if ("gini" == index) {zx<-rc*sqrt(1.5*nu);estat<-"Gini's gamma"}
			 		if ("kendall" == index) {zx<-rc*sqrt((C6*nu*nm1)/(4*nu+10));estat<-"Kendall's tau"}	
		 			if ("r4" == index) {zx<-rc*sqrt(1.00762*nm1);estat<-"r4"}
					if (type=="two-sided") {Pv<-2*pnorm(-abs(zx), mean = 0, sd=1)}
		 			if (!(type=="two-sided")) {Pv<-pnorm(-abs(zx), mean = 0, sd=1)} 
		 			if (r>0 & type=="less") {Pv<-1-Pv}
					if (r<0 & type=="greater") {Pv<-1-Pv}
					}	   
			if (approx =="student"){eappr<-"t-Student"
					if ("spearman"== index){;estat<-"Spearman's rho"
															zx<-nm2/(1-rc^2);rcp<- abs(rc*sqrt(zx))
														   Lamx<-nm2;Lamx<-asNumeric(Lamx)}
					if ("gini" == index) {estat<-"Gini's gamma"
													kn<-n%%2
				    								z1<-3*nm1*(nu^2-kn);z2<-2*(nu^2+2+kn);Lm<-0.5*(z1/z2 -1)
		  	  	  									Lamx<-2*trunc(Lm);Lamx<-asNumeric(Lamx)
		  	   	 									rcp<- abs(rc*sqrt(2*Lm/(1-rc^2)))}
					if ("kendall" == index){estat<-"Kendall's tau"
													  z1<-(4*nu+10)/(C6*nu*nm1);Lm<-0.5*(1/z1-1)
			 										   Lamx<-2*trunc(Lm);rcp<- abs(rc*sqrt(2*Lm/(1-rc^2)))
			 										   Lamx<-noquote(formatMpfr(Lamx,digits=cifer));Lamx<-asNumeric(Lamx)}
					if ("r4" == index){estat<-"r4"
											   nm2<-trunc((nu-C33)/C32)
		 									   zx<-nm2/(1-rc^2);rcp<- abs(rc*sqrt(zx));Lamx<-nm2;Lamx<-asNumeric(Lamx)}		
		 			rcp<-noquote(formatMpfr(rcp,digits=cifer));rcp<-asNumeric(rcp)
		 			if (type=="two-sided") {Pv<-2*pt(-rcp,Lamx)}
			 		if (!(type=="two-sided")) {Pv<-pt(-abs(rcp),Lamx)}
		 			if (r>0 & type=="less") {Pv<-1-Pv}
					if (r<0 & type=="greater") {Pv<-1-Pv}
					}
			if (approx=="vggfr") {eappr<-"Vianelli GGFR"
					if ("spearman"== index){estat<-"Spearman's rho"
															mu2n<-1/nm1;mu4n=3*(C1*nu^3-C2*nu^2-C3*nu+C4)
			   	 											mu4n<-mu4n/(C1*nu*np1*nm1^3)}
					if ("gini"== index){estat<-"Gini's gamma"
							kn<-n%%2;gf=2/(3*nm1)
							if (kn==0) {ws4<-1+2/(nu^2)
											ws1<-C3*nu^7-C13*nu^6+C14*nu^5-C15*nu^4+C16*nu^3-C17*nu^2-C18*nu+C19
											ws3<-C20*nu^7*nm1*nm3}
					   		else {ws4<-1+C21/(nu^2-1)
									 ws1<-C3*nu^7-C22*nu^6+C23*nu^5-C24*nu^4+C25*nu^3-C26*nu^2-C27*nu-C28
							 		 ws3<-C29*np1^3*nu*nm1^4*nm2}
			     					mu2n<-gf*ws4;mu4n<-C21*ws1/ws3}
 	 				if  ("kendall"== index) {estat<-"Kendall's tau"
 	 													mu2n<-2*(2*nu+C5)/(C6*nu*nm1)
						 								mu4n<-C7*nu^4+C8*nu^3-C9*nu^2-C0*nu-C11
													    mu4n<-mu4n/(C12*(nu*nm1/2)^3)}
					if ("r4"== index){estat<-"r4"
											    mu2n<-C32/nm1
												Dc1<-mpfr(6.96735908,300);Dc2<-mpfr(16.72873450,300);Dc3<-mpfr(35.74780918,300)
												mu4n<-  Dc1/nm1-Dc2/nm1/nm1+Dc3/nm1/nm1/nm1}
			 		a<-Timc(nu,mu2n,mu4n);Lam<-mpfr(a[[1]],300)
		 			if (type=="two-sided") {Pv<-integrateR(Vian,L1=Lam[1],L2=Lam[2],-1,-abs(rc),rel.tol=1e-09)$value;Pv<-Pv*2}
		 			if (!(type=="two-sided")) {Pv<-integrateR(Vian,L1=Lam[1],L2=Lam[2],-1,-abs(rc),rel.tol=1e-09)$value}
		 			if (r>0 & type=="less") {Pv<-1-Pv}
		 			if (r<0 & type=="greater") {Pv<-1-Pv}
		 			Pv<-formatMpfr(Pv,digits=cifer)
		 			La1<-formatMpfr(Lam[1],digits=cifer);La2<-formatMpfr(Lam[2],digits=cifer)
		 			Lz<-as.numeric(c(noquote(La1),noquote(La2)))						 						 
					}	
				if (approx=="gaussian" | approx=="student"){
						if (print){cat("Statistic:",estat,noquote(sprintf("%.4f",r)),"\n","Approximation: ",eappr,"  Alternative :",type,"  Cont. Adj. :",CC,"\n")
				 		 			cat("n=",n,"Observed value:",noquote(sprintf("%.4f",r))," Est. p-value=",noquote(sprintf("%.5f",as.numeric(Pv))),"\n")}
						outrank<-list(n=n,Statistic=estat,Value=r,Approx=eappr,Tails=type, Cpv=as.numeric(Pv),Lpv=as.numeric(Pv))
						}
				if (approx=="vggfr"){
						if(print){cat("Statistic:",estat,noquote(sprintf("%.4f",r)),"\n","Approximation: ",eappr,"  Alternative :",type,"  Cont. Adj. :",CC,"\n")
									cat("n=",n,"Observed value:",noquote(sprintf("%.4f",r))," Est. p-value=",noquote(sprintf("%.5f",as.numeric(Pv))),"\n")}
				outrank<-list(n=n,Statistic=estat,Value=r,Approx=eappr,Tails=type, Cpv=as.numeric(Pv),Lpv=as.numeric(Pv),Lambda=Lz)
					}
	return(outrank)}	 
#
	if (approx=="exact"){eappr<-"Exact level"			
		 if  ("kendall" == index) {estat<-"Kendall's tau"
		 				Cs1<-2/(n^2-n);Cm<-(n^2-n)/2+1
		 				B<-mpfr(matrix(0,Cm,2),300)			 				
		 				fname<-paste("Kend",n,".txt",sep="")
		 				fpath <- system.file("extdata", fname, package="pvrank")
		 				A<-read.table(fpath,sep="",colClasses = "character")
		        		B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		        		B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
		 if  ("gini" == index) {estat<-"Gini's gamma"
		 				kn<-n%%2;Cs1<-2/(n^2-kn);Cm<-(n^2-kn)/2+1
		 				B<-mpfr(matrix(0,Cm,2),300)	
		 				fname<-paste("Gini",n,".txt",sep="")
		 				fpath <- system.file("extdata", fname, package="pvrank")
		 				A<-read.table(fpath,sep="",colClasses = "character")
		        		B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		        		B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
		 if ("spearman" == index & n<=26) {estat<-"Spearman's rho"
		 				Cs1<-(6/(n^3-n));Cs2<-nu*(nu+1)*(2*nu+1)/3;Cm<-(n^3-n)/6+1	
			 			B<-mpfr(matrix(0,Cm,2),300)
		 				fname<-paste("Spear",n,".txt",sep="")
		 				fpath <- system.file("extdata", fname, package="pvrank")
		 				A<-read.table(fpath,sep="",colClasses = "character")
		        		B[,1]<-mpfr(noquote(A$V1),300);B[,2]<-mpfr(noquote(A$V2),300)	
		        		B[,1]<-1-Cs1*(Cs2-2*B[,1])
		        		B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}	       			
			if ("r4" == index) {	estat<-"r4"
						fname<-paste("R4",n,".csv",sep="")
		 				fpath <- system.file("extdata", fname, package="pvrank")
						B<-read.table(fpath,header=F,sep=" ")
						B<-as.matrix(noquote(B));B[,1]<-B[,1]/10000
		    			ws1<-sum(as.numeric(B[,2]));B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])	}
		   	if (!(index=="r4")){	    			
		   					if (type=="two-sided") {j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))];j2<-which(B[,1]>=rc)[1]
		    									ws1<-2*B[j1,2];ws2<-2*B[j2,2]
		    	   								if (ws2>1) {ws2<-1}
		    	   								Cpv<-ws1;Lpv<-ws2}
							if (!(type=="two-sided")){j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))];j2<-which(B[,1]>rc)[1]
		 										  if (r>=0) {Cpv<-formatMpfr(1-B[j2,2],digits=cifer);Lpv<-formatMpfr(1-B[j1,2],digits=cifer)
		 										  				ws1<-formatMpfr(B[j1,2],digits=cifer);ws2<-formatMpfr(B[j2,2],digits=cifer)} 
		 							              else       {Cpv<-formatMpfr(B[j1,2],digits=cifer);Lpv<-formatMpfr(B[j2,2],digits=cifer)
		 							              				ws1<-formatMpfr(1-B[j2,2],digits=cifer);ws2<-formatMpfr(1-B[j1,2],digits=cifer)}
		 				   						  if (type=="greater"){Cpv<-ws2;Lpv<-ws1}
		 				   						  }
		 				   	}
		 	 if (index=="r4"){	    			
		   					if (type=="two-sided") {j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))];j2<-which(B[,1]>=rc)[1]
		    									ws1<-2*B[j1,2];ws2<-2*B[j2,2]
		    	   								if (ws2>1) {ws2<-1}
		    	   								Cpv<-ws1;Lpv<-ws2}
							if (!(type=="two-sided")){j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))];j2<-which(B[,1]>rc)[1]
		 										  if (r>=0) {Cpv<-1-B[j2,2];Lpv<-1-B[j1,2];ws1<-B[j1,2];ws2<-B[j2,2]} 
		 							              else       {Cpv<-B[j1,2];Lpv<-B[j2,2];ws1<-1-B[j2,2];ws2<-1-B[j1,2]}
		 				   						 if (type=="greater"){Cpv<-ws2;Lpv<-ws1}
		 										}
		 						}
 	 	if(print){cat("Statistic:",estat,noquote(sprintf("%.4f",r)),"\n","Approximation: ",eappr,"  Alternative :",type, "Corr. Cont.",CC,"\n")
			        cat(" Conservative p-value:",noquote(sprintf("%.5f",as.numeric(Cpv)))," Liberal p-value:",noquote(sprintf("%.5f",as.numeric(Lpv))),"\n")}
  		outrank<-list(n=n,Statistic=estat,Value=r,approx=eappr,tails=type,Cpv=as.numeric(Cpv),Lpv=as.numeric(Lpv))
  		}
return(outrank)}
