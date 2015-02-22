Timc<- function(n,mu2n,mu4n){
			lam<-vector(mode = "numeric", length = 2)	
			Eval<-vector(mode = "numeric", length = 1)
			y <-.Fortran("CRS",lam,as.double(mu2n),as.double(mu4n),Eval,as.double(n))
			return(y)}