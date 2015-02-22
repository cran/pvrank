DealwT<- function(n,x,y,indr="spearman",ties="woodbury",sizer=1000000,repl=1000,print=TRUE){
	ifault<-0;ifault<-as.integer(ifault);r<-0
	if (((ties=="midrank") | (ties=="dubois")) & (indr=="r4")){
	  	cat("Combination:",indr,ties,"\n")	  	
	  	stop("Such a feature has not yet been implemented")
	  	}
	Hilo<-vector(mode = "numeric", length = 2);Hilo<-rep(0,2)
	y <-.Fortran("DealwT",as.integer(n),as.double(x),as.double(y),as.integer(ties),as.integer(indr),as.double(sizer),as.double(repl),as.double(r),as.double(Hilo),as.integer(ifault))
	names(y) <- c("n", "x", "y","ities","indr","sizer","repl","r","Hilo","ifault")
	if(y$ifault==1) stop("When a sequence of more than 9 tied scores are present in one or both rankings, the execution is halted")	
	if(!(ties==6) & print) {cat("Ties are resolved by using the method:",ties,"\n")}
	if ((ties==2 | ties==3) & print) {cat("Gideon-Hollister bounds",round(y$Hilo,5),"\n")}
	return(y)}