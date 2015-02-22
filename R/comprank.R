# Compute various rank correlations
comprank<-function (x,y=NULL,index="spearman",tiex="woodbury",sizer=1000000,repl=1000,print=TRUE){
	 if (is.matrix(x)){
	 	     if (!is.numeric(x)) {stop("Non-numeric argument to mathematical function")}
	 	    w<-rep(NA,6);z<-rep(NA,7)
	  		w[1]<-y;w[2]<-index;w[3]<-tiex;w[4]<-sizer;w[5]<-repl;w[6]<-print
	  		for(i in 1:6){z[i+1]<-w[i]}
	  		index<-as.character(z[2]);tiex<-as.character(z[3])
	  		sizer<-noquote(z[5]);repl<-noquote(z[6]);print<-as.logical(noquote(z[7]))
	  	}
	  else {
	  			if (length(x) != length(y)) {stop("x and y must have same length.")}
	  			 if (!is.numeric(x) | !is.numeric(y) ) {stop("Non-numeric argument to mathematical function")}
	  		   }
	  ain<-c("spearman","kendall","gini","r4")
	  tos<-c("woodbury","gh","wgh","midrank","dubois")
	  index<-tolower(index);index<-match.arg(index, ain, several.ok = TRUE)
	  tiex<-tolower(tiex);tiex<-match.arg(tiex, tos, several.ok = TRUE)
	  if (((tiex=="midrank") | (tiex=="dubois")) & (index=="r4")){
	  	cat("Combination:",index,tiex,"\n")	  	
	  	stop("Such a feature has not yet been implemented")
	  	}
	  if (!(is.matrix(x))){
	  		isw<-0
	 		nx<-length(x);ny<-length(y)
	  		n<-max(nx,ny);p<-1:n
	  		if (!(nx==ny)) {stop( "x and y are not of the same length")}
      		x1<-rank(x,ties.method="average");y1<-rank(y,ties.method="average")
	  		z<-p %in% x1;n1<-length(which(!z))
	  		z<-p %in% y1;n2<-length(which(!z))
	  		isw<-n1+n2
      		ties<-match(tiex,tos)
      		if(isw==0) {ties<-6}
      		ties<-as.integer(ties);indr<-match(index,ain)
      		a<-DealwT(n,x1,y1,indr,ties,sizer,repl,print)
     		r<-a$r}
     	if (is.matrix(x)){m<-ncol(x);m1<-m-1;n<-nrow(x);p<-1:n;r<-matrix(1,m,m)
     		for(i in 1:m1){i1<-i+1
     			for(j in i1:m){
     				x1<-rank(x[,i],ties.method="average");y1<-rank(x[,j],ties.method="average")
	  				z<-p %in% x1;n1<-length(which(!z))
	  				z<-p %in% y1;n2<-length(which(!z))
	  				isw<-n1+n2;ties<-match(tiex,tos)
      				if(isw==0) {ties<-6} else {cat("Tied scores are present in one or both rankings",i,j,"\n")}
      				ties<-as.integer(ties);indr<-match(index,ain)
      				a<-DealwT(n,x1,y1,indr,ties,sizer,repl,print)
     				r[i,j]<-a$r; r[j,i]<-a$r
     				}}
     		}
return(r)}